#' Estimate Causal Effects via Targeted Maximum Likelihood Estimation (TMLE)
#'
#' Computes the TMLE estimator at each requested time point. TMLE combines
#' outcome regression (G-comp/ICE backward pass) with a targeting step using
#' the treatment/censoring model predictions (clever covariate). The result is
#' doubly robust: consistent if either the outcome model OR the nuisance models
#' are correctly specified. Inference uses the efficient influence function (EIF).
#'
#' @param obj A \code{longy_data} object with treatment, censoring (if
#'   applicable), and outcome models already fit.
#' @param regime Character. Name of the regime.
#' @param times Numeric vector. Time points at which to estimate. If NULL,
#'   estimates at all time points.
#' @param inference Character. Inference method: \code{"eif"} (default,
#'   efficient influence function), \code{"bootstrap"}, or \code{"none"}.
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param n_boot Integer. Number of bootstrap replicates (only for
#'   \code{"bootstrap"}).
#' @param g_bounds Numeric vector of length 2. Bounds for cumulative g
#'   (denominator of clever covariate). Default \code{c(0.01, 1)}.
#' @param outcome_range Numeric vector of length 2. Range for scaling
#'   continuous outcomes to \code{[0,1]}. If NULL, uses empirical range.
#'   Ignored for binary/survival outcomes.
#' @param risk_set Character. Which subjects form the training set for Q
#'   models at each backward step. \code{"all"} (default) uses all uncensored
#'   subjects. \code{"followers"} restricts to regime-followers, matching
#'   \code{fit_outcome(risk_set = "followers")}.
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class \code{"longy_result"} with elements:
#'   \describe{
#'     \item{estimates}{data.table with time, estimate, se, ci_lower, ci_upper,
#'       n_effective, n_at_risk}
#'     \item{regime}{Name of the regime}
#'     \item{estimator}{\code{"tmle"}}
#'     \item{inference}{Inference method used}
#'     \item{ci_level}{Confidence level}
#'   }
#'
#' @export
estimate_tmle <- function(obj, regime = NULL, times = NULL, inference = "eif",
                          ci_level = 0.95, n_boot = 200L,
                          g_bounds = c(0.01, 1), outcome_range = NULL,
                          risk_set = "all",
                          verbose = TRUE) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)
  inference <- match.arg(inference, c("eif", "bootstrap", "none"))

  for (rname in regime) {

  if (isTRUE(obj$crossfit$enabled)) {
    obj <- .cf_estimate_tmle(obj, regime = rname, times = times,
                              inference = inference, ci_level = ci_level,
                              n_boot = n_boot, g_bounds = g_bounds,
                              outcome_range = outcome_range,
                              risk_set = risk_set,
                              verbose = verbose)
    next
  }

  if (is.null(obj$fits$treatment[[rname]]) || length(obj$fits$treatment[[rname]]) == 0) {
    stop(sprintf("Treatment model not fit for regime '%s'. Run fit_treatment() first.", rname),
         call. = FALSE)
  }
  if (is.null(obj$fits$outcome[[rname]]) || length(obj$fits$outcome[[rname]]) == 0) {
    stop(sprintf("Outcome model not fitted for regime '%s'. Run fit_outcome() first.", rname),
         call. = FALSE)
  }

  nodes <- obj$nodes
  dt <- obj$data
  id_col <- nodes$id
  time_col <- nodes$time
  a_col <- nodes$treatment
  y_col <- nodes$outcome
  all_time_vals <- obj$meta$time_values
  reg <- obj$regimes[[rname]]

  is_binary <- nodes$outcome_type %in% c("binary", "survival")
  is_survival <- nodes$outcome_type == "survival"

  # Read outcome model settings (covariates, learners, bounds, sl_fn)
  outcome_settings <- obj$fits$outcome[[rname]]
  covariates <- outcome_settings$covariates
  learners <- outcome_settings$learners
  q_bounds <- outcome_settings$bounds
  sl_fn <- if (!is.null(outcome_settings$sl_fn)) outcome_settings$sl_fn else "SuperLearner"
  min_obs <- 50L  # match fit_outcome default

  # Determine target times
  if (!is.null(times)) {
    target_times <- sort(all_time_vals[all_time_vals %in% times])
  } else {
    target_times <- all_time_vals
  }

  # Compute cumulative g (shared helper from weights.R)
  g_cum_dt <- .compute_cumulative_g(obj, regime = rname, g_bounds = g_bounds)

  # --- Outcome scaling for continuous outcomes ---
  # Scale Y to [0,1] so we can use quasibinomial for all types
  if (!is_binary) {
    all_y <- dt[[y_col]]
    all_y <- all_y[!is.na(all_y)]
    if (is.null(outcome_range)) {
      y_min <- min(all_y)
      y_max <- max(all_y)
    } else {
      y_min <- outcome_range[1]
      y_max <- outcome_range[2]
    }
    y_range_width <- y_max - y_min
    if (y_range_width == 0) y_range_width <- 1  # degenerate case
  } else {
    y_min <- 0
    y_max <- 1
    y_range_width <- 1
  }

  # Build tracking columns for risk set computation
  dt <- .add_tracking_columns(dt, nodes, reg)
  regime_vals <- .evaluate_regime(reg, dt)
  dt[, .longy_regime_a := regime_vals]

  # Pre-compute lagged regime values for counterfactual treatment history
  k <- if (is.null(nodes$lag_k)) 0 else nodes$lag_k
  max_regime_lags <- if (is.infinite(k)) length(all_time_vals) - 1 else min(k, length(all_time_vals) - 1)
  if (max_regime_lags > 0) {
    for (j in seq_len(max_regime_lags)) {
      lag_regime_col <- paste0(".longy_lag_regime_", a_col, "_", j)
      dt[, (lag_regime_col) := shift(.longy_regime_a, n = j, type = "lag"), by = c(id_col)]
    }
  }

  # For survival: precompute first event time per subject
  if (is_survival) {
    event_rows <- dt[!is.na(dt[[y_col]]) & as.numeric(dt[[y_col]]) == 1]
    if (nrow(event_rows) > 0) {
      first_ev <- event_rows[, list(.longy_first_event = min(get(time_col))),
                              by = c(id_col)]
      dt[first_ev, .longy_first_event := i..longy_first_event, on = id_col]
    }
    if (!".longy_first_event" %in% names(dt)) {
      dt[, .longy_first_event := NA_real_]
    }
  }

  # For competing risks: precompute first competing event time per subject
  has_competing <- !is.null(nodes$competing)
  if (has_competing) {
    d_col <- nodes$competing
    comp_rows <- dt[!is.na(dt[[d_col]]) & as.numeric(dt[[d_col]]) == 1]
    if (nrow(comp_rows) > 0) {
      first_comp <- comp_rows[, list(.longy_first_competing = min(get(time_col))),
                               by = c(id_col)]
      dt[first_comp, .longy_first_competing := i..longy_first_competing,
         on = id_col]
    }
    if (!".longy_first_competing" %in% names(dt)) {
      dt[, .longy_first_competing := NA_real_]
    }
  }

  # Epsilon for bounding Q to avoid log(0)
  eps <- q_bounds[1]

  # Store results per target time
  est_list <- vector("list", length(target_times))

  for (target_idx in seq_along(target_times)) {
    target_t <- target_times[target_idx]

    if (verbose) .vmsg("TMLE target time %d:", target_t)

    # Initialize Q: Y at target time, NA elsewhere
    dt[, .longy_Q := NA_real_]
    dt[dt[[time_col]] == target_t,
       .longy_Q := {
         y_raw <- as.numeric(get(y_col))
         if (!is_binary) (y_raw - y_min) / y_range_width else y_raw
       }]

    # Backward pass: from target_t down to min time
    backward_times <- rev(all_time_vals[all_time_vals <= target_t])
    n_back <- length(backward_times)

    # Store Q* at each time step for EIF computation
    Q_star_list <- vector("list", n_back)
    names(Q_star_list) <- as.character(backward_times)

    for (i in seq_along(backward_times)) {
      tt <- backward_times[i]
      dt_t <- dt[dt[[time_col]] == tt, ]

      # Uncensored subjects through C(t) (Convention B: C before Y).
      still_in <- dt_t$.longy_cum_uncens == 1L
      still_in[is.na(still_in)] <- FALSE

      # Followers-only restriction: limit training to regime-consistent subjects
      if (risk_set == "followers") {
        consist_prev <- dt_t$.longy_consist_prev == 1L
        consist_prev[is.na(consist_prev)] <- FALSE
        still_in <- still_in & consist_prev
      }

      # For survival: exclude absorbed subjects (event strictly before tt)
      if (is_survival) {
        fe <- dt_t$.longy_first_event
        absorbed_primary <- !is.na(fe) & fe < tt
      } else {
        absorbed_primary <- rep(FALSE, nrow(dt_t))
      }
      # For competing risks: exclude subjects with competing event before tt
      if (has_competing) {
        fc <- dt_t$.longy_first_competing
        absorbed_competing <- !is.na(fc) & fc < tt
      } else {
        absorbed_competing <- rep(FALSE, nrow(dt_t))
      }
      at_risk <- still_in & !absorbed_primary & !absorbed_competing
      n_at_risk <- sum(at_risk)
      n_primary <- sum(absorbed_primary)
      n_competing <- sum(absorbed_competing)
      primary_absorbed_ids <- dt_t[[id_col]][absorbed_primary]
      competing_absorbed_ids <- dt_t[[id_col]][absorbed_competing]

      if (n_at_risk == 0 && n_primary == 0 && n_competing == 0) {
        if (verbose) .vmsg("  TMLE time %d: 0 at risk, skipping", tt)
        next
      }

      # Model fitting on at-risk subjects only
      risk_ids <- dt_t[[id_col]][at_risk]
      Q_bar <- numeric(0)
      method <- "none"
      n_train <- 0L

      if (n_at_risk > 0) {
        Q_at_t <- dt_t$.longy_Q[at_risk]
        has_Q <- !is.na(Q_at_t)
        n_train <- sum(has_Q)
      }

      if (n_train > 0) {
        time_index <- match(tt, all_time_vals)
        lag_covs <- .get_lag_covariates(nodes, time_index)
        all_covs <- c(covariates, lag_covs)
        X_risk <- as.data.frame(dt_t[at_risk, all_covs, with = FALSE])
        X_train <- X_risk[has_Q, , drop = FALSE]
        Y_train <- Q_at_t[has_Q]

        # TMLE Q-models always use quasibinomial: predictions feed into
        # qlogis(Q_bar) for fluctuation, so they MUST be in (0,1).
        # gaussian can produce predictions outside [0,1] that clip to
        # bounds, creating extreme logit values and huge epsilon.
        # .safe_sl() handles quasibinomial→binomial swap for SL compatibility.
        step_family <- stats::quasibinomial()

        # Fit Q model
        if (n_train >= min_obs && length(unique(Y_train)) > 1) {
          ctx <- sprintf("TMLE-Q, target=%d, time=%d, n_train=%d", target_t, tt, n_train)
          fit <- .safe_sl(Y = Y_train, X = X_train,
                          family = step_family,
                          learners = learners, cv_folds = 10L,
                          sl_fn = sl_fn, context = ctx, verbose = FALSE)
          method <- fit$method

          # Counterfactual prediction: set A (current + lagged) to regime values
          X_cf <- X_risk
          if (a_col %in% all_covs) {
            X_cf[[a_col]] <- dt_t$.longy_regime_a[at_risk]
          }
          trt_lag_prefix <- paste0(".longy_lag_", a_col, "_")
          for (lc in all_covs[startsWith(all_covs, trt_lag_prefix)]) {
            regime_lc <- sub(trt_lag_prefix,
                             paste0(".longy_lag_regime_", a_col, "_"), lc, fixed = TRUE)
            if (regime_lc %in% names(dt_t)) {
              X_cf[[lc]] <- dt_t[[regime_lc]][at_risk]
            }
          }

          Q_bar <- .predict_from_fit(fit, X_cf)
        } else {
          Q_bar <- rep(mean(Y_train), n_at_risk)
          method <- "marginal"
        }

        # Bound Q_bar
        Q_bar <- .bound(Q_bar, eps, 1 - eps)
      } # end if (n_train > 0)

      # --- TMLE fluctuation (at-risk subjects only) ---
      fluct_epsilon <- 0
      Q_star <- numeric(0)

      if (n_train > 0) {
        # Fluctuation set: regime-consistent AND uncensored through t
        cum_consist <- dt_t$.longy_cum_consist[at_risk] == 1L
        cum_consist[is.na(cum_consist)] <- FALSE
        in_fluct <- cum_consist

        # Get g_cum and g_r for at-risk subjects (Convention B: g_cum includes g_c(t)).
        g_at_t <- merge(
          data.table::data.table(.tmp_id = risk_ids, .tmp_order = seq_along(risk_ids)),
          g_cum_dt[g_cum_dt$.time == tt, c(id_col, ".g_cum", ".g_r"), with = FALSE],
          by.x = ".tmp_id", by.y = id_col, all.x = TRUE
        )
        data.table::setorder(g_at_t, .tmp_order)
        g_cum_vals <- g_at_t$.g_cum
        g_cum_vals[is.na(g_cum_vals)] <- 1

        g_r_vals <- g_at_t$.g_r
        g_r_vals[is.na(g_r_vals)] <- 1

        # g_r only enters at target time (i==1) where actual Y is observed.
        # At intermediate times, pseudo-outcomes are model predictions, not
        # observed data, so the observation mechanism doesn't apply.
        g_denom <- if (i == 1) g_cum_vals * g_r_vals else g_cum_vals

        # Pseudo-outcome: Q values for at-risk subjects
        pseudo_out <- dt_t$.longy_Q[at_risk]

        fluct <- .tmle_fluctuate(Q_bar = Q_bar,
                                 pseudo_outcome = pseudo_out,
                                 g_cum = g_denom,
                                 in_fluctuation_set = in_fluct,
                                 bounds = c(eps, 1 - eps))
        Q_star <- fluct$Q_star
        fluct_epsilon <- fluct$epsilon
      }

      if (verbose) {
        .vmsg("  TMLE time %d: n_at_risk=%d, n_train=%d, n_primary_abs=%d, n_competing_abs=%d, eps=%.4f, method=%s",
              tt, n_at_risk, n_train, n_primary, n_competing, fluct_epsilon, method)
      }

      # --- Predict Q* for newly-censored subjects (needed for EIF) ---
      # Subjects who were at-risk at the previous time but dropped out now
      # (newly censored). They need Q* predictions so the EIF can compute
      # H_{s-1} * (Q*_s - Q*_{s-1}) correctly at the previous step.
      # Without this, Q*_s defaults to 0 in the EIF, causing enormous SEs.
      cens_ids <- dt_t[[id_col]][integer(0)]  # empty, same type as id column
      Q_star_cens <- numeric(0)

      if (n_train > 0) {
        # Detect subjects who were at-risk at t-1 but not at t
        # (Convention B: at-risk = uncensored through C(t), so newly-censored
        # means uncensored through C(t-1) but censored at C(t))
        newly_cens <- !at_risk & !absorbed_primary & !absorbed_competing &
          (dt_t$.longy_uncens_prev == 1L | dt_t$.longy_consist_prev == 1L)
        newly_cens[is.na(newly_cens)] <- FALSE
        # Exclude subjects already at-risk (they have Q* already)
        newly_cens <- newly_cens & !at_risk
        n_newly_cens <- sum(newly_cens)

        if (n_newly_cens > 0) {
          # all_covs already computed above (covariates + lag columns)
          X_cens <- as.data.frame(dt_t[newly_cens, all_covs, with = FALSE])
          X_cens_cf <- X_cens
          if (a_col %in% all_covs) {
            X_cens_cf[[a_col]] <- dt_t$.longy_regime_a[newly_cens]
          }
          trt_lag_prefix <- paste0(".longy_lag_", a_col, "_")
          for (lc in all_covs[startsWith(all_covs, trt_lag_prefix)]) {
            regime_lc <- sub(trt_lag_prefix,
                             paste0(".longy_lag_regime_", a_col, "_"), lc, fixed = TRUE)
            if (regime_lc %in% names(dt_t)) {
              X_cens_cf[[lc]] <- dt_t[[regime_lc]][newly_cens]
            }
          }

          if (method != "marginal" && method != "none") {
            Q_bar_cens <- .predict_from_fit(fit, X_cens_cf)
          } else {
            marg_val <- mean(Q_at_t[has_Q])
            Q_bar_cens <- rep(marg_val, n_newly_cens)
          }

          Q_bar_cens <- .bound(Q_bar_cens, eps, 1 - eps)
          Q_star_cens <- .expit(stats::qlogis(Q_bar_cens) + fluct_epsilon)
          Q_star_cens <- .bound(Q_star_cens, eps, 1 - eps)
          cens_ids <- dt_t[[id_col]][newly_cens]
        }
      }

      # Combine at-risk Q* with hard-coded values for absorbed subjects
      # and predictions for newly-censored subjects
      # Primary absorbed (prior Y=1) → Q*=1; Competing absorbed (prior D=1) → Q*=0
      all_ids <- c(risk_ids, primary_absorbed_ids, competing_absorbed_ids, cens_ids)
      all_Q_star <- c(Q_star, rep(1, n_primary), rep(0, n_competing), Q_star_cens)

      # Store Q* for EIF (includes absorbed and newly-censored subjects)
      if (length(all_ids) > 0) {
        Q_star_list[[as.character(tt)]] <- data.table::data.table(
          .tmp_id = all_ids,
          .Q_star = all_Q_star
        )
        data.table::setnames(Q_star_list[[as.character(tt)]], ".tmp_id", id_col)
      }

      # Propagate: set Q at prev_t = Q*_t (TARGETED, not raw)
      # Includes newly-censored subjects so they get pseudo-outcomes at prev_t
      prev_times <- all_time_vals[all_time_vals < tt]
      if (length(prev_times) > 0 && length(all_ids) > 0) {
        prev_t <- max(prev_times)
        pred_dt_prop <- data.table::data.table(
          .tmp_id = all_ids,
          .tmp_pred = all_Q_star
        )
        data.table::setnames(pred_dt_prop, ".tmp_id", id_col)
        prev_rows <- dt[[time_col]] == prev_t
        prev_ids <- dt[[id_col]][prev_rows]
        match_idx <- match(prev_ids, pred_dt_prop[[id_col]])
        has_match <- !is.na(match_idx)
        dt$.longy_Q[which(prev_rows)[has_match]] <-
          pred_dt_prop$.tmp_pred[match_idx[has_match]]
      }
    }

    # Estimate: psi_hat = mean(Q*_0)
    earliest_t <- backward_times[n_back]
    Q_star_0 <- Q_star_list[[as.character(earliest_t)]]

    if (is.null(Q_star_0) || nrow(Q_star_0) == 0) {
      est_list[[target_idx]] <- data.table::data.table(
        time = target_t, estimate = NA_real_,
        n_effective = 0, n_at_risk = 0L
      )
      next
    }

    # Back-transform if continuous
    psi_scaled <- mean(Q_star_0$.Q_star)
    psi_hat <- if (!is_binary) psi_scaled * y_range_width + y_min else psi_scaled
    n_at_risk <- nrow(Q_star_0)

    est_row <- data.table::data.table(
      time = target_t,
      estimate = psi_hat,
      n_effective = as.numeric(n_at_risk),
      n_at_risk = n_at_risk
    )

    # EIF inference
    if (inference == "eif") {
      D_i <- .compute_tmle_eif(
        Q_star_list = Q_star_list,
        g_cum_dt = g_cum_dt,
        obj = obj,
        target_t = target_t,
        psi_hat = psi_scaled,  # on scaled [0,1]
        regime = rname,
        y_min = y_min,
        y_range_width = y_range_width
      )

      # Preserve decomposition before back-transform (multiplication drops attrs)
      init_raw <- attr(D_i, "initial")
      aug_raw <- attr(D_i, "augmentation")

      # Back-transform EIF if continuous
      if (!is_binary) {
        D_i <- D_i * y_range_width
        if (!is.null(init_raw)) init_raw <- init_raw * y_range_width
        if (!is.null(aug_raw)) aug_raw <- aug_raw * y_range_width
      }

      n_i <- length(D_i)
      if (n_i >= 2 && stats::var(D_i) > 0) {
        se <- sqrt(stats::var(D_i) / n_i)
      } else {
        se <- NA_real_
      }

      # EIF variance decomposition diagnostics
      if (verbose) {
        if (!is.null(init_raw) && !is.null(aug_raw)) {
          se_init <- sqrt(stats::var(init_raw) / n_i)
          se_aug <- sqrt(stats::var(aug_raw) / n_i)
          n_aug_nonzero <- sum(abs(aug_raw) > 1e-10)
          .vmsg("  EIF diagnostics: SE(D)=%.4f, SE(Q*0-psi)=%.4f, SE(aug)=%.4f, n_aug_nonzero=%d/%d",
                se, se_init, se_aug, n_aug_nonzero, n_i)
        }
      }

      z <- stats::qnorm(1 - (1 - ci_level) / 2)
      est_row[, se := se]
      est_row[, ci_lower := if (is.na(se)) NA_real_ else psi_hat - z * se]
      est_row[, ci_upper := if (is.na(se)) NA_real_ else psi_hat + z * se]
    }

    est_list[[target_idx]] <- est_row
  }

  estimates <- data.table::rbindlist(est_list, fill = TRUE)

  # Bootstrap inference
  if (inference == "bootstrap" && n_boot > 0) {
    inf_dt <- .bootstrap_tmle_inference(
      obj, regime = rname, times = target_times,
      n_boot = n_boot, ci_level = ci_level,
      g_bounds = g_bounds, outcome_range = outcome_range,
      verbose = verbose
    )
    # Remove any existing inference columns before adding bootstrap results
    for (col in c("se", "ci_lower", "ci_upper")) {
      if (col %in% names(estimates)) estimates[, (col) := NULL]
    }
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival outcomes
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    if ("se" %in% names(estimates)) {
      z <- stats::qnorm(1 - (1 - ci_level) / 2)
      estimates$ci_lower <- estimates$estimate - z * estimates$se
      estimates$ci_upper <- estimates$estimate + z * estimates$se
    }
  }

  # Clean up
  dt[, .longy_Q := NULL]
  dt[, .longy_regime_a := NULL]
  regime_lag_cols <- grep("^\\.longy_lag_regime_", names(dt), value = TRUE)
  if (length(regime_lag_cols) > 0) dt[, (regime_lag_cols) := NULL]
  if (is_survival && ".longy_first_event" %in% names(dt)) {
    dt[, .longy_first_event := NULL]
  }
  if (has_competing && ".longy_first_competing" %in% names(dt)) {
    dt[, .longy_first_competing := NULL]
  }
  .remove_tracking_columns(obj$data)

  result <- list(
    estimates = estimates,
    regime = rname,
    estimator = "tmle",
    inference = inference,
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_tmle")]] <- result

  } # end for (rname in regime)

  obj
}

#' TMLE Fluctuation Step
#'
#' Fits a quasibinomial GLM with intercept-only, offset = logit(Q_bar),
#' weights = 1/g_cum on the fluctuation set. Returns updated Q* for ALL
#' subjects in the risk set.
#'
#' @param Q_bar Numeric vector. Initial Q predictions for all in risk set.
#' @param pseudo_outcome Numeric vector. Pseudo-outcomes (Q from s+1 or Y_T)
#'   for all in risk set. May contain NAs for subjects outside training set.
#' @param g_cum Numeric vector. Cumulative g values for all in risk set.
#' @param in_fluctuation_set Logical vector. TRUE for subjects in the
#'   fluctuation set (regime-consistent AND uncensored).
#' @param bounds Numeric(2). Bounds for Q_star.
#'
#' @return List with Q_star (numeric vector, same length as Q_bar) and
#'   epsilon (the fluctuation coefficient).
#' @noRd
.tmle_fluctuate <- function(Q_bar, pseudo_outcome, g_cum,
                            in_fluctuation_set, bounds) {
  n <- length(Q_bar)

  # Subjects in fluctuation set with non-missing pseudo-outcome
  in_set <- in_fluctuation_set & !is.na(pseudo_outcome)
  n_fluct <- sum(in_set)
  n_risk <- length(Q_bar)

  if (n_fluct < 2) {
    if (n_fluct > 0) {
      warning(sprintf(
        "Fluctuation set has only %d regime-consistent subject(s). Skipping targeting.",
        n_fluct), call. = FALSE)
    }
    return(list(Q_star = Q_bar, epsilon = 0))
  }

  if (n_risk > 0 && n_fluct < n_risk * 0.1) {
    warning(sprintf(
      "Fluctuation set (%d regime-consistent) is <10%% of risk set (%d). Targeting will have limited effect.",
      n_fluct, n_risk), call. = FALSE)
  }

  offset_vals <- stats::qlogis(Q_bar[in_set])
  y_vals <- pseudo_outcome[in_set]
  w_vals <- 1 / g_cum[in_set]

  # Bound weights to prevent Inf
  w_vals[!is.finite(w_vals)] <- max(w_vals[is.finite(w_vals)], na.rm = TRUE)

  epsilon <- tryCatch({
    fluct_fit <- stats::glm.fit(
      x = matrix(1, nrow = sum(in_set), ncol = 1),
      y = y_vals,
      offset = offset_vals,
      weights = w_vals,
      family = stats::quasibinomial()
    )
    fluct_fit$coefficients[1]
  }, error = function(e) {
    warning(sprintf(
      "TMLE fluctuation failed: %s. Setting epsilon=0 (no targeting).",
      e$message), call. = FALSE)
    0
  })

  if (!is.finite(epsilon)) {
    warning("TMLE fluctuation produced non-finite epsilon; setting epsilon=0.",
            call. = FALSE)
    epsilon <- 0
  }

  # Apply epsilon to ALL risk set subjects (counterfactual predictions)
  Q_star <- .expit(stats::qlogis(Q_bar) + epsilon)
  Q_star <- .bound(Q_star, bounds[1], bounds[2])

  list(Q_star = Q_star, epsilon = epsilon)
}

#' Compute TMLE Efficient Influence Function
#'
#' Computes the doubly-robust EIF for each subject:
#' D_i = Q*_{0,i} - psi + sum_{s=0}^{T} H_{s,i} * (Q*_{s+1,i} - Q*_{s,i})
#'
#' @param Q_star_list Named list keyed by time (as character), each element
#'   a data.table with (id_col, .Q_star).
#' @param g_cum_dt data.table from .compute_cumulative_g() with
#'   (id_col, .time, .g_cum).
#' @param obj longy_data object.
#' @param target_t Target time for this EIF computation.
#' @param psi_hat Point estimate (on scaled \code{[0,1]} for continuous).
#' @param regime Character regime name.
#' @param y_min Numeric. Minimum of outcome range (for scaling continuous Y).
#' @param y_range_width Numeric. Width of outcome range.
#'
#' @return Numeric vector of D_i values, one per subject (on scaled
#'   \code{[0,1]} for continuous).
#' @noRd
.compute_tmle_eif <- function(Q_star_list, g_cum_dt, obj, target_t,
                               psi_hat, regime, y_min = 0,
                               y_range_width = 1) {
  nodes <- obj$nodes
  dt <- obj$data
  id_col <- nodes$id
  time_col <- nodes$time
  reg <- obj$regimes[[regime]]
  all_time_vals <- obj$meta$time_values

  backward_times <- rev(all_time_vals[all_time_vals <= target_t])
  all_ids <- unique(dt[[id_col]])
  n <- length(all_ids)

  # Get Q*_0 for each subject
  earliest_t <- backward_times[length(backward_times)]
  Q0_dt <- Q_star_list[[as.character(earliest_t)]]
  if (is.null(Q0_dt)) {
    return(rep(0, n))
  }

  # Build subject-level data: Q*_0
  subj_dt <- data.table::data.table(.tmp_id = all_ids)
  data.table::setnames(subj_dt, ".tmp_id", id_col)

  Q0_merge <- data.table::copy(Q0_dt)
  data.table::setnames(Q0_merge, ".Q_star", ".Q_star_0")
  subj_dt <- merge(subj_dt, Q0_merge, by = id_col, all.x = TRUE)
  # Subjects not at risk at t=0 get Q*_0 = psi_hat (so their initial term is 0)
  subj_dt[is.na(.Q_star_0), .Q_star_0 := psi_hat]

  # Build augmentation term: sum over time steps
  augmentation <- rep(0, n)
  names(augmentation) <- as.character(all_ids)

  # Need tracking columns for regime-consistency and uncensored
  dt_track <- data.table::copy(dt)
  dt_track <- .add_tracking_columns(dt_track, nodes, reg)

  for (s_idx in seq_along(backward_times)) {
    s <- backward_times[s_idx]

    # Q*_s
    Qs_dt <- Q_star_list[[as.character(s)]]
    if (is.null(Qs_dt)) next

    # Q*_{s+1} (the pseudo-outcome used at step s)
    # For s = target_t, Q*_{s+1} = Y_T (scaled)
    if (s == target_t) {
      # Y at target_t, scaled
      dt_target <- dt[dt[[time_col]] == target_t, ]
      Qs_next <- data.table::data.table(
        .tmp_id = dt_target[[id_col]],
        .Q_star_next = as.numeric(dt_target[[nodes$outcome]])
      )
      data.table::setnames(Qs_next, ".tmp_id", id_col)
      # For continuous, we'd need to scale Y -- but this is called with
      # psi_hat already on scaled [0,1], so we need Y scaled too
      # The caller scales Y in the main function; here we read raw Y
      # and scale it consistent with estimate_tmle's scaling
      # Actually, the Q_star values are already on scaled [0,1] from
      # estimate_tmle, and psi_hat is on that scale too. Y_T was scaled
      # when initialized in estimate_tmle. So we need the scaled Y_T here.
      # Since we don't have access to the scaling params directly, the
      # caller must ensure consistency. For binary, Y is already 0/1.
      # For continuous, we need to know y_min, y_range_width. These are
      # NOT passed to this function. Let's compute from the Q_star values
      # at target_t, which are already scaled predictions.
      # Actually, for the EIF the Q_star_next at target_t should be the
      # observed (scaled) Y_T. But we only have the raw Y in dt.
      # The simplest fix: at target_t, Q_star_next for a subject is their
      # Q_star value (since at target_t, Q was initialized to Y and then
      # the model was fit and fluctuated ON that Y). So Q*_T is the
      # fluctuated version, and Q*_{T+1} = Y_T (the raw pseudo-outcome).
      # We handle this by noting that the augmentation at T is:
      # H_T * (Y_T - Q*_T)
      # For now, we use the Q_star at target_t as Q*_s, and the pseudo-outcome
      # (Y_T) is what was passed to the fluctuation. Let me just use Q_star
      # values directly.
    } else {
      # Q*_{s+1} is the Q_star from the previous backward step
      prev_back_idx <- s_idx - 1
      if (prev_back_idx >= 1) {
        s_next <- backward_times[prev_back_idx]
        Qs_next_raw <- Q_star_list[[as.character(s_next)]]
        if (is.null(Qs_next_raw)) next
        Qs_next <- data.table::copy(Qs_next_raw)
        data.table::setnames(Qs_next, ".Q_star", ".Q_star_next")
      } else {
        next
      }
    }

    # For s = target_t: augmentation_s = H_s * (Y_T_scaled - Q*_T)
    if (s == target_t) {
      # Build the augmentation for this step
      dt_s <- dt_track[dt_track[[time_col]] == s, ]
      # regime-consistent AND uncensored through C(s) (Convention B)
      cum_consist <- dt_s$.longy_cum_consist == 1L
      cum_consist[is.na(cum_consist)] <- FALSE
      cum_uncens <- dt_s$.longy_cum_uncens == 1L
      cum_uncens[is.na(cum_uncens)] <- FALSE
      indicator <- cum_consist & cum_uncens

      # Merge g_cum and g_r (Convention B: g_cum includes g_c(s))
      g_s <- g_cum_dt[g_cum_dt$.time == s, c(id_col, ".g_cum", ".g_r"), with = FALSE]

      # Scale Y_T to [0,1] for continuous (binary already 0/1)
      y_raw <- as.numeric(dt_s[[nodes$outcome]])
      y_scaled <- (y_raw - y_min) / y_range_width

      aug_dt <- data.table::data.table(
        .tmp_id = dt_s[[id_col]],
        .indicator = indicator,
        .Y_T = y_scaled
      )
      data.table::setnames(aug_dt, ".tmp_id", id_col)
      aug_dt <- merge(aug_dt, g_s, by = id_col, all.x = TRUE)
      aug_dt <- merge(aug_dt, Qs_dt, by = id_col, all.x = TRUE)

      aug_dt[is.na(.g_cum), .g_cum := 1]
      aug_dt[is.na(.g_r), .g_r := 1]
      # Exclude subjects missing Q* from augmentation (no prediction available)
      aug_dt[is.na(.Q_star), .indicator := FALSE]
      aug_dt[is.na(.Q_star), .Q_star := 0]  # avoid NaN in arithmetic

      # Note: .Y_T may be NA for unobserved subjects. Their indicator
      # should be FALSE (uncensored doesn't mean observed). But Y_T is the
      # scaled outcome. For EIF, subjects with NA Y contribute 0.
      aug_dt[is.na(.Y_T), .indicator := FALSE]

      aug_dt[, .H := ifelse(.indicator, 1 / (.g_cum * .g_r), 0)]
      aug_dt[, .aug := .H * (.Y_T - .Q_star)]

      # Add to augmentation vector
      match_idx <- match(aug_dt[[id_col]], all_ids)
      valid <- !is.na(match_idx) & is.finite(aug_dt$.aug)
      augmentation[match_idx[valid]] <- augmentation[match_idx[valid]] +
        aug_dt$.aug[valid]
    } else {
      # For s < target_t: augmentation_s = H_s * (Q*_{s+1} - Q*_s)
      dt_s <- dt_track[dt_track[[time_col]] == s, ]
      cum_consist <- dt_s$.longy_cum_consist == 1L
      cum_consist[is.na(cum_consist)] <- FALSE
      # Convention B: uncensored through C(s)
      cum_uncens <- dt_s$.longy_cum_uncens == 1L
      cum_uncens[is.na(cum_uncens)] <- FALSE
      indicator <- cum_consist & cum_uncens

      # Convention B: g_cum includes g_c(s)
      g_s <- g_cum_dt[g_cum_dt$.time == s, c(id_col, ".g_cum", ".g_r"), with = FALSE]

      aug_dt <- data.table::data.table(
        .tmp_id = dt_s[[id_col]],
        .indicator = indicator
      )
      data.table::setnames(aug_dt, ".tmp_id", id_col)
      aug_dt <- merge(aug_dt, g_s, by = id_col, all.x = TRUE)
      aug_dt <- merge(aug_dt, Qs_dt, by = id_col, all.x = TRUE)
      aug_dt <- merge(aug_dt, Qs_next, by = id_col, all.x = TRUE)

      aug_dt[is.na(.g_cum), .g_cum := 1]
      aug_dt[is.na(.g_r), .g_r := 1]
      # Exclude subjects missing Q* from augmentation — without this,
      # NA defaults to 0 creating huge spurious H * (0 - Q*) terms
      aug_dt[is.na(.Q_star), .indicator := FALSE]
      aug_dt[is.na(.Q_star_next), .indicator := FALSE]
      aug_dt[is.na(.Q_star), .Q_star := 0]          # avoid NaN in arithmetic
      aug_dt[is.na(.Q_star_next), .Q_star_next := 0]

      # At intermediate times (s < target_t), g_r does NOT enter the EIF.
      # Pseudo-outcomes are model predictions, not observed data.
      aug_dt[, .H := ifelse(.indicator, 1 / .g_cum, 0)]
      aug_dt[, .aug := .H * (.Q_star_next - .Q_star)]

      match_idx <- match(aug_dt[[id_col]], all_ids)
      valid <- !is.na(match_idx) & is.finite(aug_dt$.aug)
      augmentation[match_idx[valid]] <- augmentation[match_idx[valid]] +
        aug_dt$.aug[valid]
    }
  }

  # D_i = (Q*_0 - psi) + augmentation
  Q_star_0_vals <- subj_dt$.Q_star_0
  initial_term <- Q_star_0_vals - psi_hat
  D_i <- initial_term + augmentation

  .remove_tracking_columns(dt_track)

  # Attach decomposition as attributes for diagnostics
  attr(D_i, "initial") <- initial_term
  attr(D_i, "augmentation") <- augmentation

  D_i
}
