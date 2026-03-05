#' Fit Outcome Regression Models via Iterated Conditional Expectations (ICE)
#'
#' For each target time T, runs a backward sequential regression from T to the
#' earliest time. At each step s, a model is fit for the pseudo-outcome Q
#' (initially Y_T, then predictions from s+1) on uncensored subjects through s.
#' Counterfactual predictions are obtained by setting A to the regime value.
#' The final estimate is mean(Q_hat_0) for each target T.
#'
#' @param obj A \code{longy_data} object with at least one regime defined.
#' @param regime Character. Name of the regime to fit for.
#' @param covariates Character vector. Column names to use as predictors.
#'   If NULL, uses all baseline + timevarying covariates.
#' @param learners Character vector. SuperLearner library names.
#'   If NULL, uses glm only.
#' @param sl_control List. Additional arguments passed to SuperLearner.
#' @param adaptive_cv Logical. Use adaptive CV fold selection based on effective
#'   sample size.
#' @param min_obs Integer. Minimum observations to fit a model; below this
#'   threshold, uses marginal mean.
#' @param bounds Numeric vector of length 2. Bounds for predicted probabilities
#'   (only applied for binary/survival outcomes).
#' @param times Numeric vector. If provided, only fit models for these target
#'   times. Saves computation when estimation is only needed at specific
#'   time points.
#' @param use_ffSL Logical. If TRUE, use future-factorial SuperLearner.
#'   Default FALSE. Forced to FALSE inside parallel workers.
#' @param parallel Logical. If TRUE and a non-sequential \code{future::plan()}
#'   is active, dispatches regime x target-time backward passes in parallel.
#'   Each worker copies the full dataset. Default FALSE.
#' @param risk_set Character. Which subjects form the training set for outcome
#'   models at each backward step. \code{"all"} (default) uses all uncensored
#'   subjects. \code{"followers"} restricts to regime-followers (subjects
#'   consistent with the regime through the previous time), so the model learns
#'   \code{E(Q|H)} without extrapolation to non-followers. Predictions are still made
#'   for all subjects at the earliest time.
#' @param verbose Logical. Print progress.
#' @param refit Logical. If FALSE (default), errors when outcome is already
#'   fitted for the requested regime(s). Set to TRUE to re-fit.
#'
#' @return Modified \code{longy_data} object with outcome fits stored in
#'   \code{obj$fits$outcome}.
#' @export
fit_outcome <- function(obj, regime = NULL, covariates = NULL, learners = NULL,
                        sl_control = list(), adaptive_cv = TRUE,
                        min_obs = 50L, bounds = c(0.005, 0.995),
                        times = NULL, use_ffSL = FALSE,
                        parallel = FALSE,
                        risk_set = c("all", "followers"),
                        verbose = TRUE, refit = FALSE) {
  obj <- .as_longy_data(obj)
  learners <- .resolve_learners(learners, "outcome")
  regime <- .resolve_regimes(obj, regime)
  risk_set <- match.arg(risk_set)

  if (!refit) {
    fitted <- Filter(function(r) !is.null(obj$fits$outcome[[r]]) &&
                       length(obj$fits$outcome[[r]]) > 0, regime)
    if (length(fitted) > 0)
      stop(sprintf("Outcome already fitted for: %s. Use refit=TRUE to override.",
                   paste(fitted, collapse = ", ")), call. = FALSE)
  }

  dt <- obj$data
  nodes <- obj$nodes
  all_time_vals <- obj$meta$time_values
  id_col <- nodes$id
  time_col <- nodes$time
  a_col <- nodes$treatment
  y_col <- nodes$outcome

  # Determine outcome type flags
  is_binary <- nodes$outcome_type %in% c("binary", "survival")
  is_survival <- nodes$outcome_type == "survival"
  has_competing <- !is.null(nodes$competing)

  # Determine target times
  if (!is.null(times)) {
    target_times <- sort(all_time_vals[all_time_vals %in% times])
  } else {
    target_times <- all_time_vals
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying, nodes$treatment)
  }

  # Lag depth for regime counterfactual columns
  k <- if (is.null(nodes$lag_k)) 0 else nodes$lag_k
  max_regime_lags <- if (is.infinite(k)) length(all_time_vals) - 1 else min(k, length(all_time_vals) - 1)

  # Collect regime objects for task function
  regimes_list <- lapply(regime, function(rname) obj$regimes[[rname]])
  names(regimes_list) <- regime

  # Worker SL flag: force sequential SL inside parallel workers
  worker_ffSL <- if (parallel) FALSE else use_ffSL
  worker_verbose <- if (parallel) FALSE else verbose

  # Snapshot data for parallel safety (by-reference semantics)
  if (parallel) dt <- data.table::copy(dt)

  # Flatten regime x target_time into a single task list for better load balancing
  tasks <- expand.grid(regime_idx = seq_along(regime),
                       target_idx = seq_along(target_times),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  one_task <- function(task_row) {
    ridx <- tasks$regime_idx[task_row]
    tidx <- tasks$target_idx[task_row]
    rname <- regime[ridx]
    reg <- regimes_list[[rname]]
    target_t <- target_times[tidx]

    # Each worker needs its own copy because the backward pass mutates .longy_Q
    dt_w <- data.table::copy(dt)

    # Regime-specific prep: tracking columns + regime values
    .add_tracking_columns(dt_w, nodes, reg)
    regime_vals <- .evaluate_regime(reg, dt_w)
    dt_w[, .longy_regime_a := regime_vals]

    # Lagged regime values for counterfactual treatment history
    if (max_regime_lags > 0) {
      for (j in seq_len(max_regime_lags)) {
        lag_regime_col <- paste0(".longy_lag_regime_", a_col, "_", j)
        dt_w[, (lag_regime_col) := shift(.longy_regime_a, n = j, type = "lag"), by = c(id_col)]
      }
    }

    # Survival: precompute first observed event time per subject
    if (is_survival) {
      event_rows <- dt_w[!is.na(dt_w[[y_col]]) & as.numeric(dt_w[[y_col]]) == 1]
      if (nrow(event_rows) > 0) {
        first_ev <- event_rows[, list(.longy_first_event = min(get(time_col))),
                                by = c(id_col)]
        dt_w[first_ev, .longy_first_event := i..longy_first_event, on = id_col]
      }
      if (!".longy_first_event" %in% names(dt_w)) {
        dt_w[, .longy_first_event := NA_real_]
      }
    }

    # Competing risks: precompute first competing event time per subject
    if (has_competing) {
      d_col <- nodes$competing
      comp_rows <- dt_w[!is.na(dt_w[[d_col]]) & as.numeric(dt_w[[d_col]]) == 1]
      if (nrow(comp_rows) > 0) {
        first_comp <- comp_rows[, list(.longy_first_competing = min(get(time_col))),
                                 by = c(id_col)]
        dt_w[first_comp, .longy_first_competing := i..longy_first_competing,
             on = id_col]
      }
      if (!".longy_first_competing" %in% names(dt_w)) {
        dt_w[, .longy_first_competing := NA_real_]
      }
    }

    # Initialize Q: Y only at target time, NA elsewhere
    dt_w[, .longy_Q := NA_real_]
    dt_w[dt_w[[time_col]] == target_t, .longy_Q := as.numeric(get(y_col))]

    # Backward pass: from target_t down to min time
    backward_times <- rev(all_time_vals[all_time_vals <= target_t])
    n_back <- length(backward_times)
    sl_info_target <- vector("list", n_back)

    for (i in seq_along(backward_times)) {
      tt <- backward_times[i]

      dt_t <- dt_w[dt_w[[time_col]] == tt, ]

      # Uncensored subjects through t
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
        absorbed_primary <- still_in & !is.na(fe) & fe < tt
      } else {
        absorbed_primary <- rep(FALSE, nrow(dt_t))
      }
      # For competing risks: exclude subjects with competing event before tt
      if (has_competing) {
        fc <- dt_t$.longy_first_competing
        absorbed_competing <- still_in & !is.na(fc) & fc < tt
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
        if (worker_verbose) .vmsg("  Q target=%d time %d: 0 at risk, skipping",
                                   target_t, tt)
        next
      }

      # Model fitting on at-risk subjects only
      risk_ids <- dt_t[[id_col]][at_risk]
      preds <- numeric(0)
      method <- "none"
      sl_risk <- NULL
      sl_coef <- NULL
      n_train <- 0L
      n_clipped_lower <- 0L
      n_clipped_upper <- 0L

      # Per-step family: quasibinomial (logit link) at all steps for binary/survival
      step_family <- if (is_binary) stats::quasibinomial() else stats::gaussian()

      if (n_at_risk > 0) {
        Q_at_t <- dt_t$.longy_Q[at_risk]
        has_Q <- !is.na(Q_at_t)
        n_train <- sum(has_Q)
      }

      if (n_train > 0) {
        time_index <- match(tt, all_time_vals)
        lag_covs <- .get_lag_covariates(nodes, time_index)
        all_covs <- c(covariates, lag_covs)

        if (worker_verbose) {
          q_base <- nodes$baseline
          q_tv <- nodes$timevarying
          q_extra <- setdiff(covariates, c(q_base, q_tv))
          tv_label <- if (length(q_extra) > 0) {
            c(q_tv, q_extra)
          } else {
            q_tv
          }
          rs_label <- if (risk_set == "followers") " (followers)" else ""
          .vmsg("  Q target=%d time=%d%s: n_train=%d, family=%s",
                target_t, tt, rs_label, n_train, step_family$family)
          .vmsg_covariates(q_base, tv_label, lag_covs)
        }

        X_risk <- as.data.frame(dt_t[at_risk, all_covs, with = FALSE])
        X_train <- X_risk[has_Q, , drop = FALSE]
        Y_train <- Q_at_t[has_Q]

        ow <- NULL
        if (!is.null(nodes$sampling_weights)) {
          ow_risk <- dt_t[[nodes$sampling_weights]][at_risk]
          ow <- ow_risk[has_Q]
        }

        if (n_train >= min_obs && length(unique(Y_train)) > 1) {
        cv_folds <- 10L
        if (adaptive_cv) {
          cv_info <- .adaptive_cv_folds(Y_train, binary = (i == 1 && is_binary))
          cv_folds <- cv_info$V
        }
        ctx <- sprintf("Q(%s), target=%d, time=%d, n_train=%d, mean_Q=%.3f",
                       rname, target_t, tt, n_train, mean(Y_train))
        fit <- .safe_sl(Y = Y_train, X = X_train, family = step_family,
                        learners = learners, cv_folds = cv_folds,
                        obs_weights = ow, use_ffSL = worker_ffSL,
                        context = ctx, verbose = worker_verbose)
        method <- fit$method
        sl_risk <- fit$sl_risk
        sl_coef <- fit$sl_coef

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

        preds <- .predict_from_fit(fit, X_cf)
        } else {
          if (!is.null(ow)) {
            marg <- stats::weighted.mean(Y_train, ow)
          } else {
            marg <- mean(Y_train)
          }
          preds <- rep(marg, n_at_risk)
          method <- "marginal"
        }

        if (is_binary) {
          n_clipped_lower <- sum(preds < bounds[1])
          n_clipped_upper <- sum(preds > bounds[2])
          preds <- .bound(preds, bounds[1], bounds[2])
        } else {
          n_clipped_lower <- 0L
          n_clipped_upper <- 0L
        }
      } # end if (n_train > 0)

      all_ids <- c(risk_ids, primary_absorbed_ids, competing_absorbed_ids)
      all_preds <- c(preds, rep(1, n_primary), rep(0, n_competing))

      if (worker_verbose) {
        .vmsg("  Q(%s) target=%d time %d: n_at_risk=%d, n_train=%d, n_primary_abs=%d, n_competing_abs=%d, method=%s, clipped=%d/%d (lo/hi)",
              rname, target_t, tt, n_at_risk, n_train, n_primary, n_competing, method,
              n_clipped_lower, n_clipped_upper)
      }

      sl_info_target[[i]] <- list(time = tt, target_time = target_t,
                                  method = method,
                                  family = step_family$family,
                                  sl_risk = sl_risk, sl_coef = sl_coef,
                                  n_risk = n_at_risk, n_train = n_train,
                                  n_clipped_lower = n_clipped_lower,
                                  n_clipped_upper = n_clipped_upper,
                                  n_preds = length(preds))

      # Propagate pseudo-outcomes backward
      prev_times <- all_time_vals[all_time_vals < tt]
      if (length(prev_times) > 0 && length(all_ids) > 0) {
        prev_t <- max(prev_times)
        pred_dt_prop <- data.table::data.table(
          .tmp_id = all_ids,
          .tmp_pred = all_preds
        )
        data.table::setnames(pred_dt_prop, ".tmp_id", id_col)
        prev_rows <- dt_w[[time_col]] == prev_t
        prev_ids <- dt_w[[id_col]][prev_rows]
        match_idx <- match(prev_ids, pred_dt_prop[[id_col]])
        has_match <- !is.na(match_idx)
        dt_w$.longy_Q[which(prev_rows)[has_match]] <- pred_dt_prop$.tmp_pred[match_idx[has_match]]
      }

      # If this is the final backward step (earliest time), store predictions
      if (tt == backward_times[n_back] && length(all_ids) > 0) {
        return_preds <- data.table::data.table(
          .id = all_ids,
          .target_time = target_t,
          .Q_hat = all_preds
        )
      }
    }

    preds_out <- if (exists("return_preds", inherits = FALSE)) return_preds else NULL
    sl_info_out <- sl_info_target[!vapply(sl_info_target, is.null, logical(1))]
    list(regime = rname, predictions = preds_out, sl_info = sl_info_out)
  }

  # Strip obj from closure to avoid serializing the full longy_data object
  if (parallel) {
    one_task <- .clean_closure(one_task, c(
      "tasks", "regime", "regimes_list", "target_times",
      "dt", "nodes", "all_time_vals", "id_col", "time_col",
      "a_col", "y_col", "is_binary", "is_survival", "has_competing",
      "max_regime_lags", "covariates", "risk_set",
      "worker_ffSL", "worker_verbose",
      "min_obs", "bounds", "adaptive_cv", "learners"
    ))
  }

  if (verbose && parallel)
    .vmsg("  Q: dispatching %d regime x target tasks...", nrow(tasks))

  task_results <- .parallel_or_sequential(
    seq_len(nrow(tasks)), one_task, parallel = parallel,
    verbose = FALSE
  )

  # Regroup results by regime
  for (rname in regime) {
    rname_tasks <- Filter(function(x) identical(x$regime, rname), task_results)
    predictions_list <- lapply(rname_tasks, `[[`, "predictions")
    sl_info_list <- lapply(rname_tasks, `[[`, "sl_info")

    # Combine predictions
    non_null <- !vapply(predictions_list, is.null, logical(1))
    if (!any(non_null)) {
      warning(sprintf("No observations at risk for any time point in outcome model (regime=%s).",
                      rname), call. = FALSE)
    }
    all_preds <- data.table::rbindlist(predictions_list[non_null])
    data.table::setnames(all_preds, ".id", id_col)

    # Flatten sl_info
    sl_info <- unlist(sl_info_list, recursive = FALSE)
    sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

    obj$fits$outcome[[rname]] <- list(
      regime = rname,
      predictions = all_preds,
      covariates = covariates,
      learners = learners,
      bounds = bounds,
      use_ffSL = use_ffSL,
      sl_info = sl_info,
      family = "gaussian",
      risk_set = risk_set
    )
  }

  obj
}
