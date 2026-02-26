#' Cross-Fitting Infrastructure
#'
#' Sample-splitting (cross-fitting) for nuisance model fitting and TMLE
#' estimation. When enabled, nuisance models are fit on training folds and
#' predictions are made on validation folds, eliminating overfitting bias
#' in EIF-based inference.
#'
#' @name crossfit
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Cross-fitted treatment model
# ---------------------------------------------------------------------------

#' Cross-fitted treatment model (g_A)
#'
#' For each fold k, at each time t: fit g_A on training fold, predict on
#' validation fold. Marginal rates are computed from the FULL risk set.
#'
#' @inheritParams fit_treatment
#' @return Modified \code{longy_data} object with cross-fitted treatment
#'   predictions stored.
#' @noRd
.cf_fit_treatment <- function(obj, regime, covariates = NULL, learners = NULL,
                               sl_control = list(), adaptive_cv = TRUE,
                               min_obs = 50L, min_events = 20L,
                               bounds = c(0.005, 0.995),
                               times = NULL, sl_fn = "SuperLearner",
                               verbose = TRUE, risk_set = "all") {
  reg <- obj$regimes[[regime]]
  dt <- obj$data
  nodes <- obj$nodes
  id_col <- nodes$id
  time_vals <- obj$meta$time_values
  if (!is.null(times)) {
    time_vals <- time_vals[time_vals <= max(times)]
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  fold_col <- obj$crossfit$fold_id
  n_folds <- obj$crossfit$n_folds

  dt <- .add_tracking_columns(dt, nodes, reg)

  results <- vector("list", length(time_vals))

  for (i in seq_along(time_vals)) {
    tt <- time_vals[i]
    dt_t <- dt[dt[[nodes$time]] == tt, ]

    # Risk set: uncensored through t-1
    if (i == 1) {
      still_in <- rep(TRUE, nrow(dt_t))
    } else {
      still_in <- dt_t$.longy_uncens_prev == 1L
      still_in[is.na(still_in)] <- FALSE
    }

    # Followers-only restriction: limit to regime-consistent subjects
    if (risk_set == "followers") {
      consist_prev <- dt_t$.longy_consist_prev == 1L
      consist_prev[is.na(consist_prev)] <- FALSE
      still_in <- still_in & consist_prev
    }

    n_risk <- sum(still_in)
    if (n_risk == 0) {
      if (verbose) .vmsg("  CF g_A time %d: 0 at risk, skipping", tt)
      next
    }

    ids_risk <- dt_t[[id_col]][still_in]
    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(covariates, lag_covs)
    X_risk <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
    Y_risk <- dt_t[[nodes$treatment]][still_in]
    folds_risk <- dt_t[[fold_col]][still_in]

    # Sampling weights for at-risk subjects
    ow_risk <- NULL
    if (!is.null(nodes$sampling_weights)) {
      ow_risk <- dt_t[[nodes$sampling_weights]][still_in]
    }

    # Marginal rate from FULL risk set (population constant)
    marg_a <- if (!is.null(ow_risk)) {
      stats::weighted.mean(Y_risk, ow_risk)
    } else {
      mean(Y_risk)
    }

    # Cross-fitted predictions
    p_a <- rep(NA_real_, n_risk)

    for (k in seq_len(n_folds)) {
      val_idx <- which(folds_risk == k)
      train_idx <- which(folds_risk != k)

      if (length(val_idx) == 0) next

      Y_train <- Y_risk[train_idx]
      X_train <- X_risk[train_idx, , drop = FALSE]
      X_val <- X_risk[val_idx, , drop = FALSE]

      ow_train <- if (!is.null(ow_risk)) ow_risk[train_idx] else NULL

      n_minority_train <- min(sum(Y_train == 1), sum(Y_train == 0))
      minority_rate_train <- min(mean(Y_train), 1 - mean(Y_train))
      rare_events_train <- n_minority_train < min_events && minority_rate_train < 0.01
      if (length(train_idx) >= min_obs && length(unique(Y_train)) > 1 && !rare_events_train) {
        cv_folds <- 10L
        if (adaptive_cv) {
          cv_info <- .adaptive_cv_folds(Y_train)
          cv_folds <- cv_info$V
        }
        ctx <- sprintf("CF g_A, time=%d, fold=%d/%d, n_train=%d",
                       tt, k, n_folds, length(Y_train))
        fit <- .safe_sl(Y = Y_train, X = X_train, learners = learners,
                        cv_folds = cv_folds, obs_weights = ow_train,
                        sl_fn = sl_fn, context = ctx, verbose = FALSE)
        preds_val <- .predict_from_fit(fit, X_val)
      } else {
        marg_train <- if (!is.null(ow_train)) {
          stats::weighted.mean(Y_train, ow_train)
        } else {
          mean(Y_train)
        }
        preds_val <- rep(marg_train, length(val_idx))
      }

      p_a[val_idx] <- .bound(preds_val, bounds[1], bounds[2])
    }

    results[[i]] <- data.table::data.table(
      .id = ids_risk,
      .time = tt,
      .n_risk = n_risk,
      .treatment = Y_risk,
      .marg_a = marg_a,
      .p_a = p_a,
      .method = "crossfit"
    )

    if (verbose) {
      rs_label <- if (risk_set == "followers") {
        sprintf(" (%s followers)", regime)
      } else {
        ""
      }
      .vmsg("  CF g_A time %d%s: n_risk=%d, marg=%.3f, K=%d",
            tt, rs_label, n_risk, marg_a, n_folds)
    }
  }

  non_null <- !vapply(results, is.null, logical(1))
  if (!any(non_null)) {
    warning("No observations at risk for any time point in CF treatment model.",
            call. = FALSE)
  }
  results <- data.table::rbindlist(results[non_null])
  data.table::setnames(results, ".id", nodes$id)

  obj$fits$treatment[[regime]] <- list(
    regime = regime,
    predictions = results,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    sl_fn = sl_fn,
    sl_info = list(),
    crossfit = TRUE,
    risk_set = risk_set
  )

  .remove_tracking_columns(obj$data)
  obj
}

# ---------------------------------------------------------------------------
# Cross-fitted censoring model
# ---------------------------------------------------------------------------

#' Cross-fitted censoring model (g_C)
#' @inheritParams fit_censoring
#' @noRd
.cf_fit_censoring <- function(obj, regime, covariates = NULL, learners = NULL,
                               sl_control = list(), adaptive_cv = TRUE,
                               min_obs = 50L, min_events = 20L,
                               bounds = c(0.005, 0.995),
                               times = NULL, sl_fn = "SuperLearner",
                               verbose = TRUE) {
  nodes <- obj$nodes

  if (is.null(nodes$censoring) || length(nodes$censoring) == 0) {
    if (verbose) .vmsg("  No censoring columns defined, skipping CF fit_censoring()")
    obj$fits$censoring[[regime]] <- list()
    return(obj)
  }

  reg <- obj$regimes[[regime]]
  dt <- obj$data
  id_col <- nodes$id
  time_vals <- obj$meta$time_values
  if (!is.null(times)) {
    time_vals <- time_vals[time_vals <= max(times)]
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  fold_col <- obj$crossfit$fold_id
  n_folds <- obj$crossfit$n_folds

  dt <- .add_tracking_columns(dt, nodes, reg)

  # Initialize per-regime censoring list
  obj$fits$censoring[[regime]] <- list()

  for (cvar in nodes$censoring) {
    if (verbose) .vmsg("  CF fitting censoring model for '%s'...", cvar)

    results <- vector("list", length(time_vals))

    for (i in seq_along(time_vals)) {
      tt <- time_vals[i]
      dt_t <- dt[dt[[nodes$time]] == tt, ]

      if (i == 1) {
        still_in <- rep(TRUE, nrow(dt_t))
      } else {
        still_in <- dt_t$.longy_uncens_prev == 1L
        still_in[is.na(still_in)] <- FALSE
      }

      n_risk <- sum(still_in)
      if (n_risk == 0) {
        if (verbose) .vmsg("    CF g_C(%s) time %d: 0 at risk, skipping", cvar, tt)
        next
      }

      ids_risk <- dt_t[[id_col]][still_in]
      lag_covs <- .get_lag_covariates(nodes, i)
      all_covs <- c(covariates, lag_covs)
      X_risk <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
      Y_risk <- 1L - dt_t[[cvar]][still_in]  # P(uncensored)
      folds_risk <- dt_t[[fold_col]][still_in]

      ow_risk <- NULL
      if (!is.null(nodes$sampling_weights)) {
        ow_risk <- dt_t[[nodes$sampling_weights]][still_in]
      }

      marg_c <- if (!is.null(ow_risk)) {
        stats::weighted.mean(Y_risk, ow_risk)
      } else {
        mean(Y_risk)
      }

      p_c <- rep(NA_real_, n_risk)

      for (k in seq_len(n_folds)) {
        val_idx <- which(folds_risk == k)
        train_idx <- which(folds_risk != k)

        if (length(val_idx) == 0) next

        Y_train <- Y_risk[train_idx]
        X_train <- X_risk[train_idx, , drop = FALSE]
        X_val <- X_risk[val_idx, , drop = FALSE]

        ow_train <- if (!is.null(ow_risk)) ow_risk[train_idx] else NULL

        n_minority_train <- min(sum(Y_train == 1), sum(Y_train == 0))
        minority_rate_train <- min(mean(Y_train), 1 - mean(Y_train))
        rare_events_train <- n_minority_train < min_events && minority_rate_train < 0.01
        if (length(train_idx) >= min_obs && length(unique(Y_train)) > 1 && !rare_events_train) {
          cv_folds <- 10L
          if (adaptive_cv) {
            cv_info <- .adaptive_cv_folds(Y_train)
            cv_folds <- cv_info$V
          }
          ctx <- sprintf("CF g_C(%s), time=%d, fold=%d/%d, n_train=%d",
                         cvar, tt, k, n_folds, length(Y_train))
          fit <- .safe_sl(Y = Y_train, X = X_train, learners = learners,
                          cv_folds = cv_folds, obs_weights = ow_train,
                          sl_fn = sl_fn, context = ctx, verbose = FALSE)
          preds_val <- .predict_from_fit(fit, X_val)
        } else {
          marg_train <- if (!is.null(ow_train)) {
            stats::weighted.mean(Y_train, ow_train)
          } else {
            mean(Y_train)
          }
          preds_val <- rep(marg_train, length(val_idx))
        }

        p_c[val_idx] <- .bound(preds_val, bounds[1], bounds[2])
      }

      results[[i]] <- data.table::data.table(
        .id = ids_risk,
        .time = tt,
        .n_risk = n_risk,
        .censored = dt_t[[cvar]][still_in],
        .marg_c = marg_c,
        .p_c = p_c,
        .method = "crossfit"
      )

      if (verbose) {
        .vmsg("    CF g_C(%s) time %d: n_risk=%d, marg_uncens=%.3f, K=%d",
              cvar, tt, n_risk, marg_c, n_folds)
      }
    }

    non_null <- !vapply(results, is.null, logical(1))
    if (!any(non_null)) {
      warning(sprintf(
        "No observations at risk for any time point in CF censoring model for '%s'.",
        cvar), call. = FALSE)
    }
    results <- data.table::rbindlist(results[non_null])
    data.table::setnames(results, ".id", nodes$id)

    obj$fits$censoring[[regime]][[cvar]] <- list(
      predictions = results,
      covariates = covariates,
      learners = learners,
      bounds = bounds,
      sl_fn = sl_fn,
      sl_info = list(),
      crossfit = TRUE
    )
  }

  .remove_tracking_columns(obj$data)
  obj
}

# ---------------------------------------------------------------------------
# Cross-fitted observation model
# ---------------------------------------------------------------------------

#' Cross-fitted observation model (g_R)
#' @inheritParams fit_observation
#' @noRd
.cf_fit_observation <- function(obj, regime, covariates = NULL, learners = NULL,
                                 sl_control = list(), adaptive_cv = TRUE,
                                 min_obs = 50L, min_events = 20L,
                                 bounds = c(0.005, 0.995),
                                 times = NULL, sl_fn = "SuperLearner",
                                 verbose = TRUE) {
  nodes <- obj$nodes

  if (is.null(nodes$observation)) {
    if (verbose) .vmsg("  No observation column defined, skipping CF fit_observation()")
    return(obj)
  }

  reg <- obj$regimes[[regime]]
  dt <- obj$data
  id_col <- nodes$id
  time_vals <- obj$meta$time_values
  if (!is.null(times)) {
    time_vals <- time_vals[time_vals <= max(times)]
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  fold_col <- obj$crossfit$fold_id
  n_folds <- obj$crossfit$n_folds

  dt <- .add_tracking_columns(dt, nodes, reg)

  results <- vector("list", length(time_vals))

  for (i in seq_along(time_vals)) {
    tt <- time_vals[i]
    dt_t <- dt[dt[[nodes$time]] == tt, ]

    if (i == 1) {
      still_in <- rep(TRUE, nrow(dt_t))
    } else {
      still_in <- dt_t$.longy_uncens_prev == 1L
      still_in[is.na(still_in)] <- FALSE
    }

    # Must be uncensored at t (observation requires being present)
    still_in <- still_in & dt_t$.longy_uncens == 1L

    n_risk <- sum(still_in)
    if (n_risk == 0) {
      if (verbose) .vmsg("  CF g_R time %d: 0 at risk, skipping", tt)
      next
    }

    ids_risk <- dt_t[[id_col]][still_in]
    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(covariates, lag_covs)
    X_risk <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
    Y_risk <- dt_t[[nodes$observation]][still_in]
    folds_risk <- dt_t[[fold_col]][still_in]

    ow_risk <- NULL
    if (!is.null(nodes$sampling_weights)) {
      ow_risk <- dt_t[[nodes$sampling_weights]][still_in]
    }

    marg_r <- if (!is.null(ow_risk)) {
      stats::weighted.mean(Y_risk, ow_risk)
    } else {
      mean(Y_risk)
    }

    p_r <- rep(NA_real_, n_risk)

    for (k in seq_len(n_folds)) {
      val_idx <- which(folds_risk == k)
      train_idx <- which(folds_risk != k)

      if (length(val_idx) == 0) next

      Y_train <- Y_risk[train_idx]
      X_train <- X_risk[train_idx, , drop = FALSE]
      X_val <- X_risk[val_idx, , drop = FALSE]

      ow_train <- if (!is.null(ow_risk)) ow_risk[train_idx] else NULL

      n_minority_train <- min(sum(Y_train == 1), sum(Y_train == 0))
      minority_rate_train <- min(mean(Y_train), 1 - mean(Y_train))
      rare_events_train <- n_minority_train < min_events && minority_rate_train < 0.01
      if (length(train_idx) >= min_obs && length(unique(Y_train)) > 1 && !rare_events_train) {
        cv_folds <- 10L
        if (adaptive_cv) {
          cv_info <- .adaptive_cv_folds(Y_train)
          cv_folds <- cv_info$V
        }
        ctx <- sprintf("CF g_R, time=%d, fold=%d/%d, n_train=%d",
                       tt, k, n_folds, length(Y_train))
        fit <- .safe_sl(Y = Y_train, X = X_train, learners = learners,
                        cv_folds = cv_folds, obs_weights = ow_train,
                        sl_fn = sl_fn, context = ctx, verbose = FALSE)
        preds_val <- .predict_from_fit(fit, X_val)
      } else {
        marg_train <- if (!is.null(ow_train)) {
          stats::weighted.mean(Y_train, ow_train)
        } else {
          mean(Y_train)
        }
        preds_val <- rep(marg_train, length(val_idx))
      }

      p_r[val_idx] <- .bound(preds_val, bounds[1], bounds[2])
    }

    results[[i]] <- data.table::data.table(
      .id = ids_risk,
      .time = tt,
      .n_risk = n_risk,
      .observed = Y_risk,
      .marg_r = marg_r,
      .p_r = p_r,
      .method = "crossfit"
    )

    if (verbose) {
      .vmsg("  CF g_R time %d: n_risk=%d, marg=%.3f, K=%d",
            tt, n_risk, marg_r, n_folds)
    }
  }

  non_null <- !vapply(results, is.null, logical(1))
  if (!any(non_null)) {
    warning("No observations at risk for any time point in CF observation model.",
            call. = FALSE)
  }
  results <- data.table::rbindlist(results[non_null])
  data.table::setnames(results, ".id", nodes$id)

  obj$fits$observation[[regime]] <- list(
    regime = regime,
    predictions = results,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    sl_fn = sl_fn,
    sl_info = list(),
    crossfit = TRUE
  )

  .remove_tracking_columns(obj$data)
  obj
}

# ---------------------------------------------------------------------------
# Cross-fitted Q step helper (for TMLE backward pass)
# ---------------------------------------------------------------------------

#' Cross-fitted Q model at a single backward step
#'
#' Fits Q on each training fold and predicts for the corresponding validation
#' fold. Returns cross-fitted Q_bar for all risk-set subjects.
#'
#' @param Y_train_all Numeric vector. Pseudo-outcomes for all trainable subjects
#'   in the risk set (non-NA Q).
#' @param X_risk data.frame. Covariates for ALL subjects in the risk set.
#' @param has_Q Logical vector (length = nrow(X_risk)). TRUE for subjects with
#'   non-missing Q.
#' @param folds_risk Integer vector. Fold assignments for risk-set subjects.
#' @param n_folds Integer. Number of folds.
#' @param family Model family.
#' @param learners SuperLearner library.
#' @param cv_folds Integer. Number of CV folds within each training fold.
#' @param sl_fn Character. SuperLearner implementation.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param regime_cf Named list. Column names mapping to counterfactual regime
#'   values (vectors of length nrow(X_risk)) for treatment and lagged treatment
#'   columns.
#'
#' @return Numeric vector of cross-fitted Q_bar, length = nrow(X_risk).
#' @noRd
.cf_fit_q_step <- function(Y_train_all, X_risk, has_Q, folds_risk, n_folds,
                            family, learners, cv_folds = 10L,
                            sl_fn = "SuperLearner", min_obs = 50L,
                            regime_cf = NULL) {
  n_risk <- nrow(X_risk)
  Q_bar <- rep(NA_real_, n_risk)

  for (k in seq_len(n_folds)) {
    val_idx <- which(folds_risk == k)
    train_idx <- which(folds_risk != k)

    if (length(val_idx) == 0) next

    # Training data: train fold subjects with non-NA Q
    train_has_Q <- has_Q[train_idx]
    train_Q_idx <- train_idx[train_has_Q]

    if (length(train_Q_idx) == 0) {
      # Fall back to overall mean of all subjects with Q
      Q_bar[val_idx] <- mean(Y_train_all)
      next
    }

    Y_train <- Y_train_all[train_has_Q[seq_along(train_idx)]]
    # Recompute: Y for training fold subjects that have Q
    # train_Q_idx gives positions in the risk set
    Y_k <- numeric(0)
    # We need to map from risk-set positions to the Y_train_all vector
    # Y_train_all contains values for all has_Q == TRUE subjects in risk set
    # Need to find which of those are in the training fold
    all_Q_idx <- which(has_Q)
    train_Q_mask <- folds_risk[all_Q_idx] != k
    Y_k <- Y_train_all[train_Q_mask]
    X_k <- X_risk[all_Q_idx[train_Q_mask], , drop = FALSE]

    if (length(Y_k) < 2) {
      Q_bar[val_idx] <- mean(Y_train_all)
      next
    }

    # Counterfactual covariates for validation fold
    X_val_cf <- X_risk[val_idx, , drop = FALSE]
    if (!is.null(regime_cf)) {
      for (cf_col in names(regime_cf)) {
        if (cf_col %in% names(X_val_cf)) {
          X_val_cf[[cf_col]] <- regime_cf[[cf_col]][val_idx]
        }
      }
    }

    if (length(Y_k) >= min_obs && length(unique(Y_k)) > 1) {
      ctx <- sprintf("CF Q-step, fold=%d/%d, n_train=%d", k, n_folds, length(Y_k))
      fit <- .safe_sl(Y = Y_k, X = X_k, family = family,
                      learners = learners, cv_folds = cv_folds,
                      sl_fn = sl_fn, context = ctx, verbose = FALSE)
      preds <- .predict_from_fit(fit, X_val_cf)
    } else {
      preds <- rep(mean(Y_k), length(val_idx))
    }

    Q_bar[val_idx] <- preds
  }

  Q_bar
}

# ---------------------------------------------------------------------------
# Cross-fitted TMLE estimator
# ---------------------------------------------------------------------------

#' Cross-fitted TMLE estimation
#'
#' Implements CV-TMLE (Zheng & van der Laan 2011) using pooled fluctuation.
#' Nuisance g models are already cross-fitted. Q models are cross-fitted
#' within each backward step. The single fluctuation epsilon is fit on all
#' n subjects (safe because epsilon is scalar).
#'
#' @inheritParams estimate_tmle
#' @return An S3 object of class \code{"longy_result"}.
#' @noRd
.cf_estimate_tmle <- function(obj, regime, times = NULL, inference = "eif",
                               ci_level = 0.95, n_boot = 200L,
                               g_bounds = c(0.01, 1), outcome_range = NULL,
                               verbose = TRUE) {
  nodes <- obj$nodes
  dt <- obj$data
  id_col <- nodes$id
  time_col <- nodes$time
  a_col <- nodes$treatment
  y_col <- nodes$outcome
  all_time_vals <- obj$meta$time_values
  reg <- obj$regimes[[regime]]

  is_binary <- nodes$outcome_type %in% c("binary", "survival")
  is_survival <- nodes$outcome_type == "survival"

  outcome_settings <- obj$fits$outcome[[regime]]
  covariates <- outcome_settings$covariates
  learners <- outcome_settings$learners
  q_bounds <- outcome_settings$bounds
  sl_fn <- if (!is.null(outcome_settings$sl_fn)) outcome_settings$sl_fn else "SuperLearner"
  min_obs <- 50L

  fold_col <- obj$crossfit$fold_id
  n_folds <- obj$crossfit$n_folds

  if (!is.null(times)) {
    target_times <- sort(all_time_vals[all_time_vals %in% times])
  } else {
    target_times <- all_time_vals
  }

  # Compute cross-fitted cumulative g (uses already cross-fitted g_A, g_C, g_R)
  g_cum_dt <- .compute_cumulative_g(obj, regime = regime, g_bounds = g_bounds)

  # Outcome scaling
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
    if (y_range_width == 0) y_range_width <- 1
  } else {
    y_min <- 0
    y_max <- 1
    y_range_width <- 1
  }

  dt <- .add_tracking_columns(dt, nodes, reg)
  regime_vals <- .evaluate_regime(reg, dt)
  dt[, .longy_regime_a := regime_vals]

  # Pre-compute lagged regime values for counterfactual treatment history
  k_val <- if (is.null(nodes$lag_k)) 0 else nodes$lag_k
  max_regime_lags <- if (is.infinite(k_val)) length(all_time_vals) - 1 else min(k_val, length(all_time_vals) - 1)
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

  eps <- q_bounds[1]

  est_list <- vector("list", length(target_times))

  for (target_idx in seq_along(target_times)) {
    target_t <- target_times[target_idx]

    if (verbose) .vmsg("CF-TMLE target time %d:", target_t)

    # Initialize Q
    dt[, .longy_Q := NA_real_]
    dt[dt[[time_col]] == target_t,
       .longy_Q := {
         y_raw <- as.numeric(get(y_col))
         if (!is_binary) (y_raw - y_min) / y_range_width else y_raw
       }]

    backward_times <- rev(all_time_vals[all_time_vals <= target_t])
    n_back <- length(backward_times)

    Q_star_list <- vector("list", n_back)
    names(Q_star_list) <- as.character(backward_times)

    for (i in seq_along(backward_times)) {
      tt <- backward_times[i]
      dt_t <- dt[dt[[time_col]] == tt, ]

      # Uncensored subjects through C(t) (Convention B: C before Y).
      still_in <- dt_t$.longy_cum_uncens == 1L
      still_in[is.na(still_in)] <- FALSE

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
        if (verbose) .vmsg("  CF-TMLE time %d: 0 at risk, skipping", tt)
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

      fluct_epsilon <- 0
      Q_star <- numeric(0)

      if (n_train > 0) {
        time_index <- match(tt, all_time_vals)
        lag_covs <- .get_lag_covariates(nodes, time_index)
        all_covs <- c(covariates, lag_covs)
        X_risk <- as.data.frame(dt_t[at_risk, all_covs, with = FALSE])
        folds_risk <- dt_t[[fold_col]][at_risk]

        # Cross-fitted Q_bar
        Y_train_all <- Q_at_t[has_Q]

        # Counterfactual regime values for at-risk subjects (current + lagged)
        regime_cf <- list()
        if (a_col %in% all_covs) {
          regime_cf[[a_col]] <- dt_t$.longy_regime_a[at_risk]
        }
        trt_lag_prefix <- paste0(".longy_lag_", a_col, "_")
        for (lc in all_covs[startsWith(all_covs, trt_lag_prefix)]) {
          regime_lc <- sub(trt_lag_prefix,
                           paste0(".longy_lag_regime_", a_col, "_"), lc, fixed = TRUE)
          if (regime_lc %in% names(dt_t)) {
            regime_cf[[lc]] <- dt_t[[regime_lc]][at_risk]
          }
        }

        # TMLE Q-models always use quasibinomial: predictions feed into
        # qlogis(Q_bar) for fluctuation, so they MUST be in (0,1).
        # gaussian can produce predictions outside [0,1] that clip to
        # bounds, creating extreme logit values and huge epsilon.
        # .safe_sl() handles quasibinomial→binomial swap for SL compatibility.
        step_family <- stats::quasibinomial()

        Q_bar <- .cf_fit_q_step(
          Y_train_all = Y_train_all,
          X_risk = X_risk,
          has_Q = has_Q,
          folds_risk = folds_risk,
          n_folds = n_folds,
          family = step_family,
          learners = learners,
          cv_folds = 10L,
          sl_fn = sl_fn,
          min_obs = min_obs,
          regime_cf = regime_cf
        )

        # Bound Q_bar
        Q_bar <- .bound(Q_bar, eps, 1 - eps)

        # --- Pooled fluctuation (at-risk subjects only) ---
        cum_consist <- dt_t$.longy_cum_consist[at_risk] == 1L
        cum_consist[is.na(cum_consist)] <- FALSE
        in_fluct <- cum_consist

        # Get cross-fitted g_cum and g_r (Convention B: g_cum includes g_c(t))
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

        pseudo_out <- dt_t$.longy_Q[at_risk]

        # Pooled fluctuation: single epsilon over at-risk subjects
        fluct <- .tmle_fluctuate(Q_bar = Q_bar,
                                 pseudo_outcome = pseudo_out,
                                 g_cum = g_denom,
                                 in_fluctuation_set = in_fluct,
                                 bounds = c(eps, 1 - eps))
        Q_star <- fluct$Q_star
        fluct_epsilon <- fluct$epsilon
      } # end if (n_train > 0)

      if (verbose) {
        .vmsg("  CF-TMLE time %d: n_at_risk=%d, n_train=%d, n_primary_abs=%d, n_competing_abs=%d, eps=%.4f",
              tt, n_at_risk, n_train, n_primary, n_competing, fluct_epsilon)
      }

      # --- Predict Q* for newly-censored subjects (needed for EIF) ---
      # Subjects who were at-risk at the previous time but dropped out now.
      # Without Q* predictions for them, the EIF defaults missing Q*_{s+1}
      # to 0, causing enormous spurious augmentation.
      cens_ids <- dt_t[[id_col]][integer(0)]  # empty, same type as id column
      Q_star_cens <- numeric(0)

      if (n_train > 0) {
        # Detect subjects newly censored at t (Convention B: at-risk = uncensored
        # through C(t), so newly-censored = uncensored through C(t-1) but not C(t))
        newly_cens <- !at_risk & !absorbed_primary & !absorbed_competing &
          (dt_t$.longy_uncens_prev == 1L | dt_t$.longy_consist_prev == 1L)
        newly_cens[is.na(newly_cens)] <- FALSE
        newly_cens <- newly_cens & !at_risk
        n_newly_cens <- sum(newly_cens)

        if (n_newly_cens > 0) {
          # all_covs already computed above (covariates + lag columns)
          X_cens <- as.data.frame(dt_t[newly_cens, all_covs, with = FALSE])
          X_cens_cf <- X_cens
          if (a_col %in% all_covs) {
            X_cens_cf[[a_col]] <- dt_t$.longy_regime_a[newly_cens]
          }
          trt_lag_prefix_c <- paste0(".longy_lag_", a_col, "_")
          for (lc in all_covs[startsWith(all_covs, trt_lag_prefix_c)]) {
            regime_lc <- sub(trt_lag_prefix_c,
                             paste0(".longy_lag_regime_", a_col, "_"), lc, fixed = TRUE)
            if (regime_lc %in% names(dt_t)) {
              X_cens_cf[[lc]] <- dt_t[[regime_lc]][newly_cens]
            }
          }

          # Fit one model on all at-risk training data (non-cross-fitted)
          # for newly-censored prediction — sufficient for EIF correction
          all_Q_idx <- which(has_Q)
          Y_all <- Q_at_t[has_Q]
          X_all <- as.data.frame(dt_t[at_risk, all_covs, with = FALSE])[all_Q_idx, , drop = FALSE]

          if (length(Y_all) >= min_obs && length(unique(Y_all)) > 1) {
            cens_fit <- .safe_sl(Y = Y_all, X = X_all,
                                 family = step_family,
                                 learners = learners,
                                 cv_folds = min(10L, length(Y_all)),
                                 sl_fn = sl_fn,
                                 context = "CF newly-censored Q",
                                 verbose = FALSE)
            Q_bar_cens <- .predict_from_fit(cens_fit, X_cens_cf)
          } else {
            Q_bar_cens <- rep(mean(Y_all), n_newly_cens)
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

      if (length(all_ids) > 0) {
        Q_star_list[[as.character(tt)]] <- data.table::data.table(
          .tmp_id = all_ids,
          .Q_star = all_Q_star
        )
        data.table::setnames(Q_star_list[[as.character(tt)]], ".tmp_id", id_col)
      }

      # Propagate Q*
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
        rows_to_update <- which(prev_rows)[has_match]
        vals_to_set <- pred_dt_prop$.tmp_pred[match_idx[has_match]]
        data.table::set(dt, i = rows_to_update, j = ".longy_Q",
                        value = vals_to_set)
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

    psi_scaled <- mean(Q_star_0$.Q_star)
    psi_hat <- if (!is_binary) psi_scaled * y_range_width + y_min else psi_scaled
    n_at_risk <- nrow(Q_star_0)

    est_row <- data.table::data.table(
      time = target_t,
      estimate = psi_hat,
      n_effective = as.numeric(n_at_risk),
      n_at_risk = n_at_risk
    )

    # EIF inference (uses cross-fitted g and Q*)
    if (inference == "eif") {
      D_i <- .compute_tmle_eif(
        Q_star_list = Q_star_list,
        g_cum_dt = g_cum_dt,
        obj = obj,
        target_t = target_t,
        psi_hat = psi_scaled,
        regime = regime,
        y_min = y_min,
        y_range_width = y_range_width
      )

      if (!is_binary) {
        D_i <- D_i * y_range_width
      }

      n_i <- length(D_i)
      if (n_i >= 2 && stats::var(D_i) > 0) {
        se <- sqrt(stats::var(D_i) / n_i)
      } else {
        se <- NA_real_
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
      obj, regime = regime, times = target_times,
      n_boot = n_boot, ci_level = ci_level,
      g_bounds = g_bounds, outcome_range = outcome_range,
      verbose = verbose
    )
    for (col in c("se", "ci_lower", "ci_upper")) {
      if (col %in% names(estimates)) estimates[, (col) := NULL]
    }
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival
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
    regime = regime,
    estimator = "tmle",
    inference = inference,
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(regime, "_tmle")]] <- result
  obj
}
