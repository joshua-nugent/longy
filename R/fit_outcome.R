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
#' @param verbose Logical. Print progress.
#'
#' @return Modified \code{longy_data} object with outcome fits stored in
#'   \code{obj$fits$outcome}.
#' @export
fit_outcome <- function(obj, regime, covariates = NULL, learners = NULL,
                        sl_control = list(), adaptive_cv = TRUE,
                        min_obs = 50L, bounds = c(0.005, 0.995),
                        times = NULL, verbose = TRUE) {
  stopifnot(inherits(obj, "longy_data"))

  if (!regime %in% names(obj$regimes)) {
    stop(sprintf("Regime '%s' not found. Use define_regime() first.", regime),
         call. = FALSE)
  }

  reg <- obj$regimes[[regime]]
  dt <- obj$data
  nodes <- obj$nodes
  all_time_vals <- obj$meta$time_values
  id_col <- nodes$id
  time_col <- nodes$time
  a_col <- nodes$treatment
  y_col <- nodes$outcome

  # Determine family based on outcome type
  is_binary <- nodes$outcome_type %in% c("binary", "survival")
  family <- if (is_binary) stats::binomial() else stats::gaussian()

  # Determine target times
  if (!is.null(times)) {
    target_times <- sort(all_time_vals[all_time_vals %in% times])
  } else {
    target_times <- all_time_vals
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying, nodes$treatment)
  }

  # Build cumulative regime-consistency and uncensored tracking
  dt <- .add_tracking_columns(dt, nodes, reg)

  # Get regime values for counterfactual prediction
  regime_vals <- .evaluate_regime(reg, dt)
  dt[, .longy_regime_a := regime_vals]

  # Outer loop: one backward pass per target time (ICE)
  predictions_list <- vector("list", length(target_times))
  sl_info_list <- vector("list", length(target_times))

  for (target_idx in seq_along(target_times)) {
    target_t <- target_times[target_idx]

    # Initialize Q: Y only at target time, NA elsewhere
    dt[, .longy_Q := NA_real_]
    dt[dt[[time_col]] == target_t, .longy_Q := as.numeric(get(y_col))]

    # Backward pass: from target_t down to min time
    backward_times <- rev(all_time_vals[all_time_vals <= target_t])
    n_back <- length(backward_times)
    sl_info_target <- vector("list", n_back)

    for (i in seq_along(backward_times)) {
      tt <- backward_times[i]
      dt_t <- dt[dt[[time_col]] == tt, ]

      # Risk set: ALL uncensored subjects through t
      still_in <- dt_t$.longy_cum_uncens == 1L
      still_in[is.na(still_in)] <- FALSE

      n_risk <- sum(still_in)
      if (n_risk == 0) {
        if (verbose) .vmsg("  Q target=%d time %d: 0 at risk, skipping",
                           target_t, tt)
        next
      }

      # Training subset: risk set AND Q is non-missing
      Q_at_t <- dt_t$.longy_Q[still_in]
      has_Q <- !is.na(Q_at_t)
      n_train <- sum(has_Q)

      if (n_train == 0) {
        if (verbose) .vmsg("  Q target=%d time %d: n_risk=%d, 0 non-missing Q, skipping",
                           target_t, tt, n_risk)
        next
      }

      # Covariates for subjects in risk set
      risk_ids <- dt_t[[id_col]][still_in]
      X_risk <- as.data.frame(dt_t[still_in, covariates, with = FALSE])

      # Training data
      X_train <- X_risk[has_Q, , drop = FALSE]
      Y_train <- Q_at_t[has_Q]

      # Extract sampling weights for training subjects (NULL if none)
      ow <- NULL
      if (!is.null(nodes$sampling_weights)) {
        ow_risk <- dt_t[[nodes$sampling_weights]][still_in]
        ow <- ow_risk[has_Q]
      }

      # Fit model
      if (n_train >= min_obs && length(unique(Y_train)) > 1) {
        cv_folds <- 10L
        if (adaptive_cv && is_binary) {
          cv_info <- .adaptive_cv_folds(Y_train)
          cv_folds <- cv_info$V
        }
        fit <- .safe_sl(Y = Y_train, X = X_train, family = family,
                        learners = learners, cv_folds = cv_folds,
                        obs_weights = ow, verbose = verbose)
        method <- fit$method
        sl_risk <- fit$sl_risk
        sl_coef <- fit$sl_coef

        # Counterfactual prediction: set A to regime value for ALL in risk set
        X_cf <- X_risk
        if (a_col %in% covariates) {
          X_cf[[a_col]] <- dt_t$.longy_regime_a[still_in]
        }

        # Predict using the fitted model
        if (method == "SuperLearner") {
          preds <- as.numeric(stats::predict(fit$fit, newdata = X_cf)$pred)
        } else if (method == "glm") {
          preds <- as.numeric(stats::predict(fit$fit, newdata = X_cf,
                                              type = "response"))
        } else {
          marg <- mean(Y_train)
          preds <- rep(marg, n_risk)
        }
      } else {
        # Marginal fallback
        if (!is.null(ow)) {
          marg <- stats::weighted.mean(Y_train, ow)
        } else {
          marg <- mean(Y_train)
        }
        preds <- rep(marg, n_risk)
        method <- "marginal"
        sl_risk <- NULL
        sl_coef <- NULL
      }

      # Bound predictions for binary/survival
      if (is_binary) {
        preds <- .bound(preds, bounds[1], bounds[2])
      }

      if (verbose) {
        .vmsg("  Q target=%d time %d: n_risk=%d, n_train=%d, mean_Q=%.3f, method=%s",
              target_t, tt, n_risk, n_train, mean(preds), method)
      }

      sl_info_target[[i]] <- list(time = tt, target_time = target_t,
                                  method = method,
                                  sl_risk = sl_risk, sl_coef = sl_coef,
                                  n_risk = n_risk, n_train = n_train)

      # Propagate pseudo-outcomes backward: set Q at time tt-1 for subjects
      # in the risk set at tt. Their Q at tt-1 becomes the prediction from tt.
      prev_times <- all_time_vals[all_time_vals < tt]
      if (length(prev_times) > 0) {
        prev_t <- max(prev_times)
        # For each subject in the risk set at tt, set their Q at prev_t
        pred_dt_prop <- data.table::data.table(
          .tmp_id = risk_ids,
          .tmp_pred = preds
        )
        data.table::setnames(pred_dt_prop, ".tmp_id", id_col)
        # Match by id in the prev_t rows of dt
        prev_rows <- dt[[time_col]] == prev_t
        prev_ids <- dt[[id_col]][prev_rows]
        match_idx <- match(prev_ids, pred_dt_prop[[id_col]])
        has_match <- !is.na(match_idx)
        dt$.longy_Q[which(prev_rows)[has_match]] <- pred_dt_prop$.tmp_pred[match_idx[has_match]]
      }

      # If this is the final backward step (earliest time), store predictions
      if (tt == backward_times[n_back]) {
        predictions_list[[target_idx]] <- data.table::data.table(
          .id = risk_ids,
          .target_time = target_t,
          .Q_hat = preds
        )
      }
    }

    sl_info_list[[target_idx]] <- sl_info_target[
      !vapply(sl_info_target, is.null, logical(1))
    ]
  }

  # Combine predictions
  non_null <- !vapply(predictions_list, is.null, logical(1))
  if (!any(non_null)) {
    warning("No observations at risk for any time point in outcome model.",
            call. = FALSE)
  }
  all_preds <- data.table::rbindlist(predictions_list[non_null])
  data.table::setnames(all_preds, ".id", id_col)

  # Flatten sl_info
  sl_info <- unlist(sl_info_list, recursive = FALSE)
  sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

  obj$fits$outcome <- list(
    regime = regime,
    predictions = all_preds,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    sl_info = sl_info,
    family = if (is_binary) "binomial" else "gaussian"
  )

  # Clean up tracking columns
  dt[, .longy_Q := NULL]
  dt[, .longy_regime_a := NULL]
  .remove_tracking_columns(obj$data)

  obj
}
