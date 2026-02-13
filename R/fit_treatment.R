#' Fit Treatment Models (g_A)
#'
#' Fits propensity score models P(A(t) = 1 | past) at each time point for
#' subjects in the risk set (regime-consistent and uncensored through t-1).
#'
#' @param obj A `longy_data` object with at least one regime defined.
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
#' @param bounds Numeric vector of length 2. Bounds for predicted probabilities.
#' @param verbose Logical. Print progress.
#'
#' @return Modified `longy_data` object with treatment fits stored.
#' @export
fit_treatment <- function(obj, regime, covariates = NULL, learners = NULL,
                          sl_control = list(), adaptive_cv = TRUE,
                          min_obs = 50L, bounds = c(0.005, 0.995),
                          verbose = TRUE) {
  stopifnot(inherits(obj, "longy_data"))

  if (!regime %in% names(obj$regimes)) {
    stop(sprintf("Regime '%s' not found. Use define_regime() first.", regime),
         call. = FALSE)
  }

  reg <- obj$regimes[[regime]]
  dt <- obj$data
  nodes <- obj$nodes
  time_vals <- obj$meta$time_values

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  # Build cumulative regime-consistency and uncensored tracking
  dt <- .add_tracking_columns(dt, nodes, reg)

  results <- vector("list", length(time_vals))
  sl_info <- vector("list", length(time_vals))

  for (i in seq_along(time_vals)) {
    tt <- time_vals[i]
    dt_t <- dt[dt[[nodes$time]] == tt, ]

    # Risk set: regime-consistent through t-1 AND uncensored through t-1
    if (i == 1) {
      still_in <- rep(TRUE, nrow(dt_t))
    } else {
      still_in <- dt_t$.longy_consist_prev == 1L & dt_t$.longy_uncens_prev == 1L
      still_in[is.na(still_in)] <- FALSE
    }

    n_risk <- sum(still_in)
    if (n_risk == 0) {
      if (verbose) .vmsg("  g_A time %d: 0 at risk, skipping", tt)
      next
    }

    X <- as.data.frame(dt_t[still_in, covariates, with = FALSE])
    Y <- dt_t[[nodes$treatment]][still_in]

    # Fit model or use marginal
    if (length(unique(Y)) > 1 && n_risk >= min_obs) {
      cv_folds <- 10L
      if (adaptive_cv) {
        cv_info <- .adaptive_cv_folds(Y)
        cv_folds <- cv_info$V
      }
      fit <- .safe_sl(Y = Y, X = X, learners = learners,
                      cv_folds = cv_folds, verbose = verbose)
      p_a <- .bound(fit$predictions, bounds[1], bounds[2])
      method <- fit$method
      sl_risk <- fit$sl_risk
      sl_coef <- fit$sl_coef
    } else {
      p_a <- .bound(rep(mean(Y), n_risk), bounds[1], bounds[2])
      method <- "marginal"
      sl_risk <- NULL
      sl_coef <- NULL
    }

    marg_a <- mean(Y)

    results[[i]] <- data.table::data.table(
      .id = dt_t[[nodes$id]][still_in],
      .time = tt,
      .n_risk = n_risk,
      .treatment = Y,
      .marg_a = marg_a,
      .p_a = p_a,
      .method = method
    )

    sl_info[[i]] <- list(time = tt, method = method,
                         sl_risk = sl_risk, sl_coef = sl_coef)

    if (verbose) {
      .vmsg("  g_A time %d: n_risk=%d, marg=%.3f, method=%s",
            tt, n_risk, marg_a, method)
    }
  }

  results <- data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
  data.table::setnames(results, ".id", nodes$id)

  sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

  obj$fits$treatment <- list(
    regime = regime,
    predictions = results,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    sl_info = sl_info
  )

  # Clean up tracking columns
  .remove_tracking_columns(obj$data)

  obj
}

#' Add cumulative regime-consistency and uncensored tracking columns
#'
#' These are used by all three nuisance model fitting functions.
#' @param dt data.table
#' @param nodes List of node names
#' @param regime Regime list
#' @return Modified data.table (by reference)
#' @noRd
.add_tracking_columns <- function(dt, nodes, regime) {
  id_col <- nodes$id
  time_col <- nodes$time
  a_col <- nodes$treatment
  c_cols <- nodes$censoring

  data.table::setkeyv(dt, c(id_col, time_col))

  # Regime-consistent indicator at each time
  if (regime$type == "static") {
    if (regime$value == 1L) {
      dt[, .longy_regime_consist := as.integer(get(a_col) == 1L)]
    } else {
      dt[, .longy_regime_consist := as.integer(get(a_col) == 0L)]
    }
  } else {
    # Dynamic/stochastic: evaluate regime
    regime_vals <- .evaluate_regime(regime, dt)
    if (regime$type == "dynamic") {
      dt[, .longy_regime_consist := as.integer(get(a_col) == regime_vals)]
    } else {
      # Stochastic: always "consistent" (probabilities handled in weights)
      dt[, .longy_regime_consist := 1L]
    }
  }

  # Cumulative consistency through current time
  dt[, .longy_cum_consist := cumprod(.longy_regime_consist), by = c(id_col)]

  # Cumulative uncensored through current time
  if (!is.null(c_cols) && length(c_cols) > 0) {
    # Uncensored at each time = all censoring indicators are 0
    dt[, .longy_uncens := as.integer(Reduce(`+`, .SD) == 0L),
       .SDcols = c_cols]
    dt[, .longy_cum_uncens := cumprod(.longy_uncens), by = c(id_col)]
  } else {
    dt[, .longy_uncens := 1L]
    dt[, .longy_cum_uncens := 1L]
  }

  # Lagged versions (through t-1)
  dt[, .longy_consist_prev := data.table::shift(.longy_cum_consist, 1L,
                                                 fill = 1L, type = "lag"),
     by = c(id_col)]
  dt[, .longy_uncens_prev := data.table::shift(.longy_cum_uncens, 1L,
                                                fill = 1L, type = "lag"),
     by = c(id_col)]

  dt
}

#' Remove tracking columns added by .add_tracking_columns
#' @noRd
.remove_tracking_columns <- function(dt) {
  track_cols <- grep("^\\.longy_", names(dt), value = TRUE)
  if (length(track_cols) > 0) {
    dt[, (track_cols) := NULL]
  }
  invisible(dt)
}
