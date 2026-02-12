#' Fit Censoring Models (g_C)
#'
#' Fits models for P(C(t) = 0 | past) at each time point. The risk set
#' additionally requires regime-consistency AT time t (the subject has already
#' received A(t) before C(t) is determined).
#'
#' @param obj A `longy_data` object with at least one regime defined.
#' @param regime Character. Name of the regime.
#' @param covariates Character vector. Predictor columns. If NULL, uses all
#'   baseline + timevarying covariates.
#' @param learners Character vector. SuperLearner library.
#' @param sl_control List. Additional SuperLearner arguments.
#' @param adaptive_cv Logical. Adaptive CV fold selection.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param bounds Numeric(2). Prediction bounds.
#' @param verbose Logical. Progress messages.
#'
#' @return Modified `longy_data` object with censoring fits stored.
#' @export
fit_censoring <- function(obj, regime, covariates = NULL, learners = NULL,
                          sl_control = list(), adaptive_cv = TRUE,
                          min_obs = 50L, bounds = c(0.005, 0.995),
                          verbose = TRUE) {
  stopifnot(inherits(obj, "longy_data"))

  nodes <- obj$nodes

  # If no censoring columns, nothing to do
  if (is.null(nodes$censoring) || length(nodes$censoring) == 0) {
    if (verbose) .vmsg("  No censoring columns defined, skipping fit_censoring()")
    obj$fits$censoring <- list()
    return(obj)
  }

  if (!regime %in% names(obj$regimes)) {
    stop(sprintf("Regime '%s' not found.", regime), call. = FALSE)
  }

  reg <- obj$regimes[[regime]]
  dt <- obj$data
  time_vals <- obj$meta$time_values

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  dt <- .add_tracking_columns(dt, nodes, reg)

  # Fit per censoring source
  for (cvar in nodes$censoring) {
    if (verbose) .vmsg("  Fitting censoring model for '%s'...", cvar)

    results <- vector("list", length(time_vals))

    for (i in seq_along(time_vals)) {
      tt <- time_vals[i]
      dt_t <- dt[dt[[nodes$time]] == tt, ]

      # Risk set: consistent through t-1, uncensored through t-1,
      # AND treatment at t consistent with regime
      if (i == 1) {
        still_in <- rep(TRUE, nrow(dt_t))
      } else {
        still_in <- dt_t$.longy_consist_prev == 1L & dt_t$.longy_uncens_prev == 1L
        still_in[is.na(still_in)] <- FALSE
      }

      # Also require regime-consistent treatment AT time t
      still_in <- still_in & dt_t$.longy_regime_consist == 1L

      n_risk <- sum(still_in)
      if (n_risk == 0) {
        if (verbose) .vmsg("    g_C(%s) time %d: 0 at risk, skipping", cvar, tt)
        next
      }

      X <- as.data.frame(dt_t[still_in, covariates, with = FALSE])
      # Model P(uncensored) = P(C=0) => Y = 1 - C
      Y <- 1L - dt_t[[cvar]][still_in]

      if (length(unique(Y)) > 1 && n_risk >= min_obs) {
        cv_folds <- 10L
        if (adaptive_cv) {
          cv_info <- .adaptive_cv_folds(Y)
          cv_folds <- cv_info$V
        }
        fit <- .safe_sl(Y = Y, X = X, learners = learners,
                        cv_folds = cv_folds, verbose = verbose)
        p_c <- .bound(fit$predictions, bounds[1], bounds[2])
        method <- fit$method
      } else {
        p_c <- rep(mean(Y), n_risk)
        method <- "marginal"
      }

      marg_c <- mean(Y)

      results[[i]] <- data.table::data.table(
        .id = dt_t[[nodes$id]][still_in],
        .time = tt,
        .n_risk = n_risk,
        .censored = dt_t[[cvar]][still_in],
        .marg_c = marg_c,
        .p_c = p_c,
        .method = method
      )

      if (verbose) {
        .vmsg("    g_C(%s) time %d: n_risk=%d, marg_uncens=%.3f, method=%s",
              cvar, tt, n_risk, marg_c, method)
      }
    }

    results <- data.table::rbindlist(
      results[!vapply(results, is.null, logical(1))]
    )
    data.table::setnames(results, ".id", nodes$id)

    obj$fits$censoring[[cvar]] <- list(
      predictions = results,
      covariates = covariates,
      learners = learners,
      bounds = bounds
    )
  }

  .remove_tracking_columns(obj$data)

  obj
}
