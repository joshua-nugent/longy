#' Fit Treatment Models (g_A)
#'
#' Fits propensity score models P(A(t) = 1 | past) at each time point for
#' subjects in the risk set (uncensored through t-1). Models are fit on the
#' full uncensored sample regardless of regime consistency, matching ltmle/stremr.
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
#' @param min_events Integer. Minimum minority-class events required to fit a
#'   model. When the minority class count is below this AND the minority rate
#'   is below 0.01, marginal fallback is used. Default 20.
#' @param bounds Numeric vector of length 2. Bounds for predicted probabilities.
#' @param times Numeric vector. If provided, only fit models through
#'   `max(times)`. Saves computation when estimation is only needed at
#'   early time points.
#' @param sl_fn Character. SuperLearner implementation: \code{"SuperLearner"}
#'   (default) or \code{"ffSL"} (future-factorial parallel).
#' @param verbose Logical. Print progress.
#' @param refit Logical. If FALSE (default), errors when treatment is already
#'   fitted for the requested regime(s). Set to TRUE to re-fit.
#' @param risk_set Character. Which subjects form the risk set for treatment
#'   model fitting. \code{"all"} (default) uses all uncensored subjects
#'   (fit once, share across regimes). \code{"followers"} restricts to
#'   regime-followers at each time, fitting separate models per regime.
#'
#' @return Modified `longy_data` object with treatment fits stored.
#' @export
fit_treatment <- function(obj, regime = NULL, covariates = NULL, learners = NULL,
                          sl_control = list(), adaptive_cv = TRUE,
                          min_obs = 50L, min_events = 20L,
                          bounds = c(0.005, 0.995),
                          times = NULL, sl_fn = "SuperLearner",
                          verbose = TRUE, refit = FALSE,
                          risk_set = c("all", "followers")) {
  obj <- .as_longy_data(obj)
  learners <- .resolve_learners(learners, "treatment")
  regime <- .resolve_regimes(obj, regime)
  risk_set <- match.arg(risk_set)

  if (!refit) {
    fitted <- Filter(function(r) !is.null(obj$fits$treatment[[r]]) &&
                       length(obj$fits$treatment[[r]]) > 0, regime)
    if (length(fitted) > 0)
      stop(sprintf("Treatment already fitted for: %s. Use refit=TRUE to override.",
                   paste(fitted, collapse = ", ")), call. = FALSE)
  }

  # When risk_set = "all", fit once on regime[1] and copy to all regimes.
  # When risk_set = "followers", fit separately per regime (each has
  # different followers, so models differ).
  fit_regimes <- if (risk_set == "followers") regime else regime[1]

  for (fit_rname in fit_regimes) {

  if (isTRUE(obj$crossfit$enabled)) {
    obj <- .cf_fit_treatment(obj, fit_rname, covariates = covariates,
                              learners = learners, sl_control = sl_control,
                              adaptive_cv = adaptive_cv, min_obs = min_obs,
                              min_events = min_events,
                              bounds = bounds, times = times, sl_fn = sl_fn,
                              verbose = verbose, risk_set = risk_set)
    fit_result <- obj$fits$treatment[[fit_rname]]
  } else {

  reg <- obj$regimes[[fit_rname]]
  dt <- obj$data
  nodes <- obj$nodes
  time_vals <- obj$meta$time_values
  if (!is.null(times)) {
    time_vals <- time_vals[time_vals <= max(times)]
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying)
  }

  # Build cumulative regime-consistency and uncensored tracking
  dt <- .add_tracking_columns(dt, nodes, reg)

  results <- vector("list", length(time_vals))
  sl_info <- vector("list", length(time_vals))
  n_marginal <- 0L

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
      if (verbose) .vmsg("  g_A time %d: 0 at risk, skipping", tt)
      next
    }

    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(covariates, lag_covs)
    X <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
    Y <- dt_t[[nodes$treatment]][still_in]

    # Extract sampling weights for at-risk subjects (NULL if none)
    ow <- NULL
    if (!is.null(nodes$sampling_weights)) {
      ow <- dt_t[[nodes$sampling_weights]][still_in]
    }

    # Fit model or use marginal
    n_minority <- min(sum(Y == 1), sum(Y == 0))
    minority_rate <- min(mean(Y), 1 - mean(Y))
    rare_events <- n_minority < min_events && minority_rate < 0.01
    if (length(unique(Y)) > 1 && n_risk >= min_obs && !rare_events) {
      cv_folds <- 10L
      if (adaptive_cv) {
        cv_info <- .adaptive_cv_folds(Y)
        cv_folds <- cv_info$V
      }
      ctx <- sprintf("g_A, time=%d, n=%d, event_rate=%.3f", tt, n_risk, mean(Y))
      fit <- .safe_sl(Y = Y, X = X, learners = learners,
                      cv_folds = cv_folds, obs_weights = ow,
                      sl_fn = sl_fn, context = ctx, verbose = verbose)
      p_a <- .bound(fit$predictions, bounds[1], bounds[2])
      method <- fit$method
      sl_risk <- fit$sl_risk
      sl_coef <- fit$sl_coef
    } else {
      marg <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)
      p_a <- .bound(rep(marg, n_risk), bounds[1], bounds[2])
      method <- "marginal"
      sl_risk <- NULL
      sl_coef <- NULL
      n_marginal <- n_marginal + 1L
      # Warn with reason for marginal fallback
      if (length(unique(Y)) <= 1) {
        warning(sprintf(
          "g_A at time %d: outcome is constant (all %s). Using marginal. Check treatment variation.",
          tt, if (all(Y == 1)) "1" else "0"), call. = FALSE)
      } else if (rare_events) {
        warning(sprintf(
          "g_A at time %d: only %d minority-class events (rate=%.3f, min_events=%d). Using marginal.",
          tt, n_minority, minority_rate, min_events), call. = FALSE)
      } else {
        warning(sprintf(
          "g_A at time %d: only %d at risk (min_obs=%d). Using marginal (=%.3f). Consider reducing min_obs.",
          tt, n_risk, min_obs, marg), call. = FALSE)
      }
    }

    marg_a <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)

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
      rs_label <- if (risk_set == "followers") {
        sprintf(" (%s followers)", fit_rname)
      } else {
        ""
      }
      .vmsg("  g_A time %d%s: n_risk=%d, marg=%.3f, method=%s",
            tt, rs_label, n_risk, marg_a, method)
    }
  }

  non_null <- !vapply(results, is.null, logical(1))
  n_fitted <- sum(non_null)
  if (!any(non_null)) {
    warning("No observations at risk for any time point in treatment (g_A) model.",
            call. = FALSE)
  } else if (n_marginal > 0 && n_marginal >= n_fitted * 0.5) {
    warning(sprintf(
      "g_A: marginal fallback used at %d/%d time points. Model may be unreliable.",
      n_marginal, n_fitted), call. = FALSE)
  }
  results <- data.table::rbindlist(results[non_null])
  data.table::setnames(results, ".id", nodes$id)

  sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

  fit_result <- list(
    predictions = results,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    sl_fn = sl_fn,
    sl_info = sl_info,
    risk_set = risk_set
  )

  # Clean up tracking columns
  .remove_tracking_columns(obj$data)

  } # end if/else crossfit

  # Store fit result
  if (risk_set == "all") {
    for (rname in regime) {
      obj$fits$treatment[[rname]] <- fit_result
      obj$fits$treatment[[rname]]$regime <- rname
    }
  } else {
    obj$fits$treatment[[fit_rname]] <- fit_result
    obj$fits$treatment[[fit_rname]]$regime <- fit_rname
  }

  } # end for fit_rname

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
#'
#' Preserves \code{.longy_fold} (cross-fitting fold assignment) and
#' \code{.longy_lag_*} (covariate history lag columns).
#' @noRd
.remove_tracking_columns <- function(dt) {
  track_cols <- grep("^\\.longy_", names(dt), value = TRUE)
  track_cols <- setdiff(track_cols, ".longy_fold")
  track_cols <- track_cols[!grepl("^\\.longy_lag_", track_cols)]
  if (length(track_cols) > 0) {
    dt[, (track_cols) := NULL]
  }
  invisible(dt)
}
