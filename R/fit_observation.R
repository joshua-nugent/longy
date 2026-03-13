#' Fit Observation Models (g_R)
#'
#' Models P(R(t) = 1 | past) -- the probability that the outcome is measured at
#' time t. This handles intermittent missingness, which is distinct from
#' absorbing censoring.
#'
#' The risk set is all subjects uncensored through t (not regime-concordant).
#' Models are fit on the full uncensored sample.
#'
#' longy assumes the within-period ordering L(t) -> A(t) -> C(t) -> Y(t), so
#' current-time treatment A(t) is included as a default covariate. Override via
#' the \code{covariates} argument.
#'
#' @param obj A `longy_data` object.
#' @param regime Character. Name of the regime.
#' @param covariates Character vector. Predictor columns.
#' @param learners Character vector. SuperLearner library.
#' @param sl_control List. Additional arguments passed to SuperLearner.
#'   Elements named \code{cvControl} are merged with the default
#'   \code{cvControl}. \code{cvControl$V} sets the number of CV folds when
#'   \code{adaptive_cv = FALSE}; specifying \code{V} with
#'   \code{adaptive_cv = TRUE} is an error.
#' @param adaptive_cv Logical. Adaptive CV fold selection. When TRUE (default),
#'   CV folds are chosen automatically. When FALSE, uses
#'   \code{sl_control$cvControl$V} if specified, or 10.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param min_events Integer. Minimum minority-class events required to fit a
#'   model. When the minority class count is below this AND the minority rate
#'   is below 0.01, marginal fallback is used. Default 20.
#' @param bounds Numeric(2). Prediction bounds.
#' @param times Numeric vector. If provided, only fit through `max(times)`.
#' @param use_ffSL Logical. If TRUE, use future-factorial SuperLearner.
#'   Default FALSE. Forced to FALSE inside parallel workers.
#' @param parallel Logical. If TRUE and a non-sequential \code{future::plan()}
#'   is active, dispatches time-point models in parallel. Default FALSE.
#' @param verbose Logical. Progress messages.
#' @param refit Logical. If FALSE (default), errors when observation is already
#'   fitted for the requested regime(s). Set to TRUE to re-fit.
#'
#' @return Modified `longy_data` object with observation fits stored.
#' @export
fit_observation <- function(obj, regime = NULL, covariates = NULL, learners = NULL,
                            sl_control = list(), adaptive_cv = TRUE,
                            min_obs = 50L, min_events = 20L,
                            bounds = c(0.005, 0.995),
                            times = NULL, use_ffSL = FALSE,
                            parallel = FALSE,
                            verbose = TRUE, refit = FALSE) {
  obj <- .as_longy_data(obj)
  learners <- .resolve_learners(learners, "observation")
  regime <- .resolve_regimes(obj, regime)

  if (adaptive_cv && !is.null(sl_control$cvControl$V))
    stop("Cannot specify sl_control$cvControl$V when adaptive_cv=TRUE. ",
         "Set adaptive_cv=FALSE to use a fixed number of CV folds.",
         call. = FALSE)

  if (!refit) {
    fitted <- Filter(function(r) !is.null(obj$fits$observation[[r]]) &&
                       length(obj$fits$observation[[r]]) > 0, regime)
    if (length(fitted) > 0)
      stop(sprintf("Observation already fitted for: %s. Use refit=TRUE to override.",
                   paste(fitted, collapse = ", ")), call. = FALSE)
  }

  nodes <- obj$nodes

  # If no observation column, skip
  if (is.null(nodes$observation)) {
    if (verbose) .vmsg("  No observation column defined, skipping fit_observation()")
    return(obj)
  }

  # g_R models are fit on observed data regardless of regime (risk set uses
  # uncensored status only, not regime consistency). Fit once, then replicate
  # the result for all requested regimes.
  if (isTRUE(obj$crossfit$enabled)) {
    sl_fn_cf <- if (use_ffSL) "ffSL" else "SuperLearner"
    obj <- .cf_fit_observation(obj, regime[1], covariates = covariates,
                                learners = learners, sl_control = sl_control,
                                adaptive_cv = adaptive_cv, min_obs = min_obs,
                                min_events = min_events,
                                bounds = bounds, times = times, sl_fn = sl_fn_cf,
                                verbose = verbose)
    fit_result <- obj$fits$observation[[regime[1]]]
  } else {

  reg <- obj$regimes[[regime[1]]]
  dt <- obj$data
  time_vals <- obj$meta$time_values
  if (!is.null(times)) {
    time_vals <- time_vals[time_vals <= max(times)]
  }

  if (is.null(covariates)) {
    covariates <- c(nodes$baseline, nodes$timevarying, nodes$treatment)
  }

  dt <- .add_tracking_columns(dt, nodes, reg)
  on.exit(.remove_tracking_columns(obj$data), add = TRUE)

  # Snapshot data for parallel safety
  if (parallel) dt <- data.table::copy(dt)

  # Worker SL flag: force sequential SL inside parallel workers
  worker_ffSL <- if (parallel) FALSE else use_ffSL

  one_timepoint <- function(i) {
    tt <- time_vals[i]
    dt_t <- dt[dt[[nodes$time]] == tt, ]

    # Risk set: uncensored through t-1 (full sample, not regime-concordant)
    if (i == 1) {
      still_in <- rep(TRUE, nrow(dt_t))
    } else {
      still_in <- dt_t$.longy_uncens_prev == 1L
      still_in[is.na(still_in)] <- FALSE
    }

    # Must be uncensored at t (observation requires being present)
    still_in <- still_in & dt_t$.longy_uncens == 1L
    still_in[is.na(still_in)] <- FALSE

    # Exclude subjects with NA observation indicator (absorbed in survival data)
    still_in <- still_in & !is.na(dt_t[[nodes$observation]])

    n_risk <- sum(still_in)
    if (n_risk == 0) {
      return(list(result = NULL, sl_info = NULL, n_marginal = 0L))
    }

    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(covariates, lag_covs)
    X <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
    Y <- dt_t[[nodes$observation]][still_in]

    # Extract sampling weights for at-risk subjects (NULL if none)
    ow <- NULL
    if (!is.null(nodes$sampling_weights)) {
      ow <- dt_t[[nodes$sampling_weights]][still_in]
    }

    task_n_marginal <- 0L
    n_minority <- min(sum(Y == 1), sum(Y == 0))
    minority_rate <- min(mean(Y), 1 - mean(Y))
    rare_events <- n_minority < min_events && minority_rate < 0.01
    if (length(unique(Y)) > 1 && n_risk >= min_obs && !rare_events) {
      obs_rate <- mean(Y)
      ctx <- sprintf("g_R, time=%d, n=%d, obs_rate=%.3f", tt, n_risk, obs_rate)
      cv_folds <- if (!is.null(sl_control$cvControl$V)) sl_control$cvControl$V else 10L
      if (adaptive_cv) {
        cv_info <- .adaptive_cv_folds(Y)
        cv_folds <- cv_info$V
      }
      fit <- .safe_sl(Y = Y, X = X, learners = learners,
                      cv_folds = cv_folds, obs_weights = ow,
                      sl_control = sl_control,
                      use_ffSL = worker_ffSL, context = ctx,
                      verbose = !parallel && verbose)
      p_r <- .bound(fit$predictions, bounds[1], bounds[2])
      method <- fit$method
      sl_risk <- fit$sl_risk
      sl_coef <- fit$sl_coef
    } else {
      marg <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)
      p_r <- .bound(rep(marg, n_risk), bounds[1], bounds[2])
      method <- "marginal"
      sl_risk <- NULL
      sl_coef <- NULL
      task_n_marginal <- 1L
      if (length(unique(Y)) <= 1) {
        if (mean(Y) == 1) {
          warning(sprintf(
            "g_R at time %d: all subjects observed (obs_rate=1.000). Using marginal.",
            tt), call. = FALSE)
        } else {
          warning(sprintf(
            "g_R at time %d: no subjects observed (obs_rate=0.000). Using marginal.",
            tt), call. = FALSE)
        }
      } else if (rare_events) {
        warning(sprintf(
          "g_R at time %d: only %d minority-class events (rate=%.3f, min_events=%d). Using marginal.",
          tt, n_minority, minority_rate, min_events), call. = FALSE)
      } else {
        warning(sprintf(
          "g_R at time %d: only %d at risk (min_obs=%d). Using marginal (=%.3f). Consider reducing min_obs.",
          tt, n_risk, min_obs, marg), call. = FALSE)
      }
    }

    marg_r <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)

    result_dt <- data.table::data.table(
      .id = dt_t[[nodes$id]][still_in],
      .time = tt,
      .n_risk = n_risk,
      .observed = Y,
      .marg_r = marg_r,
      .p_r = p_r,
      .method = method
    )

    sl_info_entry <- list(time = tt, method = method,
                          sl_risk = sl_risk, sl_coef = sl_coef)

    if (!parallel && verbose) {
      .vmsg("  g_R time %d: n_risk=%d, marg=%.3f, method=%s",
            tt, n_risk, marg_r, method)
      .vmsg_covariates(nodes$baseline, nodes$timevarying, lag_covs)
    }

    list(result = result_dt, sl_info = sl_info_entry,
         n_marginal = task_n_marginal)
  }

  # Strip obj from closure to avoid serializing the full longy_data object
  if (parallel) {
    one_timepoint <- .clean_closure(one_timepoint, c(
      "time_vals", "dt", "nodes", "parallel",
      "learners", "adaptive_cv", "worker_ffSL", "verbose",
      "bounds", "min_obs", "min_events", "covariates", "sl_control"
    ))
  }

  if (verbose && parallel)
    .vmsg("  g_R: dispatching %d time points...", length(time_vals))

  task_results <- .parallel_or_sequential(
    seq_along(time_vals), one_timepoint, parallel = parallel,
    verbose = FALSE
  )

  # Collect results
  results <- lapply(task_results, `[[`, "result")
  sl_info <- lapply(task_results, `[[`, "sl_info")
  n_marginal <- sum(vapply(task_results, function(x) x$n_marginal, integer(1)))

  non_null <- !vapply(results, is.null, logical(1))
  n_fitted <- sum(non_null)
  if (!any(non_null)) {
    warning("No observations at risk for any time point in observation (g_R) model.",
            call. = FALSE)
    results <- data.table::data.table()
  } else {
    if (n_marginal > 0 && n_marginal >= n_fitted * 0.5) {
      warning(sprintf(
        "g_R: marginal fallback used at %d/%d time points. Model may be unreliable.",
        n_marginal, n_fitted), call. = FALSE)
    }
    results <- data.table::rbindlist(results[non_null])
    data.table::setnames(results, ".id", nodes$id)
  }

  sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

  fit_result <- list(
    predictions = results,
    covariates = covariates,
    learners = learners,
    bounds = bounds,
    use_ffSL = use_ffSL,
    sl_control = sl_control,
    adaptive_cv = adaptive_cv,
    sl_info = sl_info
  )

  .remove_tracking_columns(obj$data)

  } # end if/else crossfit

  # Store for all regimes (identical model, only regime label differs)
  for (rname in regime) {
    obj$fits$observation[[rname]] <- fit_result
    obj$fits$observation[[rname]]$regime <- rname
  }

  obj
}
