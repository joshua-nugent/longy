#' Fit Censoring Models (g_C)
#'
#' Fits models for P(C(t) = 0 | past) at each time point, separately for
#' each censoring cause. The risk set is all subjects uncensored through t-1.
#' Models are fit on the full uncensored sample regardless of regime consistency.
#'
#' longy assumes the within-period ordering L(t) -> A(t) -> C(t) -> Y(t), so
#' current-time treatment A(t) is included as a default covariate. This matches
#' the ltmle vignette ordering. Override via the \code{covariates} argument.
#'
#' This function loops over the internal binary censoring columns stored in
#' \code{obj$nodes$censoring} (e.g. \code{".cens_censored"}, \code{".cens_death"}).
#' These are created automatically by \code{\link{longy_data}} from the
#' user-provided character censoring column. If there are no censoring columns
#' (censoring is NULL), the function returns immediately.
#'
#' @param obj A \code{longy_data} object with at least one regime defined.
#' @param regime Character. Name of the regime.
#' @param covariates Character vector. Predictor columns. If NULL, uses all
#'   baseline + timevarying covariates.
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
#' @param times Numeric vector. If provided, only fit through \code{max(times)}.
#' @param use_ffSL Logical. If TRUE, use future-factorial SuperLearner.
#'   Default FALSE. Forced to FALSE inside parallel workers.
#' @param parallel Logical. If TRUE and a non-sequential \code{future::plan()}
#'   is active, dispatches cause x time tasks in parallel. Default FALSE.
#' @param verbose Logical. Progress messages.
#' @param refit Logical. If FALSE (default), errors when censoring is already
#'   fitted for the requested regime(s). Set to TRUE to re-fit.
#'
#' @return Modified \code{longy_data} object with censoring fits stored in
#'   \code{obj$fits$censoring}, a named list keyed by internal column name
#'   (e.g. \code{".cens_censored"}).
#' @export
fit_censoring <- function(obj, regime = NULL, covariates = NULL, learners = NULL,
                          sl_control = list(), adaptive_cv = TRUE,
                          min_obs = 50L, min_events = 20L,
                          bounds = c(0.005, 0.995),
                          times = NULL, use_ffSL = FALSE,
                          parallel = FALSE,
                          verbose = TRUE, refit = FALSE) {
  obj <- .as_longy_data(obj)
  learners <- .resolve_learners(learners, "censoring")
  regime <- .resolve_regimes(obj, regime)

  if (adaptive_cv && !is.null(sl_control$cvControl$V))
    stop("Cannot specify sl_control$cvControl$V when adaptive_cv=TRUE. ",
         "Set adaptive_cv=FALSE to use a fixed number of CV folds.",
         call. = FALSE)

  if (!refit) {
    fitted <- Filter(function(r) !is.null(obj$fits$censoring[[r]]) &&
                       length(obj$fits$censoring[[r]]) > 0, regime)
    if (length(fitted) > 0)
      stop(sprintf("Censoring already fitted for: %s. Use refit=TRUE to override.",
                   paste(fitted, collapse = ", ")), call. = FALSE)
  }

  nodes <- obj$nodes

  # If no censoring columns, nothing to do
  if (is.null(nodes$censoring) || length(nodes$censoring) == 0) {
    if (verbose) .vmsg("  No censoring columns defined, skipping fit_censoring()")
    return(obj)
  }

  # g_C models are fit on observed data regardless of regime (risk set uses
  # uncensored status only, not regime consistency). Fit once, then replicate
  # the result for all requested regimes.
  if (isTRUE(obj$crossfit$enabled)) {
    sl_fn_cf <- if (use_ffSL) "ffSL" else "SuperLearner"
    obj <- .cf_fit_censoring(obj, regime[1], covariates = covariates,
                              learners = learners, sl_control = sl_control,
                              adaptive_cv = adaptive_cv, min_obs = min_obs,
                              min_events = min_events,
                              bounds = bounds, times = times, sl_fn = sl_fn_cf,
                              verbose = verbose)
    fit_result <- obj$fits$censoring[[regime[1]]]
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

  # Flatten cause x time into a single task list for better load balancing
  tasks <- expand.grid(cvar_idx = seq_along(nodes$censoring),
                       time_idx = seq_along(time_vals),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  one_cens_task <- function(task_row) {
    cvar <- nodes$censoring[tasks$cvar_idx[task_row]]
    i <- tasks$time_idx[task_row]
    tt <- time_vals[i]
    dt_t <- dt[dt[[nodes$time]] == tt, ]

    # Risk set: uncensored through t-1 (full sample, not regime-concordant)
    if (i == 1) {
      still_in <- rep(TRUE, nrow(dt_t))
    } else {
      still_in <- dt_t$.longy_uncens_prev == 1L
      still_in[is.na(still_in)] <- FALSE
    }

    # Exclude subjects with NA censoring indicator (absorbed in survival data)
    still_in <- still_in & !is.na(dt_t[[cvar]])

    n_risk <- sum(still_in)
    if (n_risk == 0) {
      return(list(cvar = cvar, result = NULL, sl_info = NULL, n_marginal = 0L,
                  marginal_reason = NULL))
    }

    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(covariates, lag_covs)
    X <- as.data.frame(dt_t[still_in, all_covs, with = FALSE])
    # Model P(uncensored) = P(C=0) => Y = 1 - C
    Y <- 1L - dt_t[[cvar]][still_in]

    # Extract sampling weights for at-risk subjects (NULL if none)
    ow <- NULL
    if (!is.null(nodes$sampling_weights)) {
      ow <- dt_t[[nodes$sampling_weights]][still_in]
    }

    task_n_marginal <- 0L
    task_marginal_reason <- NULL
    n_minority <- min(sum(Y == 1), sum(Y == 0))
    minority_rate <- min(mean(Y), 1 - mean(Y))
    rare_events <- n_minority < min_events && minority_rate < 0.01
    if (length(unique(Y)) > 1 && n_risk >= min_obs && !rare_events) {
      cens_rate <- 1 - mean(Y)
      ctx <- sprintf("g_C(%s), time=%d, n=%d, censoring_rate=%.3f",
                     cvar, tt, n_risk, cens_rate)
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
      p_c <- .bound(fit$predictions, bounds[1], bounds[2])
      method <- fit$method
      sl_risk <- fit$sl_risk
      sl_coef <- fit$sl_coef
      sl_clip_log <- fit$clip_log
    } else {
      marg <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)
      p_c <- .bound(rep(marg, n_risk), bounds[1], bounds[2])
      method <- "marginal"
      sl_risk <- NULL
      sl_coef <- NULL
      sl_clip_log <- NULL
      task_n_marginal <- 1L
      # Record reason for marginal fallback (summarized after loop)
      if (length(unique(Y)) <= 1) {
        task_marginal_reason <- "constant"
      } else if (rare_events) {
        task_marginal_reason <- "minority-class"
      } else {
        task_marginal_reason <- "low_n"
      }
    }

    marg_c <- if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y)

    result_dt <- data.table::data.table(
      .id = dt_t[[nodes$id]][still_in],
      .time = tt,
      .n_risk = n_risk,
      .censored = dt_t[[cvar]][still_in],
      .marg_c = marg_c,
      .p_c = p_c,
      .method = method
    )

    sl_info_entry <- list(time = tt, method = method,
                          sl_risk = sl_risk, sl_coef = sl_coef,
                          clip_log = sl_clip_log)

    if (!parallel && verbose) {
      .vmsg("    g_C(%s) time %d: n_risk=%d, marg_uncens=%.3f, method=%s",
            cvar, tt, n_risk, marg_c, method)
      .vmsg_covariates(nodes$baseline, nodes$timevarying, lag_covs,
                        indent = "      ")
    }

    list(cvar = cvar, result = result_dt, sl_info = sl_info_entry,
         n_marginal = task_n_marginal,
         marginal_reason = task_marginal_reason)
  }

  # Strip obj from closure to avoid serializing the full longy_data object
  if (parallel) {
    one_cens_task <- .clean_closure(one_cens_task, c(
      "tasks", "nodes", "time_vals", "dt", "parallel",
      "learners", "adaptive_cv", "worker_ffSL", "verbose",
      "bounds", "min_obs", "min_events", "covariates", "sl_control"
    ))
  }

  if (verbose && parallel)
    .vmsg("  g_C: dispatching %d cause x time tasks...", nrow(tasks))

  task_results <- .parallel_or_sequential(
    seq_len(nrow(tasks)), one_cens_task, parallel = parallel,
    verbose = FALSE
  )

  # Regroup by censoring cause
  fit_result <- list()
  for (cvar in nodes$censoring) {
    cvar_tasks <- Filter(function(x) identical(x$cvar, cvar), task_results)
    results <- lapply(cvar_tasks, `[[`, "result")
    sl_info <- lapply(cvar_tasks, `[[`, "sl_info")
    n_marginal <- sum(vapply(cvar_tasks, function(x) x$n_marginal, integer(1)))
    marginal_reasons <- unlist(lapply(cvar_tasks, `[[`, "marginal_reason"))

    non_null <- !vapply(results, is.null, logical(1))
    n_fitted <- sum(non_null)
    if (!any(non_null)) {
      warning(sprintf(
        "No observations at risk for any time point in censoring (g_C) model for '%s'.",
        cvar), call. = FALSE)
      next
    } else if (n_marginal > 0) {
      reason_tbl <- table(marginal_reasons)
      reason_str <- paste(sprintf("%s(%d)", names(reason_tbl), reason_tbl),
                          collapse = ", ")
      warning(sprintf(
        "g_C(%s): marginal fallback at %d/%d time points [%s]. Run sl_diagnostics() for details.",
        cvar, n_marginal, n_fitted, reason_str), call. = FALSE)
    }
    results <- data.table::rbindlist(results[non_null])
    data.table::setnames(results, ".id", nodes$id)

    sl_info <- sl_info[!vapply(sl_info, is.null, logical(1))]

    fit_result[[cvar]] <- list(
      predictions = results,
      covariates = covariates,
      learners = learners,
      bounds = bounds,
      use_ffSL = use_ffSL,
      sl_control = sl_control,
      adaptive_cv = adaptive_cv,
      sl_info = sl_info
    )
  }

  .remove_tracking_columns(obj$data)

  } # end if/else crossfit

  # Store for all regimes (identical model)
  for (rname in regime) {
    obj$fits$censoring[[rname]] <- fit_result
  }

  obj
}
