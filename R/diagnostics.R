#' Weight Diagnostics
#'
#' Summarizes weight distributions by time point, including mean, median,
#' maximum, quantiles, and effective sample size.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with weights computed. For \code{longy_results}, uses the first
#'   element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param by_time Logical. If TRUE (default), summarize by time point.
#'
#' @return A data.table with weight summary statistics.
#' @export
weight_diagnostics <- function(obj, regime = NULL, by_time = TRUE) {
  obj <- .as_longy_data(obj)
  if (length(obj$weights) == 0) {
    stop("Weights not computed. Run compute_weights() first.", call. = FALSE)
  }
  if (is.null(regime)) regime <- names(obj$weights)[1]
  if (is.null(obj$weights[[regime]])) {
    stop(sprintf("No weights found for regime '%s'.", regime), call. = FALSE)
  }

  w_dt <- obj$weights[[regime]]$weights_dt

  if (by_time) {
    diag <- w_dt[, list(
      n = .N,
      mean_weight = round(mean(.final_weight), 2),
      median_weight = round(stats::median(.final_weight), 2),
      sd_weight = round(stats::sd(.final_weight), 2),
      min_weight = round(min(.final_weight), 2),
      max_weight = round(max(.final_weight), 2),
      p01 = round(stats::quantile(.final_weight, 0.01), 2),
      p05 = round(stats::quantile(.final_weight, 0.05), 2),
      p95 = round(stats::quantile(.final_weight, 0.95), 2),
      p99 = round(stats::quantile(.final_weight, 0.99), 2),
      ess = round(.ess(.final_weight), 1),
      cv = round(stats::sd(.final_weight) / mean(.final_weight), 2),
      mean_sw_a = round(mean(.sw_a), 2),
      mean_sw_c = round(mean(.sw_c), 2),
      mean_csw_ac = round(mean(.csw_ac), 2),
      mean_sw_r = round(mean(.sw_r), 2)
    ), by = .time]
    data.table::setnames(diag, ".time", "time")
    data.table::setkey(diag, time)
  } else {
    diag <- w_dt[, list(
      n = .N,
      mean_weight = round(mean(.final_weight), 2),
      median_weight = round(stats::median(.final_weight), 2),
      sd_weight = round(stats::sd(.final_weight), 2),
      min_weight = round(min(.final_weight), 2),
      max_weight = round(max(.final_weight), 2),
      p01 = round(stats::quantile(.final_weight, 0.01), 2),
      p99 = round(stats::quantile(.final_weight, 0.99), 2),
      ess = round(.ess(.final_weight), 1),
      cv = round(stats::sd(.final_weight) / mean(.final_weight), 2)
    )]
  }

  diag
}

#' Positivity Diagnostics
#'
#' Identifies observations with extreme propensity scores that may indicate
#' positivity violations.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with treatment model fit. For \code{longy_results}, uses the first
#'   element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param threshold Numeric. Flag observations with predicted probabilities
#'   below this value (or above 1 - threshold). Default 0.025.
#'
#' @return A data.table with flagged observations and their propensity scores.
#' @export
positivity_diagnostics <- function(obj, regime = NULL, threshold = 0.025) {
  obj <- .as_longy_data(obj)
  if (length(obj$fits$treatment) == 0) {
    stop("Treatment model not fit. Run fit_treatment() first.", call. = FALSE)
  }
  if (is.null(regime)) regime <- names(obj$fits$treatment)[1]
  if (is.null(obj$fits$treatment[[regime]])) {
    stop(sprintf("No treatment fit found for regime '%s'.", regime), call. = FALSE)
  }

  gA <- obj$fits$treatment[[regime]]$predictions
  id_col <- obj$nodes$id

  # Flag extreme propensity scores
  flagged <- gA[gA$.p_a < threshold | gA$.p_a > (1 - threshold), ]

  if (nrow(flagged) == 0) {
    message("No positivity violations detected at threshold ", threshold)
    return(invisible(data.table::data.table()))
  }

  # Summary by time
  summary_dt <- flagged[, list(
    n_flagged = .N,
    min_p_a = min(.p_a),
    max_p_a = max(.p_a),
    mean_p_a = mean(.p_a)
  ), by = .time]
  data.table::setnames(summary_dt, ".time", "time")
  data.table::setkey(summary_dt, time)

  message(sprintf("Positivity: %d observations flagged (threshold = %.3f) across %d time points",
                  nrow(flagged), threshold, nrow(summary_dt)))

  summary_dt
}

#' SuperLearner Diagnostics
#'
#' Extracts SuperLearner fit information from all fitted nuisance models,
#' including method used (SL, glm, marginal), CV risk, and learner coefficients.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with at least one model fit. For \code{longy_results}, uses the
#'   first element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param model Character. Which model(s) to include: \code{"all"} (default),
#'   \code{"treatment"}, \code{"censoring"}, \code{"observation"}, or
#'   \code{"outcome"}.
#'
#' @return A data.table with columns: \code{model} (model type),
#'   \code{submodel} (censoring cause or NA), \code{time} (time point),
#'   \code{method} (fitting method used), \code{sl_risk} (CV risk of the
#'   SuperLearner, or NA), \code{sl_coef} (list-column of learner weights,
#'   or NULL).
#' @export
sl_diagnostics <- function(obj, regime = NULL, model = "all") {
  obj <- .as_longy_data(obj)

  model <- match.arg(model, c("all", "treatment", "censoring",
                               "observation", "outcome"))
  include <- if (model == "all") {
    c("treatment", "censoring", "observation", "outcome")
  } else {
    model
  }

  # Resolve regime: default to first available
  if (is.null(regime)) {
    available <- names(obj$fits$treatment)
    if (length(available) == 0) available <- names(obj$fits$outcome)
    if (length(available) == 0) available <- names(obj$fits$observation)
    regime <- available[1]
  }

  rows <- list()

  # If no regime resolved (no fits at all), skip to empty result
  if (is.na(regime) || is.null(regime)) {
    message("No model fits found.")
    return(invisible(data.table::data.table(
      model = character(), submodel = character(), time = integer(),
      method = character(), sl_risk = numeric(), sl_coef = list()
    )))
  }

  # Treatment
  trt_fit <- obj$fits$treatment[[regime]]
  if ("treatment" %in% include && !is.null(trt_fit$sl_info)) {
    for (entry in trt_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "treatment",
        submodel = NA_character_,
        time = entry$time,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  # Censoring (one sub-list per cause)
  cens_fits <- obj$fits$censoring[[regime]]
  if ("censoring" %in% include && length(cens_fits) > 0) {
    for (cvar in names(cens_fits)) {
      sl_info <- cens_fits[[cvar]]$sl_info
      if (is.null(sl_info)) next
      for (entry in sl_info) {
        rows[[length(rows) + 1L]] <- list(
          model = "censoring",
          submodel = cvar,
          time = entry$time,
          method = entry$method,
          sl_risk = entry$sl_risk %||% NA_real_,
          sl_coef = list(entry$sl_coef)
        )
      }
    }
  }

  # Observation
  obs_fit <- obj$fits$observation[[regime]]
  if ("observation" %in% include && !is.null(obs_fit$sl_info)) {
    for (entry in obs_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "observation",
        submodel = NA_character_,
        time = entry$time,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  # Outcome
  out_fit <- obj$fits$outcome[[regime]]
  if ("outcome" %in% include && !is.null(out_fit$sl_info)) {
    for (entry in out_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "outcome",
        submodel = NA_character_,
        time = entry$time,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  if (length(rows) == 0) {
    message("No model fits found.")
    return(invisible(data.table::data.table(
      model = character(), submodel = character(), time = integer(),
      method = character(), sl_risk = numeric(), sl_coef = list()
    )))
  }

  data.table::rbindlist(rows)
}
