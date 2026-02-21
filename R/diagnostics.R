# Extract longy_data from any longy object type
# Accepts longy_data, longy_result, or longy_results (takes first element)
.extract_longy_data <- function(obj) {
  if (inherits(obj, "longy_data")) return(obj)
  if (inherits(obj, "longy_result")) return(obj$obj)
  if (inherits(obj, "longy_results")) return(obj[[1]]$obj)
  stop("Expected a longy_data, longy_result, or longy_results object.",
       call. = FALSE)
}

#' Weight Diagnostics
#'
#' Summarizes weight distributions by time point, including mean, median,
#' maximum, quantiles, and effective sample size.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with weights computed. For \code{longy_results}, uses the first
#'   element.
#' @param by_time Logical. If TRUE (default), summarize by time point.
#'
#' @return A data.table with weight summary statistics.
#' @export
weight_diagnostics <- function(obj, by_time = TRUE) {
  obj <- .extract_longy_data(obj)
  if (is.null(obj$weights)) {
    stop("Weights not computed. Run compute_weights() first.", call. = FALSE)
  }

  w_dt <- obj$weights$weights_dt

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
#' @param threshold Numeric. Flag observations with predicted probabilities
#'   below this value (or above 1 - threshold). Default 0.025.
#'
#' @return A data.table with flagged observations and their propensity scores.
#' @export
positivity_diagnostics <- function(obj, threshold = 0.025) {
  obj <- .extract_longy_data(obj)
  if (is.null(obj$fits$treatment)) {
    stop("Treatment model not fit. Run fit_treatment() first.", call. = FALSE)
  }

  gA <- obj$fits$treatment$predictions
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
sl_diagnostics <- function(obj, model = "all") {
  obj <- .extract_longy_data(obj)

  model <- match.arg(model, c("all", "treatment", "censoring",
                               "observation", "outcome"))
  include <- if (model == "all") {
    c("treatment", "censoring", "observation", "outcome")
  } else {
    model
  }

  rows <- list()

  # Treatment
  if ("treatment" %in% include && !is.null(obj$fits$treatment$sl_info)) {
    for (entry in obj$fits$treatment$sl_info) {
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
  if ("censoring" %in% include && length(obj$fits$censoring) > 0) {
    for (cvar in names(obj$fits$censoring)) {
      sl_info <- obj$fits$censoring[[cvar]]$sl_info
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
  if ("observation" %in% include && !is.null(obj$fits$observation$sl_info)) {
    for (entry in obj$fits$observation$sl_info) {
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
  if ("outcome" %in% include && !is.null(obj$fits$outcome$sl_info)) {
    for (entry in obj$fits$outcome$sl_info) {
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
