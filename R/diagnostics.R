#' Weight Diagnostics
#'
#' Summarizes weight distributions by time point, including mean, median,
#' maximum, quantiles, and effective sample size.
#'
#' @param obj A `longy_data` object with weights computed.
#' @param by_time Logical. If TRUE (default), summarize by time point.
#'
#' @return A data.table with weight summary statistics.
#' @export
weight_diagnostics <- function(obj, by_time = TRUE) {
  stopifnot(inherits(obj, "longy_data"))
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
#' @param obj A `longy_data` object with treatment model fit.
#' @param threshold Numeric. Flag observations with predicted probabilities
#'   below this value (or above 1 - threshold). Default 0.025.
#'
#' @return A data.table with flagged observations and their propensity scores.
#' @export
positivity_diagnostics <- function(obj, threshold = 0.025) {
  stopifnot(inherits(obj, "longy_data"))
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
