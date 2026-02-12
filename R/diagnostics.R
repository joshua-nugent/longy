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
      mean_weight = mean(.final_weight),
      median_weight = stats::median(.final_weight),
      sd_weight = stats::sd(.final_weight),
      min_weight = min(.final_weight),
      max_weight = max(.final_weight),
      p01 = stats::quantile(.final_weight, 0.01),
      p05 = stats::quantile(.final_weight, 0.05),
      p95 = stats::quantile(.final_weight, 0.95),
      p99 = stats::quantile(.final_weight, 0.99),
      ess = .ess(.final_weight),
      mean_sw_a = mean(.sw_a),
      mean_sw_c = mean(.sw_c),
      mean_csw_ac = mean(.csw_ac),
      mean_sw_r = mean(.sw_r)
    ), by = .time]
    data.table::setnames(diag, ".time", "time")
    data.table::setkey(diag, time)
  } else {
    diag <- w_dt[, list(
      n = .N,
      mean_weight = mean(.final_weight),
      median_weight = stats::median(.final_weight),
      sd_weight = stats::sd(.final_weight),
      min_weight = min(.final_weight),
      max_weight = max(.final_weight),
      p01 = stats::quantile(.final_weight, 0.01),
      p99 = stats::quantile(.final_weight, 0.99),
      ess = .ess(.final_weight)
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
