#' Compute Unadjusted (Crude) Estimates
#'
#' Computes naive means of the outcome among uncensored, observed subjects
#' who naturally followed the specified regime through each time point.
#' No models, weights, or causal adjustment are applied. This serves as a
#' reference to show how much causal adjustment matters.
#'
#' Because no nuisance models are needed, this estimator can be called
#' directly after \code{\link{define_regime}} without any \code{fit_*} steps.
#'
#' @param obj A \code{longy_data} object with at least one regime defined.
#' @param regime Character. Name(s) of the regime(s). NULL = all defined.
#' @param times Numeric vector. Time points at which to estimate. NULL = all.
#' @param ci_level Numeric. Confidence level (default 0.95).
#'
#' @return Modified \code{longy_data} object with unadjusted results stored in
#'   \code{obj$results} (keyed as \code{{regime}_unadjusted}).
#'
#' @export
estimate_unadjusted <- function(obj, regime = NULL, times = NULL,
                                ci_level = 0.95) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)

  for (rname in regime) {

  nodes <- obj$nodes
  id_col <- nodes$id
  time_col <- nodes$time
  y_col <- nodes$outcome
  r_col <- nodes$observation

  # Add tracking columns (regime consistency + censoring status)
  .add_tracking_columns(obj$data, nodes, obj$regimes[[rname]])

  # Determine time points
  available_times <- sort(unique(obj$data[[time_col]]))
  eval_times <- if (is.null(times)) available_times else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s not in data, removing.",
                      paste(bad_times, collapse = ", ")))
    }
    intersect(times, available_times)
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)

  est_list <- vector("list", length(eval_times))
  for (k in seq_along(eval_times)) {
    tt <- eval_times[k]
    dt_t <- obj$data[obj$data[[time_col]] == tt, ]

    # Subset: naturally followed regime AND uncensored through t
    keep <- dt_t$.longy_cum_consist == 1L & dt_t$.longy_cum_uncens == 1L

    # Also require observed at t if observation node exists
    if (!is.null(r_col)) {
      keep <- keep & dt_t[[r_col]] == 1L
    }

    # Require non-NA outcome
    yi <- dt_t[[y_col]][keep]
    yi <- yi[!is.na(yi)]
    n <- length(yi)

    if (n > 0) {
      psi_hat <- mean(yi)
      # SE: binomial/survival -> sqrt(p(1-p)/n), continuous -> sd/sqrt(n)
      if (nodes$outcome_type %in% c("binary", "survival")) {
        se <- sqrt(psi_hat * (1 - psi_hat) / n)
      } else {
        se <- stats::sd(yi) / sqrt(n)
      }
    } else {
      psi_hat <- NA_real_
      se <- NA_real_
    }

    est_list[[k]] <- data.table::data.table(
      time = tt,
      estimate = psi_hat,
      se = se,
      ci_lower = psi_hat - z * se,
      ci_upper = psi_hat + z * se,
      n_at_risk = n
    )
  }
  estimates <- data.table::rbindlist(est_list)

  # Isotonic smoothing for survival outcomes (enforce monotone non-decreasing)
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    estimates$ci_lower <- estimates$estimate - z * estimates$se
    estimates$ci_upper <- estimates$estimate + z * estimates$se
  }

  # Clean up tracking columns
  .remove_tracking_columns(obj$data)

  result <- list(
    estimates = estimates,
    regime = rname,
    estimator = "unadjusted",
    inference = "analytic",
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_unadjusted")]] <- result

  } # end for (rname in regime)

  obj
}
