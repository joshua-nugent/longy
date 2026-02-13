#' Estimate Causal Effects via G-Computation
#'
#' Computes the G-computation (sequential regression) estimator at each
#' requested time point, using the fitted outcome models from
#' \code{\link{fit_outcome}}. Inference is via the nonparametric bootstrap.
#'
#' @param obj A \code{longy_data} object with outcome models fitted via
#'   \code{fit_outcome()}.
#' @param regime Character. Name of the regime.
#' @param times Numeric vector. Time points at which to estimate. If NULL,
#'   estimates at all time points with outcome predictions.
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param n_boot Integer. Number of bootstrap replicates (default 200).
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class \code{"longy_result"} with elements:
#'   \describe{
#'     \item{estimates}{data.table with time, estimate, n_at_risk}
#'     \item{regime}{Name of the regime}
#'     \item{inference}{Always "bootstrap" for G-comp}
#'     \item{ci_level}{Confidence level}
#'   }
#'
#' @export
estimate_gcomp <- function(obj, regime, times = NULL,
                           ci_level = 0.95, n_boot = 200L,
                           verbose = TRUE) {
  stopifnot(inherits(obj, "longy_data"))

  if (is.null(obj$fits$outcome)) {
    stop("Outcome model not fitted. Run fit_outcome() first.", call. = FALSE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  pred_dt <- obj$fits$outcome$predictions

  # Determine time points
  available_times <- sort(unique(pred_dt$.time))
  if (is.null(times)) {
    times <- available_times
  } else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s have no outcome predictions, removing.",
                      paste(bad_times, collapse = ", ")))
      times <- intersect(times, available_times)
    }
  }

  # Point estimates: mean(Q_hat(t)) for each time
  est_list <- vector("list", length(times))
  for (k in seq_along(times)) {
    tt <- times[k]
    q_t <- pred_dt[pred_dt$.time == tt, ]
    n_at_risk <- nrow(q_t)
    psi_hat <- if (n_at_risk > 0) mean(q_t$.Q_hat) else NA_real_

    est_list[[k]] <- data.table::data.table(
      time = tt,
      estimate = psi_hat,
      n_effective = as.numeric(n_at_risk),
      n_at_risk = n_at_risk
    )
  }
  estimates <- data.table::rbindlist(est_list)

  # Bootstrap inference
  if (n_boot > 0) {
    inf_dt <- .bootstrap_gcomp_inference(
      obj, regime = regime, times = times,
      n_boot = n_boot, ci_level = ci_level, verbose = verbose
    )
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival outcomes
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    if ("se" %in% names(estimates)) {
      z <- stats::qnorm(1 - (1 - ci_level) / 2)
      estimates$ci_lower <- estimates$estimate - z * estimates$se
      estimates$ci_upper <- estimates$estimate + z * estimates$se
    }
  }

  result <- list(
    estimates = estimates,
    regime = regime,
    estimator = "gcomp",
    inference = "bootstrap",
    ci_level = ci_level,
    obj = obj
  )
  class(result) <- "longy_result"
  result
}
