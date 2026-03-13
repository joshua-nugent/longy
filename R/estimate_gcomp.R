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
estimate_gcomp <- function(obj, regime = NULL, times = NULL,
                           ci_level = 0.95, n_boot = 200L,
                           verbose = TRUE) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)

  if (ci_level <= 0 || ci_level >= 1)
    stop("ci_level must be between 0 and 1.", call. = FALSE)

  if (isTRUE(obj$crossfit$enabled))
    warning("Cross-fitting is enabled but G-comp does not currently support ",
            "cross-fitted estimation. Using full-sample outcome fit.",
            call. = FALSE)

  # Preserve user's original times so each regime gets the same input

  user_times <- times

  for (rname in regime) {

  # Reset times for each regime to avoid cross-contamination
  times <- user_times
  if (!is.null(times)) times <- sort(unique(times))

  if (is.null(obj$fits$outcome[[rname]]) || length(obj$fits$outcome[[rname]]) == 0) {
    stop(sprintf("Outcome model not fitted for regime '%s'. Run fit_outcome() first.", rname),
         call. = FALSE)
  }

  if (isTRUE(obj$fits$outcome[[rname]]$metadata_only)) {
    stop(sprintf(
      "Outcome model for regime '%s' contains metadata only (no predictions). ",
      rname),
      "Re-run fit_outcome() with metadata_only=FALSE, or call longy() with ",
      "estimator including 'gcomp'.",
      call. = FALSE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  pred_dt <- obj$fits$outcome[[rname]]$predictions

  # Determine time points
  available_times <- sort(unique(pred_dt$.target_time))
  if (is.null(times)) {
    times <- available_times
  } else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s have no outcome predictions, removing.",
                      paste(bad_times, collapse = ", ")))
      times <- intersect(times, available_times)
    }
    if (length(times) == 0)
      stop(sprintf("No valid time points for regime '%s'. Available: %s",
                   rname, paste(available_times, collapse = ", ")),
           call. = FALSE)
  }

  # Point estimates: mean(Q_hat(t)) for each time
  est_list <- vector("list", length(times))
  for (k in seq_along(times)) {
    tt <- times[k]
    q_t <- pred_dt[pred_dt$.target_time == tt, ]
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
      obj, regime = rname, times = times,
      n_boot = n_boot, ci_level = ci_level, verbose = verbose
    )
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival outcomes
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    # Shift CIs to match smoothed estimates (preserves bootstrap CI width)
    if ("ci_lower" %in% names(estimates)) {
      shift_iso <- estimates$estimate - raw_est
      estimates$ci_lower <- estimates$ci_lower + shift_iso
      estimates$ci_upper <- estimates$ci_upper + shift_iso
    }
  }

  # Clamp CIs to [0, 1] for binary/survival outcomes
  if (nodes$outcome_type %in% c("binary", "survival") &&
      "ci_lower" %in% names(estimates)) {
    estimates$ci_lower <- pmax(estimates$ci_lower, 0)
    estimates$ci_upper <- pmin(estimates$ci_upper, 1)
  }

  result <- list(
    estimates = estimates,
    regime = rname,
    estimator = "gcomp",
    inference = if (n_boot > 0) "bootstrap" else "none",
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_gcomp")]] <- result

  } # end for (rname in regime)

  obj
}
