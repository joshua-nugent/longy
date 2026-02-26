#' Estimate Causal Effects via Inverse Probability Weighting
#'
#' Computes the Hajek (self-normalized) IPW estimator at each requested time
#' point, with standard errors and confidence intervals.
#'
#' @param obj A `longy_data` object with treatment (and optionally
#'   censoring/observation) models already fit. Weights are computed
#'   automatically if not already present, or recomputed when \code{g_bounds}
#'   differs from the stored value.
#' @param regime Character. Name of the regime.
#' @param times Numeric vector. Time points at which to estimate. If NULL,
#'   estimates at all time points with data.
#' @param inference Character. Inference method: `"ic"` (influence curve),
#'   `"bootstrap"`, `"sandwich"`, or `"none"` (point estimates only).
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param n_boot Integer. Number of bootstrap replicates (only for `"bootstrap"`).
#' @param cluster Character. Column name for clustered standard errors.
#' @param g_bounds Numeric vector of length 2. Bounds for cumulative g
#'   (denominator). Default \code{c(0.01, 1)}. If different from the stored
#'   weights, weights are recomputed automatically (no model refitting needed).
#'
#' @return Modified \code{longy_data} object with IPW results stored in
#'   \code{obj$results}.
#'
#' @export
estimate_ipw <- function(obj, regime = NULL, times = NULL, inference = "ic",
                         ci_level = 0.95, n_boot = 200L, cluster = NULL,
                         g_bounds = c(0.01, 1)) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)
  inference <- match.arg(inference, c("ic", "bootstrap", "sandwich", "none"))

  for (rname in regime) {

  # Recompute weights if g_bounds changed or weights not yet computed
  stored_w <- obj$weights[[rname]]
  needs_weights <- is.null(stored_w) || length(stored_w) == 0
  if (!needs_weights && !identical(stored_w$g_bounds, g_bounds)) {
    needs_weights <- TRUE
  }
  if (needs_weights) {
    if (is.null(obj$fits$treatment[[rname]])) {
      stop(sprintf(
        "No treatment model fitted for regime '%s'. Run fit_treatment() first.",
        rname), call. = FALSE)
    }
    # Recover weight-computation settings from stored weights, or use defaults
    stab  <- if (!is.null(stored_w$stabilized)) stored_w$stabilized else TRUE
    trunc <- if (!is.null(stored_w$truncation)) stored_w$truncation else NULL
    trunc_q <- if (!is.null(stored_w$truncation_quantile)) stored_w$truncation_quantile else NULL
    obj <- compute_weights(obj, regime = rname, stabilized = stab,
                           truncation = trunc, truncation_quantile = trunc_q,
                           g_bounds = g_bounds, recompute = TRUE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights[[rname]]$weights_dt

  # Determine time points
  available_times <- sort(unique(w_dt$.time))
  if (is.null(times)) {
    times <- available_times
  } else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s have no weight data, removing.",
                      paste(bad_times, collapse = ", ")))
      times <- intersect(times, available_times)
    }
  }

  # Point estimates
  est_list <- vector("list", length(times))
  for (k in seq_along(times)) {
    tt <- times[k]
    w_t <- w_dt[w_dt$.time == tt, ]
    dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
    merged <- merge(w_t, dt_t[, c(id_col, nodes$outcome), with = FALSE],
                    by = id_col)

    yi <- merged[[nodes$outcome]]
    wi <- merged$.final_weight
    n_at_risk <- nrow(merged)
    psi_hat <- if (n_at_risk > 0) stats::weighted.mean(yi, wi) else NA_real_
    n_eff <- if (n_at_risk > 0) .ess(wi) else 0

    est_list[[k]] <- data.table::data.table(
      time = tt,
      estimate = psi_hat,
      n_effective = n_eff,
      n_at_risk = n_at_risk
    )
  }
  estimates <- data.table::rbindlist(est_list)

  # Inference (computed on raw/unsmoothed estimates)
  if (inference == "ic") {
    inf_dt <- .ic_inference(estimates, obj, regime = rname, ci_level = ci_level,
                            cluster = cluster)
    estimates <- cbind(estimates, inf_dt)
  } else if (inference == "bootstrap") {
    inf_dt <- .bootstrap_inference(obj, regime = rname, times = times,
                                   n_boot = n_boot, ci_level = ci_level)
    estimates <- cbind(estimates, inf_dt)
  } else if (inference == "sandwich") {
    inf_dt <- .sandwich_inference(obj, regime = rname, times = times,
                                  ci_level = ci_level, cluster = cluster)
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival outcomes (enforce monotone non-decreasing)
  # Applied after inference so SEs are computed from unsmoothed ICs.
  # CIs are recomputed around the smoothed point estimate.
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
    regime = rname,
    estimator = "ipw",
    inference = inference,
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_ipw")]] <- result

  } # end for (rname in regime)

  obj
}

#' @export
print.longy_result <- function(x, ...) {
  est_label <- if (!is.null(x$estimator) && x$estimator == "gcomp") "G-comp" else
    if (!is.null(x$estimator) && x$estimator == "tmle") "TMLE" else "IPW"
  cat(sprintf("longy %s result -- regime: %s\n", est_label, x$regime))
  cat(sprintf("Inference: %s | CI level: %.0f%%\n\n",
              x$inference, x$ci_level * 100))

  est <- x$estimates
  # Format for display
  disp <- data.frame(
    time = est$time,
    estimate = round(est$estimate, 4),
    se = if ("se" %in% names(est)) round(est$se, 4) else NA,
    ci_lower = if ("ci_lower" %in% names(est)) round(est$ci_lower, 4) else NA,
    ci_upper = if ("ci_upper" %in% names(est)) round(est$ci_upper, 4) else NA,
    n_eff = round(est$n_effective, 1),
    n_risk = est$n_at_risk
  )
  print(disp, row.names = FALSE)
  invisible(x)
}

#' @export
summary.longy_result <- function(object, ...) {
  est_label <- if (!is.null(object$estimator) && object$estimator == "gcomp") "G-comp" else
    if (!is.null(object$estimator) && object$estimator == "tmle") "TMLE" else "IPW"
  cat(sprintf("=== longy %s Result Summary ===\n\n", est_label))
  print(object)

  if ("n_effective" %in% names(object$estimates)) {
    cat(sprintf("\n  Min ESS:             %.1f\n",
                min(object$estimates$n_effective, na.rm = TRUE)))
  }

  invisible(object)
}

#' Plot longy Results
#'
#' Creates a plot of IPW estimates over time with confidence intervals.
#'
#' @param x A `longy_result` object.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL (base plot).
#' @export
plot.longy_result <- function(x, ...) {
  est <- as.data.frame(x$estimates)

  # Compute y-axis range including CIs
  has_ci <- "ci_lower" %in% names(est) && "ci_upper" %in% names(est)
  if (has_ci) {
    y_range <- range(c(est$ci_lower, est$ci_upper), na.rm = TRUE)
  } else {
    y_range <- range(est$estimate, na.rm = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(est, ggplot2::aes(x = time, y = estimate)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2)

    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2
      )
    }

    p <- p +
      ggplot2::coord_cartesian(ylim = y_range) +
      ggplot2::labs(
        x = "Time", y = "Estimate",
        title = sprintf("%s Estimate -- %s",
                        if (!is.null(x$estimator) && x$estimator == "gcomp") "G-comp" else
                        if (!is.null(x$estimator) && x$estimator == "tmle") "TMLE" else "IPW",
                        x$regime)
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(p)
  }

  # Base R fallback
  plot(est$time, est$estimate, type = "b", pch = 19,
       xlab = "Time", ylab = "Estimate", ylim = y_range,
       main = sprintf("%s Estimate -- %s",
                      if (!is.null(x$estimator) && x$estimator == "gcomp") "G-comp" else
                        if (!is.null(x$estimator) && x$estimator == "tmle") "TMLE" else "IPW",
                      x$regime))
  if (has_ci) {
    graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                     angle = 90, code = 3, length = 0.05, col = "gray50")
  }
  invisible(NULL)
}
