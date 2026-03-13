# S3 methods for longy_result objects
# These are generic to all estimators (IPW, G-comp, TMLE, unadjusted)

#' Get display label for an estimator
#' @param x A longy_result object or a character estimator name.
#' @return Character label for display.
#' @noRd
.estimator_label <- function(x) {
  est <- if (is.list(x)) x$estimator else x
  switch(est,
    gcomp = "G-comp",
    tmle = "TMLE",
    unadjusted = "Unadjusted",
    "IPW"
  )
}

#' @export
print.longy_result <- function(x, ...) {
  est_label <- .estimator_label(x)
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
    n_eff = if ("n_effective" %in% names(est)) round(est$n_effective, 1) else est$n_at_risk,
    n_risk = est$n_at_risk
  )
  print(disp, row.names = FALSE)
  invisible(x)
}

#' @export
summary.longy_result <- function(object, ...) {
  est_label <- .estimator_label(object)
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
#' Creates a plot of estimates over time with confidence intervals.
#'
#' @param x A `longy_result` object.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL (base plot).
#' @export
plot.longy_result <- function(x, ...) {
  est <- as.data.frame(x$estimates)
  est_label <- .estimator_label(x)

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
        title = sprintf("%s Estimate -- %s", est_label, x$regime)
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(p)
  }

  # Base R fallback
  plot(est$time, est$estimate, type = "b", pch = 19,
       xlab = "Time", ylab = "Estimate", ylim = y_range,
       main = sprintf("%s Estimate -- %s", est_label, x$regime))
  if (has_ci) {
    graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                     angle = 90, code = 3, length = 0.05, col = "gray50")
  }
  invisible(NULL)
}
