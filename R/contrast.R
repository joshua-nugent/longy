#' Compute Treatment Effect Contrasts Between Regimes
#'
#' Compares two regimes under the same estimator to produce a treatment effect
#' estimate (e.g., risk difference, risk ratio, odds ratio). Inference uses the
#' delta method applied to stored per-subject influence curves when available.
#'
#' @param obj A \code{longy_data} object with results for at least two regimes.
#' @param regime Character(2). Names of the two regimes to contrast. The first
#'   is the "treatment" regime and the second is the "reference" regime.
#' @param ref Character(1). Reference regime name. If provided, overrides the
#'   second element of \code{regime}. Follows lmtp's convention:
#'   \code{contrast = regime - ref} on the difference scale.
#' @param estimator Character. Which estimator's results to contrast
#'   (\code{"ipw"}, \code{"tmle"}, \code{"gcomp"}, \code{"unadjusted"}).
#'   If NULL, auto-detects from available results.
#' @param scale Character. Scale of the contrast: \code{"difference"} (default,
#'   ATE/risk difference), \code{"ratio"} (risk ratio), or \code{"odds_ratio"}.
#' @param ci_level Numeric. Confidence level for delta-method CIs. Default 0.95.
#'
#' @return An S3 object of class \code{"longy_contrast"} with elements:
#'   \describe{
#'     \item{estimates}{data.table with time, estimate, se, ci_lower, ci_upper}
#'     \item{regime}{Character(2): treatment and reference regime names}
#'     \item{estimator}{Estimator used}
#'     \item{scale}{Contrast scale}
#'     \item{ci_level}{Confidence level}
#'     \item{inference}{Method used: \code{"delta_method"} or \code{"none"}}
#'   }
#'
#' @details
#' \strong{Inference:} When both regime results have stored influence curves
#' (from IC-based or EIF inference), the delta method is used:
#' \itemize{
#'   \item \strong{Difference:} IC_contrast = IC_1 - IC_0; SE = sqrt(var(IC) / n)
#'   \item \strong{Ratio:} IC_ratio = (IC_1 * psi_0 - IC_0 * psi_1) / psi_0^2
#'   \item \strong{Odds ratio:} Delta method on log-odds scale
#' }
#' When ICs are not available (e.g., bootstrap or G-comp), only point estimates
#' are returned with a message.
#'
#' @examples
#' \dontrun{
#' obj <- longy(data, id = "id", time = "time", outcome = "Y",
#'              treatment = "A", censoring = "C",
#'              baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'              regimes = list(always = 1L, never = 0L),
#'              estimator = "tmle")
#'
#' # Risk difference (ATE)
#' contrast(obj, regime = c("always", "never"))
#'
#' # Equivalently, using ref=
#' contrast(obj, regime = "always", ref = "never")
#'
#' # Risk ratio
#' contrast(obj, regime = c("always", "never"), scale = "ratio")
#' }
#'
#' @export
contrast <- function(obj, regime, ref = NULL, estimator = NULL,
                     scale = c("difference", "ratio", "odds_ratio"),
                     ci_level = 0.95) {
  obj <- .as_longy_data(obj)
  scale <- match.arg(scale)

  if (ci_level <= 0 || ci_level >= 1)
    stop("ci_level must be between 0 and 1.", call. = FALSE)

  # Resolve regime pair
  if (!is.null(ref)) {
    if (length(regime) != 1)
      stop("When 'ref' is provided, 'regime' must be a single regime name.",
           call. = FALSE)
    regime <- c(regime, ref)
  }
  if (length(regime) != 2)
    stop("'regime' must be a character vector of length 2 (treatment, reference).",
         call. = FALSE)
  if (regime[1] == regime[2])
    stop("Cannot contrast a regime with itself.", call. = FALSE)

  # Auto-detect estimator
  if (is.null(estimator)) {
    # Find estimators with results for BOTH regimes
    r1_keys <- grep(paste0("^", regime[1], "_"), names(obj$results), value = TRUE)
    r2_keys <- grep(paste0("^", regime[2], "_"), names(obj$results), value = TRUE)
    r1_ests <- sub(paste0("^", regime[1], "_"), "", r1_keys)
    r2_ests <- sub(paste0("^", regime[2], "_"), "", r2_keys)
    common <- intersect(r1_ests, r2_ests)
    common <- setdiff(common, "unadjusted")
    if (length(common) == 0) {
      # Fall back to unadjusted
      common <- intersect(r1_ests, r2_ests)
    }
    if (length(common) == 0)
      stop(sprintf("No common estimator results found for regimes '%s' and '%s'.",
                   regime[1], regime[2]), call. = FALSE)
    # Prefer tmle > ipw > gcomp > unadjusted
    pref_order <- c("tmle", "ipw", "gcomp", "unadjusted")
    estimator <- pref_order[pref_order %in% common][1]
  }
  estimator <- match.arg(estimator, c("ipw", "tmle", "gcomp", "unadjusted"))

  # Extract results
  key1 <- paste0(regime[1], "_", estimator)
  key2 <- paste0(regime[2], "_", estimator)
  res1 <- obj$results[[key1]]
  res2 <- obj$results[[key2]]
  if (is.null(res1))
    stop(sprintf("No %s result for regime '%s'.", toupper(estimator), regime[1]),
         call. = FALSE)
  if (is.null(res2))
    stop(sprintf("No %s result for regime '%s'.", toupper(estimator), regime[2]),
         call. = FALSE)

  # Merge estimates by time
  est1 <- res1$estimates
  est2 <- res2$estimates
  common_times <- intersect(est1$time, est2$time)
  if (length(common_times) == 0)
    stop("No overlapping time points between the two regime results.",
         call. = FALSE)

  est1 <- est1[est1$time %in% common_times, ]
  est2 <- est2[est2$time %in% common_times, ]
  data.table::setorder(est1, time)
  data.table::setorder(est2, time)

  # Compute point estimates
  psi1 <- est1$estimate
  psi2 <- est2$estimate

  contrast_est <- switch(scale,
    difference = psi1 - psi2,
    ratio = psi1 / psi2,
    odds_ratio = (psi1 / (1 - psi1)) / (psi2 / (1 - psi2))
  )

  # Delta-method inference using stored ICs
  ic1 <- res1$ic
  ic2 <- res2$ic
  has_ic <- !is.null(ic1) && !is.null(ic2)
  inference_method <- "none"

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  id_col <- obj$nodes$id
  cluster_col <- obj$nodes$cluster  # NULL if no clustering

  if (has_ic) {
    inference_method <- "delta_method"
    se_vec <- numeric(length(common_times))
    ci_lower <- numeric(length(common_times))
    ci_upper <- numeric(length(common_times))

    for (k in seq_along(common_times)) {
      tt <- common_times[k]
      ic1_t <- ic1[ic1$.time == tt, ]
      ic2_t <- ic2[ic2$.time == tt, ]

      # Align by subject ID
      merged_ic <- merge(ic1_t, ic2_t, by = c(id_col, ".time"),
                         suffixes = c("_1", "_2"))

      if (nrow(merged_ic) < 2) {
        se_vec[k] <- NA_real_
        ci_lower[k] <- NA_real_
        ci_upper[k] <- NA_real_
        next
      }

      d1 <- merged_ic$.ic_1
      d2 <- merged_ic$.ic_2

      # Delta method IC for the contrast
      ic_contrast <- switch(scale,
        difference = d1 - d2,
        ratio = (d1 * psi2[k] - d2 * psi1[k]) / psi2[k]^2,
        odds_ratio = {
          # Log-odds ratio: delta method on log scale, then exponentiate
          # d/dpsi log(psi/(1-psi)) = 1/(psi*(1-psi))
          dlogodds1 <- d1 / (psi1[k] * (1 - psi1[k]))
          dlogodds2 <- d2 / (psi2[k] * (1 - psi2[k]))
          dlogodds1 - dlogodds2
        }
      )

      # Cluster-robust SE when cluster is specified
      se_k <- .ic_se(ic_contrast, ids = merged_ic[[id_col]],
                     cluster_col = cluster_col, obj = obj)

      se_vec[k] <- se_k

      if (scale == "odds_ratio" && !is.na(se_k)) {
        # CIs on log-odds-ratio scale, then exponentiate
        log_or <- log(contrast_est[k])
        ci_lower[k] <- exp(log_or - z * se_k)
        ci_upper[k] <- exp(log_or + z * se_k)
        # SE on OR scale (for display)
        se_vec[k] <- contrast_est[k] * se_k
      } else if (!is.na(se_k)) {
        ci_lower[k] <- contrast_est[k] - z * se_k
        ci_upper[k] <- contrast_est[k] + z * se_k
      } else {
        ci_lower[k] <- NA_real_
        ci_upper[k] <- NA_real_
      }
    }
  } else {
    if (estimator %in% c("ipw", "tmle")) {
      message(sprintf(
        "No influence curves stored for %s results. Re-run with inference='ic' or 'eif' to enable delta-method contrasts.",
        toupper(estimator)))
    } else {
      message(sprintf(
        "%s does not produce influence curves. Contrast SEs require bootstrap (not yet implemented for contrasts).",
        .estimator_label(estimator)))
    }
    se_vec <- rep(NA_real_, length(common_times))
    ci_lower <- rep(NA_real_, length(common_times))
    ci_upper <- rep(NA_real_, length(common_times))
  }

  out_dt <- data.table::data.table(
    time = common_times,
    estimate = contrast_est,
    se = se_vec,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

  result <- list(
    estimates = out_dt,
    regime = regime,
    estimator = estimator,
    scale = scale,
    ci_level = ci_level,
    inference = inference_method
  )
  class(result) <- "longy_contrast"
  result
}

# --- S3 methods for longy_contrast ---

#' @export
print.longy_contrast <- function(x, ...) {
  est_label <- .estimator_label(x$estimator)
  scale_label <- switch(x$scale,
    difference = "Risk Difference",
    ratio = "Risk Ratio",
    odds_ratio = "Odds Ratio"
  )
  cat(sprintf("longy %s contrast (%s)\n", est_label, scale_label))
  cat(sprintf("  Treatment: %s | Reference: %s\n", x$regime[1], x$regime[2]))
  cat(sprintf("  Inference: %s | CI level: %.0f%%\n\n",
              x$inference, x$ci_level * 100))

  est <- x$estimates
  disp <- data.frame(
    time = est$time,
    estimate = round(est$estimate, 4),
    se = if (all(is.na(est$se))) NULL else round(est$se, 4),
    ci_lower = if (all(is.na(est$ci_lower))) NULL else round(est$ci_lower, 4),
    ci_upper = if (all(is.na(est$ci_upper))) NULL else round(est$ci_upper, 4)
  )
  print(disp, row.names = FALSE)
  invisible(x)
}

#' @export
summary.longy_contrast <- function(object, ...) {
  cat("=== longy Contrast Summary ===\n\n")
  print(object)
  invisible(object)
}

#' @export
plot.longy_contrast <- function(x, ...) {
  est <- as.data.frame(x$estimates)
  scale_label <- switch(x$scale,
    difference = "Risk Difference",
    ratio = "Risk Ratio",
    odds_ratio = "Odds Ratio"
  )
  est_label <- .estimator_label(x$estimator)
  has_ci <- !all(is.na(est$ci_lower)) && !all(is.na(est$ci_upper))

  if (has_ci) {
    y_range <- range(c(est$ci_lower, est$ci_upper), na.rm = TRUE)
  } else {
    y_range <- range(est$estimate, na.rm = TRUE)
  }

  # Add null reference line
  null_val <- switch(x$scale,
    difference = 0,
    ratio = 1,
    odds_ratio = 1
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(est, ggplot2::aes(x = time, y = estimate)) +
      ggplot2::geom_hline(yintercept = null_val, linetype = "dashed",
                          colour = "gray50") +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2)

    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2
      )
    }

    p <- p +
      ggplot2::labs(
        x = "Time", y = scale_label,
        title = sprintf("%s %s: %s vs %s",
                        est_label, scale_label, x$regime[1], x$regime[2])
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(p)
  }

  # Base R fallback
  plot(est$time, est$estimate, type = "b", pch = 19,
       xlab = "Time", ylab = scale_label, ylim = y_range,
       main = sprintf("%s %s: %s vs %s",
                      est_label, scale_label, x$regime[1], x$regime[2]))
  graphics::abline(h = null_val, lty = 2, col = "gray50")
  if (has_ci) {
    graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                     angle = 90, code = 3, length = 0.05, col = "gray50")
  }
  invisible(NULL)
}
