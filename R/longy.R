#' High-Level Wrapper for Longitudinal Causal Inference
#'
#' Runs the entire longy pipeline in one call: data setup, regime definition,
#' nuisance model fitting, weight computation, and IPW estimation.
#'
#' @param data A data.frame or data.table in long format.
#' @param id Character. Subject identifier column.
#' @param time Character. Time index column.
#' @param outcome Character. Outcome column.
#' @param treatment Character. Binary treatment column.
#' @param censoring Character vector. Censoring column(s) (absorbing). NULL if none.
#' @param observation Character. Intermittent observation column. NULL if outcome
#'   always observed.
#' @param baseline Character vector. Baseline covariate columns.
#' @param timevarying Character vector. Time-varying covariate columns.
#' @param sampling_weights Character. Column name for external sampling/survey
#'   weights. See [longy_data()] for details. NULL if none.
#' @param outcome_type Character. `"binary"`, `"continuous"`, or `"survival"`.
#' @param regimes Named list. Each element defines a regime:
#'   \itemize{
#'     \item For static: an integer (0 or 1), e.g., `list(always = 1L, never = 0L)`
#'     \item For dynamic: a function returning 0/1
#'     \item For stochastic: a function returning P(A=1)
#'   }
#' @param covariates Character vector. Predictor columns for nuisance models.
#'   If NULL, uses all baseline + timevarying.
#' @param learners Character vector. SuperLearner library names (default
#'   `c("SL.glm", "SL.mean")`). Set to NULL to use plain glm without
#'   SuperLearner.
#' @param stabilized Logical. Use stabilized weights.
#' @param truncation Numeric. Weight truncation cap.
#' @param truncation_quantile Numeric. Quantile-based weight truncation.
#' @param inference Character. `"ic"`, `"bootstrap"`, or `"sandwich"`.
#' @param ci_level Numeric. Confidence level.
#' @param n_boot Integer. Bootstrap replicates.
#' @param cluster Character. Cluster column for robust SEs.
#' @param times Numeric vector. Time points for estimation. NULL = all.
#' @param adaptive_cv Logical. Adaptive CV fold selection.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param bounds Numeric(2). Prediction probability bounds.
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class `"longy_results"` (a named list of
#'   `longy_result` objects, one per regime).
#'
#' @examples
#' \dontrun{
#' results <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C", observation = "R",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, never = 0L)
#' )
#' }
#'
#' @export
longy <- function(data,
                  id, time, outcome, treatment,
                  censoring = NULL, observation = NULL,
                  baseline = character(0), timevarying = character(0),
                  sampling_weights = NULL,
                  outcome_type = "binary",
                  regimes = list(always = 1L, never = 0L),
                  covariates = NULL,
                  learners = c("SL.glm", "SL.mean"),
                  stabilized = TRUE,
                  truncation = NULL,
                  truncation_quantile = NULL,
                  inference = "ic",
                  ci_level = 0.95,
                  n_boot = 200L,
                  cluster = NULL,
                  times = NULL,
                  adaptive_cv = TRUE,
                  min_obs = 50L,
                  bounds = c(0.005, 0.995),
                  verbose = TRUE) {

  # Step 1: Create longy_data object
  if (verbose) .vmsg("Step 1/6: Creating longy_data object...")
  obj <- longy_data(
    data = data, id = id, time = time,
    outcome = outcome, treatment = treatment,
    censoring = censoring, observation = observation,
    baseline = baseline, timevarying = timevarying,
    sampling_weights = sampling_weights,
    outcome_type = outcome_type, verbose = verbose
  )

  # Step 2: Define regimes
  if (verbose) .vmsg("Step 2/6: Defining regimes...")
  for (rname in names(regimes)) {
    rval <- regimes[[rname]]
    if (is.numeric(rval) && length(rval) == 1) {
      obj <- define_regime(obj, name = rname, static = as.integer(rval))
    } else if (is.function(rval)) {
      # Determine if dynamic or stochastic by testing output
      obj <- define_regime(obj, name = rname, dynamic = rval)
    } else {
      stop(sprintf("Regime '%s': must be an integer (0/1) or a function.", rname),
           call. = FALSE)
    }
  }

  # Fit and estimate for each regime
  all_results <- list()

  for (rname in names(regimes)) {
    if (verbose) .vmsg("\n=== Regime: %s ===", rname)

    # We need a fresh copy for each regime since fits are regime-specific
    r_obj <- obj
    r_obj$fits <- list(treatment = NULL, censoring = list(), observation = NULL)
    r_obj$weights <- NULL

    # Step 3: Fit treatment model
    if (verbose) .vmsg("Step 3/6: Fitting treatment model (g_A)...")
    r_obj <- fit_treatment(r_obj, regime = rname, covariates = covariates,
                           learners = learners, adaptive_cv = adaptive_cv,
                           min_obs = min_obs, bounds = bounds,
                           times = times, verbose = verbose)

    # Step 4: Fit censoring model
    if (verbose) .vmsg("Step 4/6: Fitting censoring model (g_C)...")
    r_obj <- fit_censoring(r_obj, regime = rname, covariates = covariates,
                           learners = learners, adaptive_cv = adaptive_cv,
                           min_obs = min_obs, bounds = bounds,
                           times = times, verbose = verbose)

    # Step 5: Fit observation model
    if (verbose) .vmsg("Step 5/6: Fitting observation model (g_R)...")
    r_obj <- fit_observation(r_obj, regime = rname, covariates = covariates,
                             learners = learners, adaptive_cv = adaptive_cv,
                             min_obs = min_obs, bounds = bounds,
                             times = times, verbose = verbose)

    # Step 6: Compute weights
    if (verbose) .vmsg("Step 6/6: Computing weights and estimating...")
    r_obj <- compute_weights(r_obj, regime = rname,
                             stabilized = stabilized,
                             truncation = truncation,
                             truncation_quantile = truncation_quantile)

    # Estimate
    result <- estimate_ipw(r_obj, regime = rname, times = times,
                           inference = inference, ci_level = ci_level,
                           n_boot = n_boot, cluster = cluster)

    all_results[[rname]] <- result
  }

  class(all_results) <- "longy_results"
  all_results
}

#' @export
print.longy_results <- function(x, ...) {
  cat(sprintf("longy results: %d regime(s)\n\n", length(x)))
  for (rname in names(x)) {
    cat(sprintf("--- %s ---\n", rname))
    print(x[[rname]])
    cat("\n")
  }
  invisible(x)
}

#' Plot longy Results Across Regimes
#'
#' Creates a plot comparing IPW estimates over time across all regimes,
#' with confidence intervals shown as ribbons.
#'
#' @param x A `longy_results` object (list of `longy_result` objects).
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL (base plot).
#' @export
plot.longy_results <- function(x, ...) {
  # Combine estimates from all regimes into one data.frame
  all_est <- lapply(names(x), function(rname) {
    est <- as.data.frame(x[[rname]]$estimates)
    est$regime <- rname
    est
  })
  combined <- do.call(rbind, all_est)

  # Compute y-axis range including CIs
  has_ci <- "ci_lower" %in% names(combined) && "ci_upper" %in% names(combined)
  if (has_ci) {
    y_range <- range(c(combined$ci_lower, combined$ci_upper), na.rm = TRUE)
  } else {
    y_range <- range(combined$estimate, na.rm = TRUE)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(combined,
           ggplot2::aes(x = time, y = estimate, colour = regime, fill = regime))

    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.15, colour = NA
      )
    }

    p <- p +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::coord_cartesian(ylim = y_range) +
      ggplot2::labs(
        x = "Time", y = "Estimate",
        colour = "Regime", fill = "Regime",
        title = "IPW Estimates by Regime"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(p)
  }

  # Base R fallback
  regimes <- names(x)
  cols <- seq_along(regimes)

  first <- TRUE
  for (i in seq_along(regimes)) {
    est <- combined[combined$regime == regimes[i], ]
    if (first) {
      plot(est$time, est$estimate, type = "b", pch = 19, col = cols[i],
           xlab = "Time", ylab = "Estimate", ylim = y_range,
           main = "IPW Estimates by Regime")
      first <- FALSE
    } else {
      graphics::lines(est$time, est$estimate, type = "b", pch = 19, col = cols[i])
    }
    if ("ci_lower" %in% names(est) && "ci_upper" %in% names(est)) {
      graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                       angle = 90, code = 3, length = 0.05, col = cols[i])
    }
  }
  graphics::legend("topleft", legend = regimes, col = cols, lty = 1, pch = 19)
  invisible(NULL)
}
