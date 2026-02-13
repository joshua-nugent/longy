#' High-Level Wrapper for Longitudinal Causal Inference
#'
#' Runs the entire longy pipeline in one call: data setup, regime definition,
#' nuisance model fitting, and estimation via IPW, G-computation, or both.
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
#' @param outcome_type Character. \code{"binary"}, \code{"continuous"}, or \code{"survival"}.
#' @param regimes Named list. Each element defines a regime:
#'   \itemize{
#'     \item For static: an integer (0 or 1), e.g., \code{list(always = 1L, never = 0L)}
#'     \item For dynamic: a function returning 0/1
#'     \item For stochastic: a function returning P(A=1)
#'   }
#' @param estimator Character. Which estimator to use: \code{"ipw"} (default),
#'   \code{"gcomp"} (G-computation via sequential regression), \code{"tmle"}
#'   (targeted maximum likelihood), or \code{"both"} (IPW + G-comp).
#'   When \code{"both"}, results are returned with \code{_ipw} and \code{_gcomp}
#'   suffixes per regime.
#' @param covariates Character vector. Predictor columns for nuisance models.
#'   If NULL, uses all baseline + timevarying.
#' @param learners Character vector. SuperLearner library names (default
#'   \code{c("SL.glm", "SL.mean")}). Set to NULL to use plain glm without
#'   SuperLearner.
#' @param stabilized Logical. Use stabilized weights (IPW only).
#' @param truncation Numeric. Weight truncation cap (IPW only).
#' @param truncation_quantile Numeric. Quantile-based weight truncation (IPW only).
#' @param inference Character. \code{"ic"}, \code{"bootstrap"}, or \code{"sandwich"}.
#'   Ignored for G-comp (always bootstrap).
#' @param ci_level Numeric. Confidence level.
#' @param n_boot Integer. Bootstrap replicates.
#' @param cluster Character. Cluster column for robust SEs (IPW only).
#' @param times Numeric vector. Time points for estimation. NULL = all.
#' @param adaptive_cv Logical. Adaptive CV fold selection.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param bounds Numeric(2). Prediction probability bounds.
#' @param g_bounds Numeric(2). Bounds for cumulative g (denominator of clever
#'   covariate). Default \code{c(0.01, 1)}. Used by TMLE and IPW.
#' @param outcome_range Numeric(2) or NULL. Range for scaling continuous
#'   outcomes to \code{[0,1]} in TMLE. If NULL, uses empirical range.
#'   Ignored for binary/survival outcomes.
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class \code{"longy_results"} (a named list of
#'   \code{longy_result} objects, one per regime). When \code{estimator = "both"},
#'   each regime produces two entries with \code{_ipw} and \code{_gcomp} suffixes.
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
                  estimator = "ipw",
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
                  g_bounds = c(0.01, 1),
                  outcome_range = NULL,
                  verbose = TRUE) {

  estimator <- match.arg(estimator, c("ipw", "gcomp", "tmle", "both"))
  do_ipw <- estimator %in% c("ipw", "both")
  do_gcomp <- estimator %in% c("gcomp", "both")
  do_tmle <- estimator == "tmle"

  # n_boot = 0 means point estimates only (no inference)
  # For TMLE, default to EIF inference (cheap) even with n_boot = 0
  if (n_boot == 0L && !do_tmle) inference <- "none"
  if (do_tmle && inference == "ic") inference <- "eif"  # map IC -> EIF for TMLE

  # Determine total steps for verbose messaging
  n_steps <- 2L  # data + regimes always
  if (do_ipw) n_steps <- n_steps + 4L  # g_A + g_C + g_R + weights/estimate
  if (do_gcomp) n_steps <- n_steps + 1L  # outcome model + estimate
  if (do_tmle) n_steps <- n_steps + 4L  # g_A + g_C + outcome + TMLE estimate

  step <- 0L

  # Step 1: Create longy_data object
  step <- step + 1L
  if (verbose) .vmsg("Step %d/%d: Creating longy_data object...", step, n_steps)
  obj <- longy_data(
    data = data, id = id, time = time,
    outcome = outcome, treatment = treatment,
    censoring = censoring, observation = observation,
    baseline = baseline, timevarying = timevarying,
    sampling_weights = sampling_weights,
    outcome_type = outcome_type, verbose = verbose
  )

  # Step 2: Define regimes
  step <- step + 1L
  if (verbose) .vmsg("Step %d/%d: Defining regimes...", step, n_steps)
  for (rname in names(regimes)) {
    rval <- regimes[[rname]]
    if (is.numeric(rval) && length(rval) == 1) {
      obj <- define_regime(obj, name = rname, static = as.integer(rval))
    } else if (is.function(rval)) {
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

    cur_step <- step  # save so IPW and gcomp branches share the base

    # --- IPW pipeline ---
    if (do_ipw) {
      if (verbose) .vmsg("Step %d/%d: Fitting treatment model (g_A)...",
                          cur_step + 1L, n_steps)
      r_obj <- fit_treatment(r_obj, regime = rname, covariates = covariates,
                             learners = learners, adaptive_cv = adaptive_cv,
                             min_obs = min_obs, bounds = bounds,
                             times = times, verbose = verbose)

      if (verbose) .vmsg("Step %d/%d: Fitting censoring model (g_C)...",
                          cur_step + 2L, n_steps)
      r_obj <- fit_censoring(r_obj, regime = rname, covariates = covariates,
                             learners = learners, adaptive_cv = adaptive_cv,
                             min_obs = min_obs, bounds = bounds,
                             times = times, verbose = verbose)

      if (verbose) .vmsg("Step %d/%d: Fitting observation model (g_R)...",
                          cur_step + 3L, n_steps)
      r_obj <- fit_observation(r_obj, regime = rname, covariates = covariates,
                               learners = learners, adaptive_cv = adaptive_cv,
                               min_obs = min_obs, bounds = bounds,
                               times = times, verbose = verbose)

      if (verbose) .vmsg("Step %d/%d: Computing weights and estimating (IPW)...",
                          cur_step + 4L, n_steps)
      r_obj <- compute_weights(r_obj, regime = rname,
                               stabilized = stabilized,
                               truncation = truncation,
                               truncation_quantile = truncation_quantile,
                               g_bounds = g_bounds)

      ipw_result <- estimate_ipw(r_obj, regime = rname, times = times,
                                 inference = inference, ci_level = ci_level,
                                 n_boot = n_boot, cluster = cluster)

      result_name <- if (estimator == "both") paste0(rname, "_ipw") else rname
      all_results[[result_name]] <- ipw_result
    }

    # --- G-comp pipeline ---
    if (do_gcomp) {
      # Use a copy that doesn't require IPW fits
      gc_obj <- obj
      gc_obj$fits <- list(treatment = NULL, censoring = list(), observation = NULL)
      gc_obj$weights <- NULL
      gc_obj$regimes <- obj$regimes

      gc_step <- cur_step + (if (do_ipw) 4L else 0L) + 1L
      if (verbose) .vmsg("Step %d/%d: Fitting outcome model and estimating (G-comp)...",
                          gc_step, n_steps)

      gc_obj <- fit_outcome(gc_obj, regime = rname, covariates = covariates,
                            learners = learners, adaptive_cv = adaptive_cv,
                            min_obs = min_obs, bounds = bounds,
                            times = times, verbose = verbose)

      gcomp_result <- estimate_gcomp(gc_obj, regime = rname, times = times,
                                     ci_level = ci_level, n_boot = n_boot,
                                     verbose = verbose)

      result_name <- if (estimator == "both") paste0(rname, "_gcomp") else rname
      all_results[[result_name]] <- gcomp_result
    }

    # --- TMLE pipeline ---
    if (do_tmle) {
      tm_obj <- obj
      tm_obj$fits <- list(treatment = NULL, censoring = list(), observation = NULL)
      tm_obj$weights <- NULL
      tm_obj$regimes <- obj$regimes

      tm_step <- cur_step + 1L
      if (verbose) .vmsg("Step %d/%d: Fitting treatment model (g_A) for TMLE...",
                          tm_step, n_steps)
      tm_obj <- fit_treatment(tm_obj, regime = rname, covariates = covariates,
                              learners = learners, adaptive_cv = adaptive_cv,
                              min_obs = min_obs, bounds = bounds,
                              times = times, verbose = verbose)

      if (verbose) .vmsg("Step %d/%d: Fitting censoring model (g_C) for TMLE...",
                          tm_step + 1L, n_steps)
      tm_obj <- fit_censoring(tm_obj, regime = rname, covariates = covariates,
                              learners = learners, adaptive_cv = adaptive_cv,
                              min_obs = min_obs, bounds = bounds,
                              times = times, verbose = verbose)

      if (verbose) .vmsg("Step %d/%d: Fitting outcome model for TMLE...",
                          tm_step + 2L, n_steps)
      tm_obj <- fit_outcome(tm_obj, regime = rname, covariates = covariates,
                            learners = learners, adaptive_cv = adaptive_cv,
                            min_obs = min_obs, bounds = bounds,
                            times = times, verbose = verbose)

      # Determine TMLE inference
      tmle_inf <- if (n_boot == 0L) "eif" else inference
      if (tmle_inf == "none") tmle_inf <- "eif"  # EIF is always cheap for TMLE

      if (verbose) .vmsg("Step %d/%d: Running TMLE estimation...",
                          tm_step + 3L, n_steps)
      tmle_result <- estimate_tmle(tm_obj, regime = rname, times = times,
                                   inference = tmle_inf, ci_level = ci_level,
                                   n_boot = n_boot, g_bounds = g_bounds,
                                   outcome_range = outcome_range,
                                   verbose = verbose)

      all_results[[rname]] <- tmle_result
    }
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
#' Creates a plot comparing estimates over time across all regimes,
#' with confidence intervals shown as ribbons.
#'
#' @param x A \code{longy_results} object (list of \code{longy_result} objects).
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
        title = "Estimates by Regime"
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
           main = "Estimates by Regime")
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
