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
#' @param censoring Character (length 1). Column name for a character/factor
#'   censoring status column. Must contain \code{"uncensored"} for non-censored
#'   rows; other values (e.g. \code{"censored"}, \code{"death"}, \code{"ltfu"})
#'   are treated as distinct censoring causes, each modeled separately. See
#'   \code{\link{longy_data}} for full details. NULL if no censoring.
#' @param observation Character. Intermittent observation column. NULL to
#'   auto-detect from NA values in the outcome column (see
#'   \code{\link{longy_data}}).
#' @param baseline Character vector. Baseline covariate columns.
#' @param timevarying Character vector. Time-varying covariate columns.
#' @param sampling_weights Character. Column name for external sampling/survey
#'   weights. See [longy_data()] for details. NULL if none.
#' @param outcome_type Character. \code{"binary"}, \code{"continuous"}, or \code{"survival"}.
#' @param competing Character. Column name for a binary absorbing competing
#'   event indicator. Used with \code{outcome_type = "survival"} to estimate
#'   cause-specific cumulative incidence. See [longy_data()] for details. NULL
#'   if no competing risks.
#' @param regimes Named list. Each element defines a regime:
#'   \itemize{
#'     \item For static: an integer (0 or 1), e.g., \code{list(always = 1L, never = 0L)}
#'     \item For dynamic: a function returning 0/1
#'     \item For stochastic: a function returning P(A=1)
#'   }
#' @param estimator Character. Which estimator(s) to use: \code{"ipw"} (default),
#'   \code{"gcomp"}, \code{"tmle"}, \code{"both"} (IPW + G-comp, backward
#'   compatible), or \code{"all"} (IPW + G-comp + TMLE). When multiple
#'   estimators are run, results are returned with \code{_ipw}, \code{_gcomp},
#'   and/or \code{_tmle} suffixes per regime. Nuisance models are shared
#'   across estimators to avoid redundant fitting.
#' @param covariates Character vector. Predictor columns for nuisance models.
#'   If NULL, uses all baseline + timevarying.
#' @param learners Character vector or named list. SuperLearner library names
#'   (default \code{c("SL.glm", "SL.mean")}). Set to NULL to use plain glm
#'   without SuperLearner. When a named list, valid keys are \code{treatment},
#'   \code{censoring}, \code{observation}, \code{outcome}, and \code{default}.
#'   Missing keys fall back to \code{default}, which itself falls back to
#'   \code{c("SL.glm", "SL.mean")}. Example:
#'   \code{list(default = c("SL.glm", "SL.gam"), outcome = c("SL.glm", "SL.gam", "SL.earth"))}.
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
#' @param min_events Integer. Minimum minority-class events required to fit a
#'   model. When the minority class count is below this AND the minority rate
#'   is below 0.01, marginal fallback is used. Default 20.
#' @param bounds Numeric(2). Prediction probability bounds.
#' @param g_bounds Numeric(2). Bounds for cumulative g (denominator of clever
#'   covariate). Default \code{c(0.01, 1)}. Used by TMLE and IPW.
#' @param outcome_range Numeric(2) or NULL. Range for scaling continuous
#'   outcomes to \code{[0,1]} in TMLE. If NULL, uses empirical range.
#'   Ignored for binary/survival outcomes.
#' @param cross_fit Integer or NULL. Number of cross-fitting folds. When
#'   non-NULL, nuisance models (g_A, g_C, g_R) are fit on training folds and
#'   predicted on validation folds, and TMLE uses cross-fitted Q models.
#'   This eliminates overfitting bias in EIF-based inference. Typical values
#'   are 5 or 10. Default NULL (no cross-fitting).
#' @param cross_fit_seed Integer or NULL. Random seed for cross-fitting fold
#'   assignment. NULL for no seed.
#' @param sl_fn Character. Which SuperLearner implementation to use:
#'   \code{"SuperLearner"} (default, sequential CV) or \code{"ffSL"}
#'   (future-factorial, parallelizes fold x algorithm combinations via
#'   \code{future.apply}).
#' @param k Integer or \code{Inf}. Number of lagged time steps of covariate
#'   history to include as additional predictors. See \code{\link{longy_data}}
#'   for details. Default \code{Inf} (all available history).
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class \code{"longy_data"} with estimation results
#'   accumulated in \code{obj$results} (a named list of \code{longy_result}
#'   objects keyed as \code{{regime}_{estimator}}). Use \code{results(obj)} to
#'   access results, or \code{obj$results$always_ipw$estimates} directly.
#'
#' @examples
#' \dontrun{
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C", observation = "R",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, never = 0L)
#' )
#' obj$results$always_ipw$estimates
#' results(obj, estimator = "ipw")
#' }
#'
#' @export
longy <- function(data,
                  id, time, outcome, treatment,
                  censoring = NULL, observation = NULL,
                  baseline = character(0), timevarying = character(0),
                  sampling_weights = NULL,
                  outcome_type = "binary",
                  competing = NULL,
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
                  min_events = 20L,
                  bounds = c(0.005, 0.995),
                  g_bounds = c(0.01, 1),
                  outcome_range = NULL,
                  cross_fit = NULL,
                  cross_fit_seed = NULL,
                  sl_fn = "ffSL",
                  k = Inf,
                  verbose = TRUE) {

  sl_fn <- match.arg(sl_fn, c("SuperLearner", "ffSL"))

  # Resolve per-model learner libraries
  if (is.list(learners) && !is.null(names(learners))) {
    default_lib <- if (!is.null(learners$default)) learners$default else c("SL.glm", "SL.mean")
    learners_treatment   <- if (!is.null(learners$treatment))   learners$treatment   else default_lib
    learners_censoring   <- if (!is.null(learners$censoring))   learners$censoring   else default_lib
    learners_observation <- if (!is.null(learners$observation))  learners$observation  else default_lib
    learners_outcome     <- if (!is.null(learners$outcome))      learners$outcome      else default_lib
  } else {
    learners_treatment <- learners_censoring <-
      learners_observation <- learners_outcome <- learners
  }

  estimator <- match.arg(estimator, c("ipw", "gcomp", "tmle", "both", "all"))
  do_ipw <- estimator %in% c("ipw", "both", "all")
  do_gcomp <- estimator %in% c("gcomp", "both", "all")
  do_tmle <- estimator %in% c("tmle", "all")
  multi <- sum(do_ipw, do_gcomp, do_tmle) > 1

  # n_boot = 0 only disables bootstrap inference; analytic methods (ic, sandwich,

  # eif) are unaffected since they don't need resampling
  if (n_boot == 0L && inference == "bootstrap") inference <- "none"

  # Determine total steps for verbose messaging
  # Shared fits: g_A + g_C fitted once if IPW or TMLE need them
  # g_R shared for IPW and TMLE; outcome fitted once if G-comp or TMLE need it
  n_steps <- 2L  # data + regimes always
  need_g <- do_ipw || do_tmle
  need_outcome <- do_gcomp || do_tmle
  if (need_g) n_steps <- n_steps + 3L  # g_A + g_C + g_R (shared)
  if (do_ipw) n_steps <- n_steps + 1L  # weights/estimate
  if (need_outcome) n_steps <- n_steps + 1L  # outcome model (shared)
  if (do_gcomp) n_steps <- n_steps + 1L  # G-comp estimate
  if (do_tmle) n_steps <- n_steps + 1L  # TMLE estimate

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
    outcome_type = outcome_type, competing = competing,
    k = k, verbose = verbose
  )

  # Set up cross-fitting if requested
  if (!is.null(cross_fit)) {
    if (verbose) .vmsg("  Setting up %d-fold cross-fitting...", cross_fit)
    obj <- set_crossfit(obj, n_folds = as.integer(cross_fit),
                        seed = cross_fit_seed)
  }

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

  # Fit all regimes at once, then estimate
  regime_names <- names(regimes)
  cur_step <- step

  # --- Shared nuisance models: g_A + g_C + g_R (IPW and/or TMLE) ---
  if (need_g) {
    if (verbose) .vmsg("Step %d/%d: Fitting treatment model (g_A)...",
                        cur_step + 1L, n_steps)
    obj <- fit_treatment(obj, regime = regime_names, covariates = covariates,
                           learners = learners_treatment, adaptive_cv = adaptive_cv,
                           min_obs = min_obs, min_events = min_events,
                           bounds = bounds,
                           times = times, sl_fn = sl_fn,
                           verbose = verbose)

    if (verbose) .vmsg("Step %d/%d: Fitting censoring model (g_C)...",
                        cur_step + 2L, n_steps)
    obj <- fit_censoring(obj, regime = regime_names, covariates = covariates,
                           learners = learners_censoring, adaptive_cv = adaptive_cv,
                           min_obs = min_obs, min_events = min_events,
                           bounds = bounds,
                           times = times, sl_fn = sl_fn,
                           verbose = verbose)

    if (verbose) .vmsg("Step %d/%d: Fitting observation model (g_R)...",
                        cur_step + 3L, n_steps)
    obj <- fit_observation(obj, regime = regime_names, covariates = covariates,
                             learners = learners_observation, adaptive_cv = adaptive_cv,
                             min_obs = min_obs, min_events = min_events,
                             bounds = bounds,
                             times = times, sl_fn = sl_fn,
                             verbose = verbose)
    cur_step <- cur_step + 3L
  }

  # --- IPW-specific: weights + estimate ---
  if (do_ipw) {
    if (verbose) .vmsg("Step %d/%d: Computing weights and estimating (IPW)...",
                        cur_step + 1L, n_steps)
    obj <- compute_weights(obj, regime = regime_names,
                             stabilized = stabilized,
                             truncation = truncation,
                             truncation_quantile = truncation_quantile,
                             g_bounds = g_bounds)

    obj <- estimate_ipw(obj, regime = regime_names, times = times,
                        inference = inference, ci_level = ci_level,
                        n_boot = n_boot, cluster = cluster,
                        g_bounds = g_bounds)
    cur_step <- cur_step + 1L
  }

  # --- Shared outcome model (G-comp and/or TMLE) ---
  if (need_outcome) {
    if (verbose) .vmsg("Step %d/%d: Fitting outcome model...",
                        cur_step + 1L, n_steps)
    obj <- fit_outcome(obj, regime = regime_names, covariates = covariates,
                         learners = learners_outcome, adaptive_cv = adaptive_cv,
                         min_obs = min_obs, bounds = bounds,
                         times = times, sl_fn = sl_fn,
                         verbose = verbose)
    cur_step <- cur_step + 1L
  }

  # --- G-comp estimate ---
  if (do_gcomp) {
    if (verbose) .vmsg("Step %d/%d: Estimating (G-comp)...",
                        cur_step + 1L, n_steps)
    obj <- estimate_gcomp(obj, regime = regime_names, times = times,
                          ci_level = ci_level, n_boot = n_boot,
                          verbose = verbose)
    cur_step <- cur_step + 1L
  }

  # --- TMLE estimate ---
  if (do_tmle) {
    # Determine TMLE inference
    tmle_inf <- if (n_boot == 0L) "eif" else inference
    if (tmle_inf %in% c("none", "ic")) tmle_inf <- "eif"

    if (verbose) .vmsg("Step %d/%d: Estimating (TMLE)...",
                        cur_step + 1L, n_steps)
    obj <- estimate_tmle(obj, regime = regime_names, times = times,
                         inference = tmle_inf, ci_level = ci_level,
                         n_boot = n_boot, g_bounds = g_bounds,
                         outcome_range = outcome_range,
                         verbose = verbose)
    cur_step <- cur_step + 1L
  }

  obj
}

#' Extract Results from a longy_data Object
#'
#' Convenience accessor to filter accumulated results by regime and/or
#' estimator.
#'
#' @param obj A \code{longy_data} object with results.
#' @param regime Character vector. Filter to these regime(s). NULL = all.
#' @param estimator Character vector. Filter to these estimator(s)
#'   (e.g. \code{"ipw"}, \code{"gcomp"}, \code{"tmle"}). NULL = all.
#'
#' @return A named list of \code{longy_result} objects matching the filters.
#'   If no results match, returns an empty list.
#'
#' @export
results <- function(obj, regime = NULL, estimator = NULL) {
  obj <- .as_longy_data(obj)
  res <- obj$results
  if (length(res) == 0) return(res)

  if (!is.null(regime)) {
    res <- Filter(function(r) r$regime %in% regime, res)
  }
  if (!is.null(estimator)) {
    res <- Filter(function(r) {
      est <- if (!is.null(r$estimator)) r$estimator else "ipw"
      est %in% estimator
    }, res)
  }
  res
}

#' Plot longy_data Results Across Regimes and Estimators
#'
#' When results are present, creates a plot comparing estimates over time.
#' When multiple estimators are present (e.g., from \code{estimator = "all"}),
#' results are faceted by estimator with regimes shown as colored lines
#' within each panel.
#'
#' @param x A \code{longy_data} object with results.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL (base plot).
#' @export
plot.longy_data <- function(x, ...) {
  if (length(x$results) == 0) {
    message("No results to plot. Run an estimator first.")
    return(invisible(NULL))
  }

  # Combine estimates, extracting regime and estimator from each result
  all_est <- lapply(names(x$results), function(rname) {
    res <- x$results[[rname]]
    est <- as.data.frame(res$estimates)

    # Derive estimator label
    est_type <- if (!is.null(res$estimator)) res$estimator else "ipw"
    est_label <- switch(est_type,
                        gcomp = "G-comp", tmle = "TMLE", "IPW")

    est$estimator_label <- est_label
    est$regime <- res$regime
    est
  })
  combined <- do.call(rbind, all_est)

  n_estimators <- length(unique(combined$estimator_label))
  use_facets <- n_estimators > 1

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
        title = if (use_facets) "Estimates by Estimator and Regime"
                else "Estimates by Regime"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    if (use_facets) {
      p <- p + ggplot2::facet_wrap(~estimator_label, scales = "free_y")
    }

    return(p)
  }

  # Base R fallback
  unique_regimes <- unique(combined$regime)
  cols <- seq_along(unique_regimes)
  names(cols) <- unique_regimes

  if (use_facets) {
    estimators <- unique(combined$estimator_label)
    n_est <- length(estimators)
    old_par <- graphics::par(mfrow = c(1, n_est))
    on.exit(graphics::par(old_par), add = TRUE)
    for (est_lab in estimators) {
      sub <- combined[combined$estimator_label == est_lab, ]
      first <- TRUE
      for (rg in unique_regimes) {
        est <- sub[sub$regime == rg, ]
        if (nrow(est) == 0) next
        if (first) {
          plot(est$time, est$estimate, type = "b", pch = 19, col = cols[rg],
               xlab = "Time", ylab = "Estimate", ylim = y_range,
               main = est_lab)
          first <- FALSE
        } else {
          graphics::lines(est$time, est$estimate, type = "b", pch = 19,
                          col = cols[rg])
        }
        if ("ci_lower" %in% names(est) && "ci_upper" %in% names(est)) {
          graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                           angle = 90, code = 3, length = 0.05, col = cols[rg])
        }
      }
    }
    graphics::legend("topleft", legend = unique_regimes, col = cols, lty = 1,
                     pch = 19, cex = 0.8)
  } else {
    first <- TRUE
    for (rg in unique_regimes) {
      est <- combined[combined$regime == rg, ]
      if (first) {
        plot(est$time, est$estimate, type = "b", pch = 19, col = cols[rg],
             xlab = "Time", ylab = "Estimate", ylim = y_range,
             main = "Estimates by Regime")
        first <- FALSE
      } else {
        graphics::lines(est$time, est$estimate, type = "b", pch = 19,
                        col = cols[rg])
      }
      if ("ci_lower" %in% names(est) && "ci_upper" %in% names(est)) {
        graphics::arrows(est$time, est$ci_lower, est$time, est$ci_upper,
                         angle = 90, code = 3, length = 0.05, col = cols[rg])
      }
    }
    graphics::legend("topleft", legend = unique_regimes, col = cols, lty = 1,
                     pch = 19)
  }
  invisible(NULL)
}
