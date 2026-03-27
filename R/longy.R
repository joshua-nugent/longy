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
#'     \item Scalar 0/1: constant static regime, e.g., \code{list(always = 1L, never = 0L)}
#'     \item Named numeric vector of 0/1: time-varying static regime, e.g.,
#'       \code{list(early = c("0" = 1, "1" = 1, "2" = 0))}
#'     \item Character string: column name with pre-computed 0/1 counterfactual
#'       treatment values, e.g., \code{list(custom = "A_counterfactual")}
#'     \item Function: dynamic regime returning 0/1 per row
#'   }
#'   For static regimes defined as a function of time, use
#'   \code{\link{define_regime}()} directly with \code{static = function(t) ...}.
#' @param estimator Character vector. Which estimator(s) to run. Any combination
#'   of \code{"ipw"}, \code{"gcomp"}, and \code{"tmle"} (e.g.
#'   \code{c("ipw", "tmle")}). Results are keyed as
#'   \code{{regime}_{estimator}}. Nuisance models are shared across estimators
#'   to avoid redundant fitting.
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
#' @param cluster Character. Column name for a clustering/grouping variable.
#'   When specified, cluster-aware CV folds are used in SuperLearner and
#'   cross-fitting, and cluster-robust SEs are computed for IPW inference.
#'   See \code{\link{longy_data}} for details. NULL if no clustering.
#' @param times Numeric vector. Time points for estimation. NULL = all.
#' @param sl_control List. Additional arguments passed to SuperLearner in all
#'   nuisance models. Elements named \code{cvControl} are merged with the
#'   default \code{cvControl} (e.g. \code{list(cvControl = list(stratifyCV = TRUE))}).
#'   \code{cvControl$V} sets the number of CV folds when \code{adaptive_cv = FALSE};
#'   specifying \code{V} with \code{adaptive_cv = TRUE} is an error.
#'   Other elements (e.g. \code{method = "method.CC_LS"}) are passed directly
#'   to SuperLearner. Default \code{list()}.
#' @param adaptive_cv Logical. Adaptive CV fold selection. When TRUE (default),
#'   the number of CV folds is chosen automatically based on effective sample
#'   size. When FALSE, uses \code{sl_control$cvControl$V} if specified, or 10.
#'   Cannot be TRUE when \code{sl_control$cvControl$V} is specified.
#' @param min_obs Integer. Minimum observations for model fitting.
#' @param min_events Integer. Minimum minority-class events required to fit a
#'   model. When the minority class count is below this AND the minority rate
#'   is below 0.01, marginal fallback is used. Default 20.
#' @param bounds Numeric(2). Bounds for the unstabilized cumulative
#'   treatment-censoring product (cumprod of g_a * g_c) in
#'   \code{compute_weights()}. Prevents extreme IPW weights from
#'   near-positivity violations. Default \code{c(0.01, 1)}. Set to NULL to
#'   skip bounding. For TMLE, use \code{g_bounds} instead.
#' @param stabilization Character. Controls the numerator of stabilized IPW
#'   weights. \code{"marginal"} (default) uses unconditional marginal rates.
#'   \code{"baseline"} fits models using only baseline covariates, which can
#'   reduce weight variability. Ignored when \code{stabilized = FALSE} or for
#'   TMLE/G-comp estimators.
#' @param numerator_learners Character vector. SuperLearner library for
#'   baseline numerator models when \code{stabilization = "baseline"}.
#'   Default \code{NULL} uses \code{"SL.glm"}.
#' @param risk_set_treatment Character. Which subjects form the risk set for
#'   treatment model fitting. \code{"all"} (default) uses all uncensored
#'   subjects (fit once, share across regimes). \code{"followers"} restricts to
#'   regime-followers at each time, fitting separate models per regime.
#' @param risk_set_outcome Character. Which subjects form the training set for
#'   outcome models (used by G-comp and TMLE). \code{"all"} (default) uses all
#'   uncensored subjects. \code{"followers"} restricts to regime-followers, so
#'   the model learns \code{E(Q|H)} without extrapolation to non-followers.
#' @param g_bounds Numeric(2). Bounds for the unstabilized cumulative
#'   propensity score used in the TMLE clever covariate denominator. Default
#'   \code{c(0.01, 1)}. Only affects TMLE. For IPW weight control, use
#'   \code{truncation} or \code{truncation_quantile} instead.
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
#' @param use_ffSL Logical. If TRUE, use future-factorial SuperLearner
#'   (parallelizes fold x algorithm combinations via \code{future.apply}).
#'   Default FALSE. Forced to FALSE inside parallel workers to prevent
#'   nested parallelism.
#' @param contrast Logical. If TRUE and \code{length(regimes) >= 2}, auto-compute
#'   pairwise risk-difference contrasts (delta-method SEs when ICs are available).
#'   Stored in \code{obj$contrasts}. Default FALSE.
#' @param parallel Logical. If TRUE and a non-sequential \code{future::plan()}
#'   is active, dispatches time-point models in parallel via
#'   \code{future.apply::future_lapply()}. Default FALSE.
#' @param k Integer or \code{Inf}. Number of lagged time steps of covariate
#'   history to include as additional predictors. See \code{\link{longy_data}}
#'   for details. Default \code{0} (no lag columns).
#' @param impute_tv Logical. LOCF imputation of time-varying covariates with
#'   missingness indicators. See \code{\link{longy_data}}. Default TRUE.
#' @param verbose Logical. Print progress.
#'
#' @return An S3 object of class \code{"longy_data"} with estimation results
#'   accumulated in \code{obj$results} (a named list of \code{longy_result}
#'   objects keyed as \code{{regime}_{estimator}}). Use \code{results(obj)} to
#'   access results, or \code{obj$results$always_ipw$estimates} directly.
#'
#' @examples
#' \dontrun{
#' # --- Basic: always vs. never treat (IPW, default) ---
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C", observation = "R",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, never = 0L)
#' )
#' obj$results$always_ipw$estimates
#' results(obj, estimator = "ipw")
#'
#' # --- Time-varying static regimes ---
#' # Compare always-treat to a "treat for first 6 months" policy.
#' # Names in the vector must match the time values in your data.
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(
#'     always     = 1L,
#'     never      = 0L,
#'     first_six  = c("0"=1, "1"=1, "2"=1, "3"=1, "4"=1, "5"=1,
#'                     "6"=0, "7"=0, "8"=0, "9"=0, "10"=0, "11"=0)
#'   ),
#'   estimator = "ipw"
#' )
#' results(obj)
#'
#' # --- Multiple estimators ---
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, never = 0L),
#'   estimator = c("ipw", "tmle")
#' )
#' results(obj, regime = "always", estimator = "tmle")
#'
#' # --- Pre-computed (shifted) regime ---
#' # Compute your own counterfactual treatment column, then pass the name.
#' # Useful for complex rules involving history, multiple covariates, etc.
#' sim_longy$A_cf <- ifelse(sim_longy$time <= 5 & sim_longy$L1 > 0, 1L, 0L)
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, custom = "A_cf"),
#'   estimator = "ipw"
#' )
#'
#' # --- Time-varying static via the modular pipeline ---
#' # Use define_regime() directly for function-based static regimes:
#' obj <- longy_data(sim_longy, id = "id", time = "time",
#'                   outcome = "Y", treatment = "A") |>
#'   define_regime("always", static = 1L) |>
#'   define_regime("first_6mo", static = function(t) as.integer(t <= 5)) |>
#'   fit_treatment() |>
#'   compute_weights() |>
#'   estimate_ipw()
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
                  bounds = c(0.01, 1),
                  stabilization = "marginal",
                  numerator_learners = NULL,
                  risk_set_treatment = "all",
                  risk_set_outcome = "all",
                  g_bounds = c(0.01, 1),
                  outcome_range = NULL,
                  cross_fit = NULL,
                  cross_fit_seed = NULL,
                  sl_control = list(),
                  use_ffSL = FALSE,
                  contrast = FALSE,
                  parallel = FALSE,
                  k = 0,
                  impute_tv = TRUE,
                  verbose = TRUE) {

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

  if (adaptive_cv && !is.null(sl_control$cvControl$V))
    stop("Cannot specify sl_control$cvControl$V when adaptive_cv=TRUE. ",
         "Set adaptive_cv=FALSE to use a fixed number of CV folds.",
         call. = FALSE)

  # Backward compatibility: "both" -> ipw + gcomp, "all" -> ipw + gcomp + tmle
  if (length(estimator) == 1 && estimator == "both") {
    estimator <- c("ipw", "gcomp")
  } else if (length(estimator) == 1 && estimator == "all") {
    estimator <- c("ipw", "gcomp", "tmle")
  }
  estimator <- match.arg(estimator, c("ipw", "gcomp", "tmle"), several.ok = TRUE)

  # IPW is not valid with competing risks — drop it with a warning
  if ("ipw" %in% estimator && !is.null(competing)) {
    estimator <- setdiff(estimator, "ipw")
    if (length(estimator) == 0) {
      stop("IPW is not supported with competing risks. ",
           "Use estimator = \"gcomp\", \"tmle\", or \"all\" instead. ",
           "See ?estimate_ipw for details.",
           call. = FALSE)
    }
    warning("IPW is not supported with competing risks and was removed from ",
            "the estimator set. Using: ", paste(estimator, collapse = ", "), ". ",
            "See ?estimate_ipw for details.",
            call. = FALSE)
  }

  do_ipw <- "ipw" %in% estimator
  do_gcomp <- "gcomp" %in% estimator
  do_tmle <- "tmle" %in% estimator
  multi <- sum(do_ipw, do_gcomp, do_tmle) > 1

  # n_boot = 0 only disables bootstrap inference; analytic methods (ic, sandwich,

  # eif) are unaffected since they don't need resampling
  if (n_boot == 0L && inference == "bootstrap") inference <- "none"

  # Determine total steps for verbose messaging
  # Shared fits: g_A + g_C fitted once if IPW or TMLE need them
  # g_R shared for IPW and TMLE; outcome fitted once if G-comp or TMLE need it
  n_steps <- 3L  # data + regimes + unadjusted (always)
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
    cluster = cluster, sampling_weights = sampling_weights,
    outcome_type = outcome_type, competing = competing,
    k = k, impute_tv = impute_tv, verbose = verbose
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
    if (is.numeric(rval) && length(rval) == 1 && is.null(names(rval))) {
      # Scalar 0/1
      obj <- define_regime(obj, name = rname, static = as.integer(rval))
    } else if (is.numeric(rval) && !is.null(names(rval))) {
      # Named vector — time-varying static
      obj <- define_regime(obj, name = rname, static = rval)
    } else if (is.character(rval) && length(rval) == 1) {
      # Character string — shifted column name
      obj <- define_regime(obj, name = rname, shifted = rval)
    } else if (is.function(rval)) {
      # Could be a dynamic regime or a static function of time.
      # Use dynamic by default; users wanting static-function should use
      # define_regime() directly or wrap in list(static = fn).
      obj <- define_regime(obj, name = rname, dynamic = rval)
    } else {
      stop(sprintf(
        "Regime '%s': must be a scalar 0/1, a named numeric vector of 0/1, a column name, or a function.",
        rname), call. = FALSE)
    }
  }

  # Fit all regimes at once, then estimate
  regime_names <- names(regimes)

  # --- Unadjusted (always runs, no models needed) ---
  step <- step + 1L
  if (verbose) .vmsg("Step %d/%d: Computing unadjusted estimates...", step, n_steps)
  obj <- estimate_unadjusted(obj, regime = regime_names, times = times,
                              ci_level = ci_level)

  cur_step <- step

  # --- Shared nuisance models: g_A + g_C + g_R (IPW and/or TMLE) ---
  if (need_g) {
    if (verbose) .vmsg("Step %d/%d: Fitting treatment model (g_A)...",
                        cur_step + 1L, n_steps)
    obj <- fit_treatment(obj, regime = regime_names, covariates = covariates,
                           learners = learners_treatment, sl_control = sl_control,
                           adaptive_cv = adaptive_cv,
                           min_obs = min_obs, min_events = min_events,
                           times = times, use_ffSL = use_ffSL,
                           parallel = parallel,
                           verbose = verbose, risk_set = risk_set_treatment)

    if (verbose) .vmsg("Step %d/%d: Fitting censoring model (g_C)...",
                        cur_step + 2L, n_steps)
    obj <- fit_censoring(obj, regime = regime_names, covariates = covariates,
                           learners = learners_censoring, sl_control = sl_control,
                           adaptive_cv = adaptive_cv,
                           min_obs = min_obs, min_events = min_events,
                           times = times, use_ffSL = use_ffSL,
                           parallel = parallel,
                           verbose = verbose)

    if (verbose) .vmsg("Step %d/%d: Fitting observation model (g_R)...",
                        cur_step + 3L, n_steps)
    obj <- fit_observation(obj, regime = regime_names, covariates = covariates,
                             learners = learners_observation, sl_control = sl_control,
                             adaptive_cv = adaptive_cv,
                             min_obs = min_obs, min_events = min_events,
                             times = times, use_ffSL = use_ffSL,
                             parallel = parallel,
                             verbose = verbose)
    cur_step <- cur_step + 3L
  }

  # --- IPW-specific: weights + estimate ---
  if (do_ipw) {
    if (verbose) .vmsg("Step %d/%d: Computing weights and estimating (IPW)...",
                        cur_step + 1L, n_steps)
    obj <- compute_weights(obj, regime = regime_names,
                             stabilized = stabilized,
                             stabilization = stabilization,
                             numerator_learners = numerator_learners,
                             bounds = bounds,
                             truncation = truncation,
                             truncation_quantile = truncation_quantile,
                             verbose = verbose)

    obj <- estimate_ipw(obj, regime = regime_names, times = times,
                        inference = inference, ci_level = ci_level,
                        n_boot = n_boot, cluster = cluster)
    cur_step <- cur_step + 1L
  }

  # --- Shared outcome model (G-comp and/or TMLE) ---
  if (need_outcome) {
    # G-comp bootstrap conflict: if parallel + gcomp + bootstrap, run
    # fit_outcome sequentially so bootstrap gets the parallel resources
    outcome_parallel <- parallel
    if (do_gcomp && n_boot > 0 && parallel) {
      outcome_parallel <- FALSE
      warning("parallel=TRUE with G-comp bootstrap: fit_outcome() runs sequentially so bootstrap gets parallel resources.",
              call. = FALSE)
    }
    # When cross-fitting is enabled and only TMLE needs the outcome model,
    # skip the expensive backward regression — .cf_estimate_tmle() will redo
    # it with cross-fitting anyway. Only store metadata.
    outcome_metadata_only <- isTRUE(obj$crossfit$enabled) && do_tmle && !do_gcomp

    if (verbose) .vmsg("Step %d/%d: Fitting outcome model%s...",
                        cur_step + 1L, n_steps,
                        if (outcome_metadata_only) " (metadata only)" else "")
    obj <- fit_outcome(obj, regime = regime_names, covariates = covariates,
                         learners = learners_outcome, sl_control = sl_control,
                         adaptive_cv = adaptive_cv,
                         min_obs = min_obs,
                         times = times, use_ffSL = use_ffSL,
                         parallel = outcome_parallel,
                         risk_set = risk_set_outcome,
                         verbose = verbose,
                         metadata_only = outcome_metadata_only)
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
                         risk_set = risk_set_outcome,
                         parallel = parallel,
                         verbose = verbose)
    cur_step <- cur_step + 1L
  }

  # --- Auto-compute pairwise contrasts ---
  if (isTRUE(contrast) && length(regime_names) >= 2) {
    obj$contrasts <- list()
    pairs <- utils::combn(regime_names, 2, simplify = FALSE)
    for (est in estimator) {
      for (pair in pairs) {
        ctr <- tryCatch(
          contrast(obj, regime = pair, estimator = est, ci_level = ci_level),
          error = function(e) NULL
        )
        if (!is.null(ctr)) {
          ctr_key <- sprintf("%s_vs_%s_%s", pair[1], pair[2], est)
          obj$contrasts[[ctr_key]] <- ctr
        }
      }
    }
    if (verbose && length(obj$contrasts) > 0)
      .vmsg("Computed %d pairwise contrast(s).", length(obj$contrasts))
  }

  obj
}

#' Add a Regime to an Existing longy_data Object
#'
#' Defines a new regime and runs estimation, reusing compatible nuisance model
#' fits from previously fitted regimes. This avoids re-fitting treatment,
#' censoring, and observation models when the underlying data and covariates
#' are the same.
#'
#' @param obj A \code{longy_data} object, typically returned by \code{longy()}
#'   or the modular pipeline, with at least one regime already fitted.
#' @param name Character. Name for the new regime (must be unique).
#' @param static,shifted,dynamic,stochastic Regime specification. Exactly one
#'   must be provided. See \code{\link{define_regime}} for details.
#' @param description Character. Optional description for the regime.
#' @param estimator Character vector. Which estimator(s) to run. If NULL
#'   (default), matches the estimators present in existing results.
#' @param reuse_nuisance Logical. If TRUE (default), copies compatible
#'   nuisance fits (g_A, g_C, g_R) from existing regimes instead of
#'   re-fitting. Treatment fits are reused when \code{risk_set} was
#'   \code{"all"}; censoring and observation fits are always reused.
#'   Set to FALSE to force re-fitting everything.
#' @param learners Character vector, named list, or NULL. SuperLearner
#'   library specification, same format as in \code{\link{longy}()}. A
#'   character vector applies the same library to all models; a named list
#'   (keys: \code{treatment}, \code{censoring}, \code{observation},
#'   \code{outcome}, \code{default}) specifies per-model libraries; NULL
#'   uses plain glm. When omitted (not explicitly passed), learners are
#'   inherited from existing fits. When explicitly provided \emph{and} a
#'   donor fit exists, the resolved per-model libraries are validated
#'   against the donor — a mismatch raises an error to prevent silent
#'   inconsistency.
#' @param covariates Character vector or NULL. Predictor columns for nuisance
#'   models. Used when fitting from scratch (no donor available). NULL uses
#'   all baseline + timevarying. Ignored when reusing donor fits.
#' @param stabilized Logical. Use stabilized weights (IPW only). Default TRUE.
#' @param bounds Numeric(2) or NULL. Bounds for the cumulative AC product
#'   (IPW only). Default \code{c(0.01, 1)}.
#' @param truncation Numeric. Weight truncation cap (IPW only).
#' @param truncation_quantile Numeric. Quantile-based weight truncation (IPW only).
#' @param inference Character. Inference method. If NULL (default), matches
#'   existing results.
#' @param ci_level Numeric. Confidence level. Default 0.95.
#' @param n_boot Integer. Bootstrap replicates. Default 200.
#' @param cluster Character. Cluster column for robust SEs (IPW only).
#' @param g_bounds Numeric(2). Bounds for TMLE clever covariate denominator.
#'   Default \code{c(0.01, 1)}.
#' @param outcome_range Numeric(2) or NULL. Range for continuous outcome
#'   scaling in TMLE.
#' @param risk_set_outcome Character. Risk set for outcome models.
#'   Default \code{"all"}.
#' @param times Numeric vector. Time points for estimation. NULL = all.
#' @param parallel Logical. Parallelize where possible. Default FALSE.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @return Modified \code{longy_data} object with the new regime added and
#'   estimation results accumulated in \code{obj$results}.
#'
#' @details
#' \code{add_regime()} is designed for the common workflow where you run
#' \code{longy()} once with initial regimes, then want to add more regimes
#' without re-fitting expensive nuisance models.
#'
#' \strong{What gets reused} (when \code{reuse_nuisance = TRUE}):
#' \itemize{
#'   \item \strong{g_C} (censoring model): always reused — identical across regimes.
#'   \item \strong{g_R} (observation model): always reused — identical across regimes.
#'   \item \strong{g_A} (treatment model): reused when the existing fit used
#'     \code{risk_set = "all"} (the default). If the existing fit used
#'     \code{risk_set = "followers"}, treatment is re-fitted for the new
#'     regime's followers.
#' }
#'
#' \strong{What is always re-fitted}:
#' \itemize{
#'   \item \strong{Outcome models}: regime-specific (counterfactual treatment
#'     differs), always fitted fresh.
#'   \item \strong{Weights}: regime-specific, always computed fresh.
#' }
#'
#' @examples
#' \dontrun{
#' # Initial analysis
#' obj <- longy(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
#'   regimes = list(always = 1L, never = 0L),
#'   estimator = c("ipw", "tmle")
#' )
#'
#' # Later, add a time-varying regime — reuses all nuisance models
#' obj <- add_regime(obj, name = "early_stop",
#'   static = c("0"=1, "1"=1, "2"=0))
#'
#' # Add a pre-computed regime
#' obj$data$A_cf <- ifelse(obj$data$time <= 1, 1L, 0L)
#' obj <- add_regime(obj, name = "first_two", shifted = "A_cf")
#'
#' # Compare all regimes
#' results(obj)
#' }
#'
#' @export
add_regime <- function(obj, name, static = NULL, shifted = NULL,
                       dynamic = NULL, stochastic = NULL,
                       description = NULL,
                       estimator = NULL,
                       reuse_nuisance = TRUE,
                       learners,
                       covariates = NULL,
                       stabilized = TRUE,
                       stabilization = "marginal",
                       numerator_learners = NULL,
                       bounds = c(0.01, 1),
                       truncation = NULL,
                       truncation_quantile = NULL,
                       inference = NULL,
                       ci_level = 0.95,
                       n_boot = 200L,
                       cluster = NULL,
                       g_bounds = c(0.01, 1),
                       outcome_range = NULL,
                       risk_set_outcome = "all",
                       times = NULL,
                       parallel = FALSE,
                       verbose = TRUE) {

  obj <- .as_longy_data(obj)

  # --- Learners: detect whether user explicitly provided them ---
  learners_provided <- !missing(learners)
  if (!learners_provided) learners <- NULL

  # Resolve per-model learner libraries (same logic as longy())
  # Used for validation against donor fits and for from-scratch fitting
  if (learners_provided) {
    lrn_treatment   <- .resolve_learners(learners, "treatment")
    lrn_censoring   <- .resolve_learners(learners, "censoring")
    lrn_observation <- .resolve_learners(learners, "observation")
    lrn_outcome     <- .resolve_learners(learners, "outcome")
  }

  # --- Define the regime ---
  obj <- define_regime(obj, name = name, static = static, shifted = shifted,
                       dynamic = dynamic, stochastic = stochastic,
                       description = description)

  # --- Detect existing estimators if not specified ---
  if (is.null(estimator)) {
    existing_estimators <- unique(vapply(obj$results, function(r) {
      if (!is.null(r$estimator)) r$estimator else "ipw"
    }, character(1)))
    # Exclude "unadjusted" from auto-detection — it always runs
    estimator <- setdiff(existing_estimators, "unadjusted")
    if (length(estimator) == 0) {
      estimator <- "ipw"
      if (verbose) .vmsg("  No existing estimators detected, defaulting to IPW")
    }
  }
  estimator <- match.arg(estimator, c("ipw", "gcomp", "tmle"), several.ok = TRUE)

  # IPW not valid with competing risks
  if ("ipw" %in% estimator && !is.null(obj$nodes$competing)) {
    estimator <- setdiff(estimator, "ipw")
    if (length(estimator) == 0) {
      stop("IPW is not supported with competing risks. ",
           "Use estimator = \"gcomp\" or \"tmle\" instead.",
           call. = FALSE)
    }
    warning("IPW is not supported with competing risks and was removed. ",
            "Using: ", paste(estimator, collapse = ", "), ".",
            call. = FALSE)
  }

  do_ipw <- "ipw" %in% estimator
  do_gcomp <- "gcomp" %in% estimator
  do_tmle <- "tmle" %in% estimator
  need_g <- do_ipw || do_tmle
  need_outcome <- do_gcomp || do_tmle

  # --- Detect existing inference method if not specified ---
  # IPW accepts: ic, bootstrap, sandwich, none
  # TMLE accepts: eif, bootstrap, none
  # Detect separately per estimator type to avoid cross-contamination
  if (is.null(inference)) {
    # Look for IPW-type inference first, then TMLE-type
    ipw_results <- Filter(function(r) identical(r$estimator, "ipw"), obj$results)
    tmle_results <- Filter(function(r) identical(r$estimator, "tmle"), obj$results)
    if (length(ipw_results) > 0) {
      ipw_inf <- ipw_results[[1]]$inference
      if (is.null(ipw_inf)) ipw_inf <- "ic"
    } else {
      ipw_inf <- "ic"
    }
    if (length(tmle_results) > 0) {
      tmle_inf <- tmle_results[[1]]$inference
      if (is.null(tmle_inf)) tmle_inf <- "eif"
    } else {
      tmle_inf <- "eif"
    }
    inference <- ipw_inf
  } else {
    ipw_inf <- inference
    tmle_inf <- if (inference %in% c("ic", "sandwich")) "eif" else inference
  }

  # --- Find a donor regime for reusing fits ---
  donor <- NULL
  if (reuse_nuisance) {
    available <- names(obj$fits$treatment)
    if (length(available) == 0) available <- names(obj$fits$censoring)
    if (length(available) == 0) available <- names(obj$fits$observation)
    if (length(available) > 0) donor <- available[1]
  }

  # --- Validate learners against donor fits if user explicitly provided them ---
  if (learners_provided && !is.null(donor)) {
    # Treatment
    donor_trt <- obj$fits$treatment[[donor]]
    if (!is.null(donor_trt) && !identical(donor_trt$learners, lrn_treatment)) {
      stop(sprintf(
        "Learner mismatch for treatment model: existing fit uses [%s] but learners specifies [%s]. ",
        paste(donor_trt$learners %||% "NULL (glm)", collapse = ", "),
        paste(lrn_treatment %||% "NULL (glm)", collapse = ", ")),
        "Use the same learners as the original longy() call, or set reuse_nuisance=FALSE to refit.",
        call. = FALSE)
    }
    # Censoring (check first cause's learners)
    donor_cens <- obj$fits$censoring[[donor]]
    if (!is.null(donor_cens) && length(donor_cens) > 0) {
      donor_cens_lrn <- donor_cens[[1]]$learners
      if (!identical(donor_cens_lrn, lrn_censoring)) {
        stop(sprintf(
          "Learner mismatch for censoring model: existing fit uses [%s] but learners specifies [%s]. ",
          paste(donor_cens_lrn %||% "NULL (glm)", collapse = ", "),
          paste(lrn_censoring %||% "NULL (glm)", collapse = ", ")),
          "Use the same learners as the original longy() call, or set reuse_nuisance=FALSE to refit.",
          call. = FALSE)
      }
    }
    # Observation
    donor_obs <- obj$fits$observation[[donor]]
    if (!is.null(donor_obs)) {
      if (!identical(donor_obs$learners, lrn_observation)) {
        stop(sprintf(
          "Learner mismatch for observation model: existing fit uses [%s] but learners specifies [%s]. ",
          paste(donor_obs$learners %||% "NULL (glm)", collapse = ", "),
          paste(lrn_observation %||% "NULL (glm)", collapse = ", ")),
          "Use the same learners as the original longy() call, or set reuse_nuisance=FALSE to refit.",
          call. = FALSE)
      }
    }
    # Outcome (if an outcome donor exists)
    donor_out_check <- NULL
    for (rn in names(obj$fits$outcome)) {
      if (!is.null(obj$fits$outcome[[rn]])) {
        donor_out_check <- obj$fits$outcome[[rn]]
        break
      }
    }
    if (!is.null(donor_out_check) &&
        !identical(donor_out_check$learners, lrn_outcome)) {
      stop(sprintf(
        "Learner mismatch for outcome model: existing fit uses [%s] but learners specifies [%s]. ",
        paste(donor_out_check$learners %||% "NULL (glm)", collapse = ", "),
        paste(lrn_outcome %||% "NULL (glm)", collapse = ", ")),
        "Use the same learners as the original longy() call, or set reuse_nuisance=FALSE to refit.",
        call. = FALSE)
    }
  }

  # --- Unadjusted (always) ---
  if (verbose) .vmsg("add_regime '%s': computing unadjusted estimate...", name)
  obj <- estimate_unadjusted(obj, regime = name, times = times,
                              ci_level = ci_level)

  # --- Nuisance models: g_A, g_C, g_R ---
  if (need_g) {

    # Treatment (g_A)
    if (!is.null(donor) && reuse_nuisance) {
      donor_trt <- obj$fits$treatment[[donor]]
      if (!is.null(donor_trt) && identical(donor_trt$risk_set, "all")) {
        if (verbose) .vmsg("  Reusing treatment model from regime '%s'", donor)
        obj$fits$treatment[[name]] <- donor_trt
        obj$fits$treatment[[name]]$regime <- name
      } else {
        # risk_set was "followers" — must refit using donor's params
        if (verbose) .vmsg("  Fitting treatment model (risk_set = 'followers')...")
        covs <- if (!is.null(donor_trt)) donor_trt$covariates else covariates
        lrn <- if (!is.null(donor_trt)) donor_trt$learners else
                 if (learners_provided) lrn_treatment else NULL
        obj <- fit_treatment(obj, regime = name, covariates = covs,
                             learners = lrn,
                             sl_control = if (!is.null(donor_trt)) donor_trt$sl_control else list(),
                             adaptive_cv = if (!is.null(donor_trt)) donor_trt$adaptive_cv else TRUE,
                             use_ffSL = if (!is.null(donor_trt)) isTRUE(donor_trt$use_ffSL) else FALSE,
                             times = times, parallel = parallel,
                             verbose = verbose, risk_set = "followers")
      }
    } else {
      # No donor — fit from scratch using user-provided learners
      if (verbose) .vmsg("  No existing fits to reuse, fitting treatment model...")
      lrn <- if (learners_provided) lrn_treatment else NULL
      obj <- fit_treatment(obj, regime = name, covariates = covariates,
                           learners = lrn, times = times,
                           parallel = parallel, verbose = verbose)
    }

    # Censoring (g_C) — always regime-independent
    if (!is.null(donor) && reuse_nuisance &&
        !is.null(obj$fits$censoring[[donor]])) {
      if (verbose) .vmsg("  Reusing censoring model from regime '%s'", donor)
      obj$fits$censoring[[name]] <- obj$fits$censoring[[donor]]
    } else {
      if (verbose) .vmsg("  Fitting censoring model...")
      lrn <- if (learners_provided) lrn_censoring else NULL
      obj <- fit_censoring(obj, regime = name, covariates = covariates,
                           learners = lrn, times = times,
                           parallel = parallel, verbose = verbose)
    }

    # Observation (g_R) — always regime-independent
    if (!is.null(donor) && reuse_nuisance &&
        !is.null(obj$fits$observation[[donor]])) {
      if (verbose) .vmsg("  Reusing observation model from regime '%s'", donor)
      obj$fits$observation[[name]] <- obj$fits$observation[[donor]]
      obj$fits$observation[[name]]$regime <- name
    } else {
      if (verbose) .vmsg("  Fitting observation model...")
      lrn <- if (learners_provided) lrn_observation else NULL
      obj <- fit_observation(obj, regime = name, covariates = covariates,
                             learners = lrn, times = times,
                             parallel = parallel, verbose = verbose)
    }
  }

  # --- IPW: weights + estimate ---
  if (do_ipw) {
    if (verbose) .vmsg("  Computing weights and estimating (IPW)...")
    obj <- compute_weights(obj, regime = name, stabilized = stabilized,
                           stabilization = stabilization,
                           numerator_learners = numerator_learners,
                           bounds = bounds,
                           truncation = truncation,
                           truncation_quantile = truncation_quantile,
                           verbose = verbose)
    obj <- estimate_ipw(obj, regime = name, times = times,
                        inference = ipw_inf, ci_level = ci_level,
                        n_boot = n_boot, cluster = cluster)
  }

  # --- Outcome model (G-comp / TMLE) ---
  if (need_outcome) {
    # Extract fitting params from an existing outcome fit if available
    donor_out <- NULL
    for (rn in names(obj$fits$outcome)) {
      if (!is.null(obj$fits$outcome[[rn]])) {
        donor_out <- obj$fits$outcome[[rn]]
        break
      }
    }
    out_covs <- if (!is.null(donor_out)) donor_out$covariates else covariates
    out_lrn <- if (!is.null(donor_out)) donor_out$learners else
                 if (learners_provided) lrn_outcome else NULL
    out_bnd <- if (!is.null(donor_out)) donor_out$bounds else c(0.005, 0.995)
    out_sl <- if (!is.null(donor_out)) donor_out$sl_control else list()
    out_acv <- if (!is.null(donor_out)) isTRUE(donor_out$adaptive_cv) else TRUE
    out_ffsl <- if (!is.null(donor_out)) isTRUE(donor_out$use_ffSL) else FALSE

    outcome_metadata_only <- isTRUE(obj$crossfit$enabled) && do_tmle && !do_gcomp

    if (verbose) .vmsg("  Fitting outcome model%s...",
                        if (outcome_metadata_only) " (metadata only)" else "")
    obj <- fit_outcome(obj, regime = name, covariates = out_covs,
                       learners = out_lrn, bounds = out_bnd,
                       sl_control = out_sl, adaptive_cv = out_acv,
                       times = times, use_ffSL = out_ffsl,
                       parallel = parallel, risk_set = risk_set_outcome,
                       verbose = verbose, metadata_only = outcome_metadata_only)
  }

  # --- G-comp estimate ---
  if (do_gcomp) {
    if (verbose) .vmsg("  Estimating (G-comp)...")
    obj <- estimate_gcomp(obj, regime = name, times = times,
                          ci_level = ci_level, n_boot = n_boot,
                          verbose = verbose)
  }

  # --- TMLE estimate ---
  if (do_tmle) {
    if (n_boot == 0L) tmle_inf <- "eif"
    if (tmle_inf %in% c("none", "ic", "sandwich")) tmle_inf <- "eif"

    if (verbose) .vmsg("  Estimating (TMLE)...")
    obj <- estimate_tmle(obj, regime = name, times = times,
                         inference = tmle_inf, ci_level = ci_level,
                         n_boot = n_boot, g_bounds = g_bounds,
                         outcome_range = outcome_range,
                         risk_set = risk_set_outcome,
                         parallel = parallel, verbose = verbose)
  }

  if (verbose) .vmsg("add_regime '%s': done.", name)
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
                        gcomp = "G-comp", tmle = "TMLE",
                        unadjusted = "Unadjusted", "IPW")

    est$estimator_label <- est_label
    est$regime <- res$regime
    est
  })
  combined <- as.data.frame(data.table::rbindlist(all_est, fill = TRUE))

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
