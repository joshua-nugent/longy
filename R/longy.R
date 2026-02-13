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
