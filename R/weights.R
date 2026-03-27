#' Enforce Cumulative Regime Consistency
#'
#' For static regimes, a subject can only contribute at time t if they followed
#' the regime at all prior times 0..t. After filtering a predictions table to
#' regime-followers, some subjects may have gaps (they deviated at some time
#' but re-matched later). This function keeps only each subject's contiguous
#' prefix from the first study time.
#'
#' @param dt A data.table with at least \code{id_col} and \code{time_col}.
#' @param id_col Character. Name of the ID column.
#' @param time_col Character. Name of the time column. Default \code{".time"}.
#' @param all_times Numeric vector. The full set of study time values (sorted).
#'
#' @return Subset of \code{dt} with inconsistent rows removed.
#' @noRd
.enforce_regime_consistency <- function(dt, id_col, time_col = ".time",
                                        all_times = NULL) {
  if (is.null(all_times)) all_times <- sort(unique(dt[[time_col]]))
  if (length(all_times) <= 1) return(dt)

  data.table::setkeyv(dt, c(id_col, time_col))

  # Map each time value to its position in the full sequence
  time_rank_map <- data.table::data.table(
    ..tr_time = all_times,
    ..expected_rank = seq_along(all_times)
  )
  data.table::setnames(time_rank_map, "..tr_time", time_col)

  # Merge to get expected rank, then compute actual rank within each subject
  dt <- merge(dt, time_rank_map, by = time_col, all.x = TRUE, sort = FALSE)
  data.table::setkeyv(dt, c(id_col, time_col))
  dt[, ..actual_rank := seq_len(.N), by = c(id_col)]

  # Keep only rows where expected rank matches actual rank (no gaps)
  dt <- dt[..expected_rank == ..actual_rank]
  dt[, c("..expected_rank", "..actual_rank") := NULL]

  dt
}

#' Compute Unstabilized Cumulative g
#'
#' Shared helper used by both \code{compute_weights()} (IPW) and
#' \code{estimate_tmle()} (TMLE). Computes the unstabilized cumulative
#' probability of following the regime and remaining uncensored through each
#' time point: g_cum(s) = prod_{j=0}^{s} g_a(j) * g_c(j).
#'
#' @param obj A \code{longy_data} object with treatment (and optionally
#'   censoring) models already fit.
#' @param regime Character. Name of the regime.
#' @param g_bounds Numeric vector of length 2. Bounds for cumulative g after
#'   cumulation. Default \code{c(0.01, 1)}.
#'
#' @return A data.table with columns: (id_col, .time, .g_a, .g_c, .g_point,
#'   .g_cum, .g_r). When no observation model is fit, \code{.g_r = 1}.
#' @noRd
.compute_cumulative_g <- function(obj, regime, g_bounds = c(0.01, 1)) {
  nodes <- obj$nodes
  id_col <- nodes$id
  reg <- obj$regimes[[regime]]


  # --- Treatment component: g_a ---
  gA <- obj$fits$treatment[[regime]]$predictions
  if (reg$type == "static") {
    # Resolve regime value d(t) at each row's time, filter to followers,
    # compute g_a = P(A = d(t) | past)
    gA_work <- data.table::copy(gA)
    gA_work[, .d := .resolve_static_at_time(reg$value, .time)]
    gA_sub <- gA_work[.treatment == .d]
    g_dt <- gA_sub[, c(id_col, ".time", ".p_a", ".d"), with = FALSE]
    g_dt[, .g_a := ifelse(.d == 1L, .p_a, 1 - .p_a)]
    g_dt[, c(".p_a", ".d") := NULL]
  } else {
    # Dynamic/stochastic: g_a = p_a for all (consistency built in)
    g_dt <- data.table::copy(gA[, c(id_col, ".time", ".p_a"), with = FALSE])
    data.table::setnames(g_dt, ".p_a", ".g_a")
  }

  # For static regimes: enforce cumulative consistency — subjects contribute
  # at time t only if they followed the regime at ALL times 0..t
  if (reg$type == "static") {
    g_dt <- .enforce_regime_consistency(g_dt, id_col, ".time", obj$meta$time_values)
  }

  # --- Censoring component: g_c ---
  cens_fits <- obj$fits$censoring[[regime]]
  if (length(cens_fits) > 0) {
    g_c_combined <- NULL
    for (cvar in names(cens_fits)) {
      gC <- cens_fits[[cvar]]$predictions
      # Keep uncensored rows (C=0), p_c = P(C=0|past)
      gC_uncens <- gC[gC$.censored == 0L, ]
      gc_col <- gC_uncens[, c(id_col, ".time", ".p_c"), with = FALSE]
      data.table::setnames(gc_col, ".p_c", paste0(".g_c_", cvar))
      if (is.null(g_c_combined)) {
        g_c_combined <- gc_col
      } else {
        g_c_combined <- merge(g_c_combined, gc_col,
                              by = c(id_col, ".time"), all = FALSE)
      }
    }
    # Product of all censoring sources
    c_cols <- paste0(".g_c_", names(cens_fits))
    g_c_combined[, .g_c := Reduce(`*`, .SD), .SDcols = c_cols]
    g_c_combined <- g_c_combined[, c(id_col, ".time", ".g_c"), with = FALSE]
    g_dt <- merge(g_dt, g_c_combined, by = c(id_col, ".time"), all = FALSE)
  } else {
    g_dt[, .g_c := 1]
  }

  # --- Combine and cumulate ---
  g_dt[, .g_point := .g_a * .g_c]
  data.table::setkeyv(g_dt, c(id_col, ".time"))
  g_dt[, .g_cum := cumprod(.g_point), by = c(id_col)]

  # Count positivity violations before bounding
  n_trunc_lower <- sum(g_dt$.g_cum < g_bounds[1])
  n_total <- nrow(g_dt)
  if (n_trunc_lower > 0) {
    pct <- 100 * n_trunc_lower / n_total
    warning(sprintf(
      "Positivity: %d subject-time obs (%.1f%%) had cumulative g truncated to g_bounds[1]=%.4f.",
      n_trunc_lower, pct, g_bounds[1]), call. = FALSE)
  }

  # Bound cumulative g
  g_dt[, .g_cum := pmax(.g_cum, g_bounds[1])]
  g_dt[, .g_cum := pmin(.g_cum, g_bounds[2])]

  # --- Observation component: g_r (point-in-time, NOT cumulated) ---
  obs_fit <- obj$fits$observation[[regime]]
  if (!is.null(obs_fit)) {
    gR <- obs_fit$predictions
    # Keep observed rows (R=1), p_r = P(R=1|past)
    gR_obs <- gR[gR$.observed == 1L, ]
    gr_col <- gR_obs[, c(id_col, ".time", ".p_r"), with = FALSE]
    data.table::setnames(gr_col, ".p_r", ".g_r")
    g_dt <- merge(g_dt, gr_col, by = c(id_col, ".time"), all.x = TRUE)
    # Subjects not in observation predictions (e.g., no intermittent missingness
    # at that time) get g_r = 1
    g_dt[is.na(.g_r), .g_r := 1]
  } else {
    g_dt[, .g_r := 1]
  }

  g_dt
}

#' Fit Baseline-Only Numerator Models
#'
#' Fits models using only baseline covariates for use as the numerator in
#' baseline-stabilized IPW weights. Called from \code{fit_treatment()},
#' \code{fit_censoring()}, or on-demand from \code{compute_weights()}.
#'
#' The risk set at each time point matches the denominator model exactly:
#' the helper reads existing denominator predictions to identify which
#' (id, time) pairs were at risk, then joins back to \code{obj$data} to
#' obtain baseline covariates.
#'
#' @param obj A \code{longy_data} object with denominator models already fit.
#' @param regime Character. Name of the regime.
#' @param model_type Character. \code{"treatment"} or \code{"censoring"}.
#' @param learners Character vector. SuperLearner library for numerator models.
#'   Default \code{NULL} uses \code{"SL.glm"} only.
#' @param verbose Logical. Print progress.
#'
#' @return Modified \code{obj} with numerator fits stored in
#'   \code{obj$fits$numerator$treatment[[regime]]} or
#'   \code{obj$fits$numerator$censoring[[regime]]}.
#' @noRd
.fit_baseline_numerator <- function(obj, regime, model_type = c("treatment", "censoring"),
                                     learners = NULL, verbose = TRUE) {
  model_type <- match.arg(model_type)
  nodes <- obj$nodes
  id_col <- nodes$id
  baseline_covs <- nodes$baseline

  if (length(baseline_covs) == 0) {
    warning("No baseline covariates available; baseline numerator equals marginal.",
            call. = FALSE)
  }

  if (is.null(learners)) learners <- "SL.glm"

  cf_enabled <- isTRUE(obj$crossfit$enabled)

  if (model_type == "treatment") {
    obj <- .fit_baseline_numerator_treatment(obj, regime, baseline_covs,
                                              learners, cf_enabled, verbose)
  } else {
    obj <- .fit_baseline_numerator_censoring(obj, regime, baseline_covs,
                                              learners, cf_enabled, verbose)
  }
  obj
}

#' @noRd
.fit_baseline_numerator_treatment <- function(obj, regime, baseline_covs,
                                               learners, cf_enabled, verbose) {
  nodes <- obj$nodes
  id_col <- nodes$id
  denom <- obj$fits$treatment[[regime]]$predictions
  dt <- obj$data
  time_vals <- sort(unique(denom$.time))

  results <- vector("list", length(time_vals))

  for (i in seq_along(time_vals)) {
    tt <- time_vals[i]
    # Risk set = subjects in the denominator predictions at this time
    denom_t <- denom[denom$.time == tt, ]
    risk_ids <- denom_t[[id_col]]
    if (length(risk_ids) == 0) next

    dt_t <- dt[dt[[nodes$time]] == tt, ]
    # Join to get baseline covariates for risk-set subjects
    dt_risk <- dt_t[dt_t[[id_col]] %in% risk_ids, ]
    # Ensure same order as denom_t
    dt_risk <- dt_risk[match(risk_ids, dt_risk[[id_col]]), ]

    lag_covs <- .get_lag_covariates(nodes, i)
    all_covs <- c(baseline_covs, lag_covs)
    # Filter to covariates that actually exist in the data
    all_covs <- intersect(all_covs, names(dt_risk))
    Y <- denom_t$.treatment

    if (length(all_covs) == 0 || length(unique(Y)) <= 1) {
      # No covariates or constant outcome: use marginal
      ow <- if (!is.null(nodes$sampling_weights)) dt_risk[[nodes$sampling_weights]] else NULL
      p_num <- rep(if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y),
                   length(Y))
    } else if (cf_enabled) {
      p_num <- .cf_fit_numerator_one_time(
        dt_risk, Y, all_covs, nodes, learners, obj$crossfit, verbose, tt, "g_A_num")
    } else {
      X <- as.data.frame(dt_risk[, all_covs, with = FALSE])
      ow <- if (!is.null(nodes$sampling_weights)) dt_risk[[nodes$sampling_weights]] else NULL
      ctx <- sprintf("g_A_num, time=%d, n=%d", tt, length(Y))
      fit <- .safe_sl(Y = Y, X = X, learners = learners,
                      cv_folds = min(10L, length(Y)), obs_weights = ow,
                      context = ctx, verbose = FALSE)
      p_num <- fit$predictions
    }

    results[[i]] <- data.table::data.table(
      .id = risk_ids,
      .time = tt,
      .p_num = p_num
    )
    data.table::setnames(results[[i]], ".id", id_col)

    if (verbose) .vmsg("  g_A_num time %d: n=%d (baseline numerator)", tt, length(Y))
  }

  results <- data.table::rbindlist(results[!vapply(results, is.null, logical(1))])

  if (is.null(obj$fits$numerator)) obj$fits$numerator <- list()
  if (is.null(obj$fits$numerator$treatment)) obj$fits$numerator$treatment <- list()
  obj$fits$numerator$treatment[[regime]] <- list(
    predictions = results,
    covariates = baseline_covs,
    learners = learners
  )
  obj
}

#' @noRd
.fit_baseline_numerator_censoring <- function(obj, regime, baseline_covs,
                                               learners, cf_enabled, verbose) {
  nodes <- obj$nodes
  id_col <- nodes$id
  cens_fits <- obj$fits$censoring[[regime]]
  dt <- obj$data

  if (length(cens_fits) == 0) return(obj)

  if (is.null(obj$fits$numerator)) obj$fits$numerator <- list()
  if (is.null(obj$fits$numerator$censoring)) obj$fits$numerator$censoring <- list()

  all_cens_results <- list()

  for (cvar in names(cens_fits)) {
    denom <- cens_fits[[cvar]]$predictions
    time_vals <- sort(unique(denom$.time))
    results <- vector("list", length(time_vals))

    for (i in seq_along(time_vals)) {
      tt <- time_vals[i]
      denom_t <- denom[denom$.time == tt, ]
      risk_ids <- denom_t[[id_col]]
      if (length(risk_ids) == 0) next

      dt_t <- dt[dt[[nodes$time]] == tt, ]
      dt_risk <- dt_t[dt_t[[id_col]] %in% risk_ids, ]
      dt_risk <- dt_risk[match(risk_ids, dt_risk[[id_col]]), ]

      lag_covs <- .get_lag_covariates(nodes, i)
      all_covs <- c(baseline_covs, lag_covs)
      all_covs <- intersect(all_covs, names(dt_risk))
      Y <- denom_t$.censored

      if (length(all_covs) == 0 || length(unique(Y)) <= 1) {
        ow <- if (!is.null(nodes$sampling_weights)) dt_risk[[nodes$sampling_weights]] else NULL
        p_num <- rep(if (!is.null(ow)) stats::weighted.mean(Y, ow) else mean(Y),
                     length(Y))
      } else if (cf_enabled) {
        p_num <- .cf_fit_numerator_one_time(
          dt_risk, Y, all_covs, nodes, learners, obj$crossfit, verbose, tt,
          sprintf("g_C_num(%s)", cvar))
      } else {
        X <- as.data.frame(dt_risk[, all_covs, with = FALSE])
        ow <- if (!is.null(nodes$sampling_weights)) dt_risk[[nodes$sampling_weights]] else NULL
        ctx <- sprintf("g_C_num(%s), time=%d, n=%d", cvar, tt, length(Y))
        fit <- .safe_sl(Y = Y, X = X, learners = learners,
                        cv_folds = min(10L, length(Y)), obs_weights = ow,
                        context = ctx, verbose = FALSE)
        p_num <- fit$predictions
      }

      results[[i]] <- data.table::data.table(
        .id = risk_ids,
        .time = tt,
        .p_num = p_num
      )
      data.table::setnames(results[[i]], ".id", id_col)

      if (verbose) .vmsg("  g_C_num(%s) time %d: n=%d (baseline numerator)",
                          cvar, tt, length(Y))
    }

    all_cens_results[[cvar]] <- list(
      predictions = data.table::rbindlist(
        results[!vapply(results, is.null, logical(1))]),
      covariates = baseline_covs,
      learners = learners
    )
  }

  obj$fits$numerator$censoring[[regime]] <- all_cens_results
  obj
}

#' Cross-fitted numerator model for one time point
#' @noRd
.cf_fit_numerator_one_time <- function(dt_risk, Y, covs, nodes, learners,
                                        crossfit_info, verbose, tt, label) {
  fold_col <- crossfit_info$fold_id
  n_folds <- crossfit_info$n_folds
  id_col <- nodes$id
  folds <- dt_risk[[fold_col]]
  X <- as.data.frame(dt_risk[, covs, with = FALSE])
  ow_all <- if (!is.null(nodes$sampling_weights)) dt_risk[[nodes$sampling_weights]] else NULL

  p_num <- rep(NA_real_, length(Y))

  for (k in seq_len(n_folds)) {
    val_idx <- which(folds == k)
    train_idx <- which(folds != k)
    if (length(val_idx) == 0) next

    Y_train <- Y[train_idx]
    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    ow_train <- if (!is.null(ow_all)) ow_all[train_idx] else NULL

    if (length(unique(Y_train)) > 1 && length(train_idx) >= 20L) {
      ctx <- sprintf("CF %s, time=%d, fold=%d/%d, n=%d",
                     label, tt, k, n_folds, length(Y_train))
      fit <- .safe_sl(Y = Y_train, X = X_train, learners = learners,
                      cv_folds = min(10L, length(Y_train)),
                      obs_weights = ow_train, context = ctx, verbose = FALSE)
      p_num[val_idx] <- .predict_from_fit(fit, X_val)
    } else {
      marg_train <- if (!is.null(ow_train)) {
        stats::weighted.mean(Y_train, ow_train)
      } else {
        mean(Y_train)
      }
      p_num[val_idx] <- rep(marg_train, length(val_idx))
    }
  }
  p_num
}

#' Compute Inverse Probability Weights
#'
#' Combines treatment, censoring, and observation model predictions into
#' stabilized inverse probability weights. Treatment and censoring weights
#' are cumulated (absorbing processes); observation weights are point-in-time
#' (intermittent process).
#'
#' The unstabilized cumulative treatment-censoring product
#' (cumprod of g_a * g_c) is bounded to \code{bounds} after cumulation.
#' This directly prevents extreme weights from near-positivity violations
#' that compound over time. Observation weights (intermittent, not cumulated)
#' are applied separately as point-in-time multipliers.
#'
#' Additional weight truncation is available via \code{truncation} (hard cap)
#' and \code{truncation_quantile} (percentile-based cap), which bound the
#' final stabilized weight directly. For TMLE, use \code{g_bounds} in
#' \code{\link{estimate_tmle}()} instead — it bounds the unstabilized
#' cumulative propensity score used in the clever covariate denominator.
#'
#' @param obj A `longy_data` object with treatment (and optionally censoring/observation)
#'   models already fit.
#' @param regime Character. Name of the regime.
#' @param stabilized Logical. Use stabilized weights (numerator = marginal probability).
#' @param stabilization Character. Controls the numerator of stabilized weights.
#'   \code{"marginal"} (default) uses unconditional marginal rates.
#'   \code{"baseline"} fits models using only baseline covariates, which can
#'   reduce weight variability when baseline covariates are strong treatment
#'   predictors. Ignored when \code{stabilized = FALSE}.
#' @param numerator_learners Character vector. SuperLearner library for
#'   baseline numerator models when \code{stabilization = "baseline"}.
#'   Default \code{NULL} uses \code{"SL.glm"} only. Ignored when
#'   \code{stabilization = "marginal"}.
#' @param bounds Numeric(2). Bounds for the unstabilized cumulative
#'   treatment-censoring product (cumprod of g_a * g_c) after cumulation.
#'   Default \code{c(0.01, 1)}. Set to NULL to skip bounding.
#' @param truncation Numeric. Hard upper bound for final weights. NULL for no
#'   truncation.
#' @param truncation_quantile Numeric in (0,1). Truncate final weights at this
#'   quantile. Applied after \code{truncation} if both specified.
#' @param recompute Logical. If FALSE (default), errors when weights are already
#'   computed for the requested regime(s). Set to TRUE to re-compute.
#' @param verbose Logical. Print progress for numerator model fitting.
#'   Default TRUE.
#'
#' @return Modified `longy_data` object with weights stored.
#' @export
compute_weights <- function(obj, regime = NULL, stabilized = TRUE,
                            stabilization = c("marginal", "baseline"),
                            numerator_learners = NULL,
                            bounds = c(0.01, 1),
                            truncation = NULL, truncation_quantile = NULL,
                            recompute = FALSE, verbose = TRUE) {
  obj <- .as_longy_data(obj)
  stabilization <- match.arg(stabilization)
  regime <- .resolve_regimes(obj, regime)

  if (!recompute) {
    computed <- Filter(function(r) !is.null(obj$weights[[r]]) &&
                         length(obj$weights[[r]]) > 0, regime)
    if (length(computed) > 0)
      stop(sprintf("Weights already computed for: %s. Use recompute=TRUE to override.",
                   paste(computed, collapse = ", ")), call. = FALSE)
  }

  for (rname in regime) {

  if (is.null(obj$fits$treatment[[rname]]) || length(obj$fits$treatment[[rname]]) == 0) {
    stop(sprintf("Treatment model not fit for regime '%s'. Run fit_treatment() first.", rname),
         call. = FALSE)
  }

  nodes <- obj$nodes
  reg <- obj$regimes[[rname]]
  id_col <- nodes$id

  # --- Treatment component: extract g_a and marginal ---
  gA <- obj$fits$treatment[[rname]]$predictions
  if (reg$type == "static") {
    gA_work <- data.table::copy(gA)
    gA_work[, .d := .resolve_static_at_time(reg$value, .time)]
    gA_follow <- gA_work[.treatment == .d]
    gA_follow[, .g_a := ifelse(.d == 1L, .p_a, 1 - .p_a)]
    gA_follow[, .marg_g_a := ifelse(.d == 1L, .marg_a, 1 - .marg_a)]
    if (stabilized) {
      gA_follow[, .sw_a := .marg_g_a / .g_a]
    } else {
      gA_follow[, .sw_a := 1 / .g_a]
    }
    gA_follow[, .d := NULL]
  } else {
    gA_follow <- data.table::copy(gA)
    gA_follow[, .g_a := .p_a]
    gA_follow[, .marg_g_a := .marg_a]
    gA_follow[, .sw_a := if (stabilized) .marg_a / .p_a else 1 / .p_a]
  }

  w <- gA_follow[, c(id_col, ".time", ".g_a", ".marg_g_a", ".sw_a"), with = FALSE]

  # For static regimes: enforce cumulative consistency — subjects contribute
  # at time t only if they followed the regime at ALL times 0..t
  if (reg$type == "static") {
    w <- .enforce_regime_consistency(w, id_col, ".time", obj$meta$time_values)
  }

  # --- Censoring component: extract g_c and marginal ---
  cens_fits <- obj$fits$censoring[[rname]]
  if (length(cens_fits) > 0) {
    for (cvar in names(cens_fits)) {
      gC <- cens_fits[[cvar]]$predictions
      # Keep uncensored rows (censored == 0)
      gC_uncens <- gC[gC$.censored == 0L, ]
      gC_uncens[, .sw_c := if (stabilized) .marg_c / .p_c else 1 / .p_c]

      cw <- gC_uncens[, c(id_col, ".time", ".p_c", ".marg_c", ".sw_c"), with = FALSE]
      data.table::setnames(cw, c(".p_c", ".marg_c", ".sw_c"),
                           c(paste0(".g_c_", cvar), paste0(".marg_c_", cvar),
                             paste0(".sw_c_", cvar)))
      w <- merge(w, cw, by = c(id_col, ".time"), all = FALSE)
    }

    # Combined censoring (product of all sources)
    c_g_cols <- paste0(".g_c_", names(cens_fits))
    c_m_cols <- paste0(".marg_c_", names(cens_fits))
    c_sw_cols <- paste0(".sw_c_", names(cens_fits))
    w[, .g_c := Reduce(`*`, .SD), .SDcols = c_g_cols]
    w[, .marg_g_c := Reduce(`*`, .SD), .SDcols = c_m_cols]
    w[, .sw_c := Reduce(`*`, .SD), .SDcols = c_sw_cols]
  } else {
    w[, .g_c := 1]
    w[, .marg_g_c := 1]
    w[, .sw_c := 1]
  }

  # --- Observation component: extract g_r and marginal ---
  obs_fit <- obj$fits$observation[[rname]]
  if (!is.null(obs_fit)) {
    gR <- obs_fit$predictions
    # Keep observed rows (observed == 1)
    gR_obs <- gR[gR$.observed == 1L, ]
    gR_obs[, .sw_r := if (stabilized) .marg_r / .p_r else 1 / .p_r]

    rw <- gR_obs[, c(id_col, ".time", ".p_r", ".marg_r", ".sw_r"), with = FALSE]
    data.table::setnames(rw, c(".p_r", ".marg_r"), c(".g_r", ".marg_g_r"))
    w <- merge(w, rw, by = c(id_col, ".time"), all = FALSE)
  } else {
    w[, .g_r := 1]
    w[, .marg_g_r := 1]
    w[, .sw_r := 1]
  }

  # --- Cumulate AC component, then bound ---
  # Compute raw point-in-time AC product, cumulate, then bound the cumulative
  # value. This directly prevents extreme weights from compounding over time.
  # Observation weights (intermittent) are applied separately, not cumulated.
  w[, .g_ac := .g_a * .g_c]
  data.table::setkeyv(w, c(id_col, ".time"))
  w[, .g_cum_ac := cumprod(.g_ac), by = c(id_col)]

  if (!is.null(bounds)) {
    n_trunc <- sum(w$.g_cum_ac < bounds[1])
    if (n_trunc > 0) {
      pct <- 100 * n_trunc / nrow(w)
      warning(sprintf(
        "Positivity: %d subject-time obs (%.1f%%) had cumulative g truncated to bounds[1]=%.4f.",
        n_trunc, pct, bounds[1]), call. = FALSE)
    }
    w[, .g_cum_ac := pmax(pmin(.g_cum_ac, bounds[2]), bounds[1])]
  }

  # --- Stabilized cumulative AC weight ---
  # .sw_a and .sw_c remain from raw probabilities (for diagnostics).
  # .csw_ac is derived from the bounded cumulative product.
  if (stabilized) {
    if (stabilization == "baseline") {
      # Baseline-stabilized: numerator = P(A=d|baseline) * P(C=0|baseline)
      # Fit on-demand if not already fitted
      if (is.null(obj$fits$numerator$treatment[[rname]])) {
        if (verbose) .vmsg("Fitting baseline numerator models for regime '%s'...", rname)
        obj <- .fit_baseline_numerator(obj, rname, "treatment",
                                        learners = numerator_learners,
                                        verbose = verbose)
      }
      if (length(obj$fits$censoring[[rname]]) > 0 &&
          is.null(obj$fits$numerator$censoring[[rname]])) {
        obj <- .fit_baseline_numerator(obj, rname, "censoring",
                                        learners = numerator_learners,
                                        verbose = verbose)
      }

      # Merge treatment numerator predictions
      num_trt <- obj$fits$numerator$treatment[[rname]]$predictions
      w <- merge(w, num_trt[, c(id_col, ".time", ".p_num"), with = FALSE],
                 by = c(id_col, ".time"), all.x = TRUE)
      # Transform by regime direction
      if (reg$type == "static") {
        w[, .d_num := .resolve_static_at_time(reg$value, .time)]
        w[, .num_g_a := ifelse(.d_num == 1L, .p_num, 1 - .p_num)]
        w[, .d_num := NULL]
      } else {
        w[, .num_g_a := .p_num]
      }
      w[, .p_num := NULL]

      # Merge censoring numerator predictions (product of all causes)
      if (!is.null(obj$fits$numerator$censoring[[rname]])) {
        num_cens <- obj$fits$numerator$censoring[[rname]]
        w[, .num_g_c := 1]
        for (cvar in names(num_cens)) {
          nc <- num_cens[[cvar]]$predictions
          # Censoring numerator: P(C=0|baseline) = 1 - P(C=1|baseline)
          nc_col <- data.table::copy(nc[, c(id_col, ".time", ".p_num"), with = FALSE])
          data.table::setnames(nc_col, ".p_num", ".p_num_c")
          w <- merge(w, nc_col, by = c(id_col, ".time"), all.x = TRUE)
          w[is.na(.p_num_c), .p_num_c := 1]
          w[, .num_g_c := .num_g_c * (1 - .p_num_c)]
          w[, .p_num_c := NULL]
        }
      } else {
        w[, .num_g_c := 1]
      }

      # Fill NA numerator predictions with marginal (subjects not in numerator risk set)
      w[is.na(.num_g_a), .num_g_a := .marg_g_a]

      w[, .num_cum_ac := cumprod(.num_g_a * .num_g_c), by = c(id_col)]
      w[, .csw_ac := .num_cum_ac / .g_cum_ac]
    } else {
      # Marginal stabilization (default)
      w[, .marg_cum_ac := cumprod(.marg_g_a * .marg_g_c), by = c(id_col)]
      w[, .csw_ac := .marg_cum_ac / .g_cum_ac]
    }
  } else {
    w[, .csw_ac := 1 / .g_cum_ac]
  }

  # Final weight: cumulative A*C * point-in-time R
  w[, .final_weight := .csw_ac * .sw_r]

  # --- Sampling weights (external, e.g., survey/population weights) ---
  if (!is.null(nodes$sampling_weights)) {
    sw_col <- nodes$sampling_weights
    sw_dt <- unique(obj$data[, c(id_col, sw_col), with = FALSE])
    w <- merge(w, sw_dt, by = id_col, all.x = TRUE)
    w[, .final_weight := .final_weight * get(sw_col)]
    w[, (sw_col) := NULL]
  }

  # --- Check for non-finite weights ---
  if (any(!is.finite(w$.final_weight))) {
    n_bad <- sum(!is.finite(w$.final_weight))
    warning(sprintf(
      "%d non-finite weight(s) detected. Consider tighter bounds or truncation.",
      n_bad), call. = FALSE)
  }

  # --- Truncation ---
  if (!is.null(truncation) && !is.null(truncation_quantile)) {
    warning("Both 'truncation' and 'truncation_quantile' specified. Using 'truncation' only.",
            call. = FALSE)
    truncation_quantile <- NULL
  }
  if (!is.null(truncation)) {
    w[, .final_weight := pmin(.final_weight, truncation)]
  }
  if (!is.null(truncation_quantile)) {
    q_val <- stats::quantile(w$.final_weight, truncation_quantile)
    w[, .final_weight := pmin(.final_weight, q_val)]
  }

  # --- ESS and extreme weight diagnostics (single consolidated warning) ---
  finite_w <- w$.final_weight[is.finite(w$.final_weight)]
  if (length(finite_w) > 0) {
    n_obs <- length(finite_w)
    ess_total <- .ess(finite_w)
    ess_pct <- 100 * ess_total / n_obs
    max_w <- max(finite_w)
    mean_w <- mean(finite_w)

    # Per-time-point ESS check
    time_ess <- w[is.finite(.final_weight),
                  list(.ess = .ess(.final_weight), .n = .N), by = ".time"]
    time_ess[, .ess_pct := 100 * .ess / .n]
    low_ess_times <- time_ess[.ess_pct < 15, .time]

    # Build single warning from all diagnostics
    if (ess_pct < 25) {
      parts <- sprintf("ESS=%.0f (%.0f%% of n=%d)", ess_total, ess_pct, n_obs)
      if (mean_w > 0 && max_w > 20 * mean_w)
        parts <- c(parts, sprintf("max weight %.0fx mean", max_w / mean_w))
      if (length(low_ess_times) > 0)
        parts <- c(parts, sprintf("worst at time(s) %s",
                                  paste(utils::head(low_ess_times, 3), collapse = ", ")))
      warning(paste0(paste(parts, collapse = "; "),
                     ". Consider truncation_quantile."), call. = FALSE)
    }
  }

  obj$weights[[rname]] <- list(
    regime = rname,
    weights_dt = w,
    stabilized = stabilized,
    stabilization = stabilization,
    numerator_learners = numerator_learners,
    bounds = bounds,
    truncation = truncation,
    truncation_quantile = truncation_quantile
  )

  } # end for (rname in regime)

  obj
}
