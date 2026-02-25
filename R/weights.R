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
    if (reg$value == 1L) {
      # P(A=1|past) for regime-followers who had A=1
      gA_sub <- gA[gA$.treatment == 1L, ]
      g_dt <- gA_sub[, c(id_col, ".time", ".p_a"), with = FALSE]
      data.table::setnames(g_dt, ".p_a", ".g_a")
    } else {
      gA_sub <- gA[gA$.treatment == 0L, ]
      g_dt <- gA_sub[, c(id_col, ".time", ".p_a"), with = FALSE]
      # g_a = P(A=0|past) = 1 - p_a
      g_dt[, .p_a := 1 - .p_a]
      data.table::setnames(g_dt, ".p_a", ".g_a")
    }
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

#' Compute Inverse Probability Weights
#'
#' Combines treatment, censoring, and observation model predictions into
#' stabilized inverse probability weights. Treatment and censoring weights
#' are cumulated (absorbing processes); observation weights are point-in-time
#' (intermittent process).
#'
#' @param obj A `longy_data` object with treatment (and optionally censoring/observation)
#'   models already fit.
#' @param regime Character. Name of the regime.
#' @param stabilized Logical. Use stabilized weights (numerator = marginal probability).
#' @param truncation Numeric. Hard upper bound for weights. NULL for no truncation.
#' @param truncation_quantile Numeric in (0,1). Truncate at this quantile.
#'   Applied after `truncation` if both specified.
#' @param g_bounds Numeric vector of length 2. Bounds for cumulative g
#'   (denominator). Default \code{c(0.01, 1)}. Stored for use by TMLE.
#' @param recompute Logical. If FALSE (default), errors when weights are already
#'   computed for the requested regime(s). Set to TRUE to re-compute.
#'
#' @return Modified `longy_data` object with weights stored.
#' @export
compute_weights <- function(obj, regime = NULL, stabilized = TRUE,
                            truncation = NULL, truncation_quantile = NULL,
                            g_bounds = c(0.01, 1), recompute = FALSE) {
  obj <- .as_longy_data(obj)
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

  # --- Treatment weights ---
  gA <- obj$fits$treatment[[rname]]$predictions
  if (reg$type == "static") {
    if (reg$value == 1L) {
      # Keep only rows where treatment == 1 (regime followers)
      gA_follow <- gA[gA$.treatment == 1L, ]
      gA_follow[, .sw_a := if (stabilized) .marg_a / .p_a else 1 / .p_a]
    } else {
      gA_follow <- gA[gA$.treatment == 0L, ]
      gA_follow[, .sw_a := if (stabilized) (1 - .marg_a) / (1 - .p_a)
                            else 1 / (1 - .p_a)]
    }
  } else {
    # Dynamic/stochastic: all rows are "followers" (consistency built into risk set)
    gA_follow <- data.table::copy(gA)
    gA_follow[, .sw_a := if (stabilized) .marg_a / .p_a else 1 / .p_a]
  }

  id_col <- nodes$id
  w <- gA_follow[, c(id_col, ".time", ".sw_a"), with = FALSE]

  # For static regimes: enforce cumulative consistency — subjects contribute
  # at time t only if they followed the regime at ALL times 0..t
  if (reg$type == "static") {
    w <- .enforce_regime_consistency(w, id_col, ".time", obj$meta$time_values)
  }

  # --- Censoring weights ---
  cens_fits <- obj$fits$censoring[[rname]]
  if (length(cens_fits) > 0) {
    for (cvar in names(cens_fits)) {
      gC <- cens_fits[[cvar]]$predictions
      # Keep uncensored rows (censored == 0)
      gC_uncens <- gC[gC$.censored == 0L, ]
      gC_uncens[, .sw_c := if (stabilized) .marg_c / .p_c else 1 / .p_c]

      cw <- gC_uncens[, c(id_col, ".time", ".sw_c"), with = FALSE]
      data.table::setnames(cw, ".sw_c", paste0(".sw_c_", cvar))
      w <- merge(w, cw, by = c(id_col, ".time"), all = FALSE)
    }

    # Combined censoring weight (product of all sources)
    c_sw_cols <- paste0(".sw_c_", names(cens_fits))
    w[, .sw_c := Reduce(`*`, .SD), .SDcols = c_sw_cols]
  } else {
    w[, .sw_c := 1]
  }

  # Combined A*C point-in-time weight
  w[, .sw_ac := .sw_a * .sw_c]

  # Cumulative A*C weight (absorbing processes)
  data.table::setkeyv(w, c(id_col, ".time"))
  w[, .csw_ac := cumprod(.sw_ac), by = c(id_col)]

  # --- Observation weights (point-in-time, NOT cumulated) ---
  obs_fit <- obj$fits$observation[[rname]]
  if (!is.null(obs_fit)) {
    gR <- obs_fit$predictions
    # Keep observed rows (observed == 1)
    gR_obs <- gR[gR$.observed == 1L, ]
    gR_obs[, .sw_r := if (stabilized) .marg_r / .p_r else 1 / .p_r]

    rw <- gR_obs[, c(id_col, ".time", ".sw_r"), with = FALSE]
    w <- merge(w, rw, by = c(id_col, ".time"), all = FALSE)
  } else {
    w[, .sw_r := 1]
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

  # --- ESS and extreme weight diagnostics ---
  finite_w <- w$.final_weight[is.finite(w$.final_weight)]
  if (length(finite_w) > 0) {
    n_obs <- length(finite_w)
    ess_total <- .ess(finite_w)
    ess_pct <- 100 * ess_total / n_obs
    if (ess_pct < 10) {
      warning(sprintf(
        "ESS=%.0f is <10%% of n=%d. Weights are highly variable. Consider truncation or tighter g_bounds.",
        ess_total, n_obs), call. = FALSE)
    } else if (ess_pct < 25) {
      warning(sprintf(
        "ESS=%.0f is %.0f%% of n=%d. Weights are moderately variable.",
        ess_total, ess_pct, n_obs), call. = FALSE)
    }

    max_w <- max(finite_w)
    mean_w <- mean(finite_w)
    if (mean_w > 0 && max_w > 20 * mean_w) {
      warning(sprintf(
        "Maximum weight (%.1f) is %.0fx the mean. Consider truncation_quantile to cap.",
        max_w, max_w / mean_w), call. = FALSE)
    }

    # Per-time-point ESS check
    time_ess <- w[is.finite(.final_weight),
                  list(.ess = .ess(.final_weight), .n = .N), by = ".time"]
    time_ess[, .ess_pct := 100 * .ess / .n]
    low_ess_times <- time_ess[.ess_pct < 15, .time]
    if (length(low_ess_times) > 0) {
      warning(sprintf(
        "Low ESS (<15%% of n) at time point(s): %s.",
        paste(low_ess_times, collapse = ", ")), call. = FALSE)
    }
  }

  obj$weights[[rname]] <- list(
    regime = rname,
    weights_dt = w,
    stabilized = stabilized,
    truncation = truncation,
    truncation_quantile = truncation_quantile,
    g_bounds = g_bounds
  )

  } # end for (rname in regime)

  obj
}
