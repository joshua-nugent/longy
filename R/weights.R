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
#'   .g_cum).
#' @noRd
.compute_cumulative_g <- function(obj, regime, g_bounds = c(0.01, 1)) {
  nodes <- obj$nodes
  id_col <- nodes$id
  reg <- obj$regimes[[regime]]


  # --- Treatment component: g_a ---
  gA <- obj$fits$treatment$predictions
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

  # --- Censoring component: g_c ---
  if (length(obj$fits$censoring) > 0) {
    g_c_combined <- NULL
    for (cvar in names(obj$fits$censoring)) {
      gC <- obj$fits$censoring[[cvar]]$predictions
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
    c_cols <- paste0(".g_c_", names(obj$fits$censoring))
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

  # Bound cumulative g
  g_dt[, .g_cum := pmax(.g_cum, g_bounds[1])]
  g_dt[, .g_cum := pmin(.g_cum, g_bounds[2])]

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
#'
#' @return Modified `longy_data` object with weights stored.
#' @export
compute_weights <- function(obj, regime, stabilized = TRUE,
                            truncation = NULL, truncation_quantile = NULL,
                            g_bounds = c(0.01, 1)) {
  stopifnot(inherits(obj, "longy_data"))

  if (is.null(obj$fits$treatment)) {
    stop("Treatment model not fit. Run fit_treatment() first.", call. = FALSE)
  }

  nodes <- obj$nodes
  reg <- obj$regimes[[regime]]

  # --- Treatment weights ---
  gA <- obj$fits$treatment$predictions
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

  # --- Censoring weights ---
  if (length(obj$fits$censoring) > 0) {
    for (cvar in names(obj$fits$censoring)) {
      gC <- obj$fits$censoring[[cvar]]$predictions
      # Keep uncensored rows (censored == 0)
      gC_uncens <- gC[gC$.censored == 0L, ]
      gC_uncens[, .sw_c := if (stabilized) .marg_c / .p_c else 1 / .p_c]

      cw <- gC_uncens[, c(id_col, ".time", ".sw_c"), with = FALSE]
      data.table::setnames(cw, ".sw_c", paste0(".sw_c_", cvar))
      w <- merge(w, cw, by = c(id_col, ".time"), all = FALSE)
    }

    # Combined censoring weight (product of all sources)
    c_sw_cols <- paste0(".sw_c_", names(obj$fits$censoring))
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
  if (!is.null(obj$fits$observation)) {
    gR <- obj$fits$observation$predictions
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

  obj$weights <- list(
    regime = regime,
    weights_dt = w,
    stabilized = stabilized,
    truncation = truncation,
    truncation_quantile = truncation_quantile,
    g_bounds = g_bounds
  )

  obj
}
