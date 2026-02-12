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
#'
#' @return Modified `longy_data` object with weights stored.
#' @export
compute_weights <- function(obj, regime, stabilized = TRUE,
                            truncation = NULL, truncation_quantile = NULL) {
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

  # --- Truncation ---
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
    truncation_quantile = truncation_quantile
  )

  obj
}
