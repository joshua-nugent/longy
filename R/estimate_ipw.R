#' Estimate Causal Effects via Inverse Probability Weighting
#'
#' Computes the Hajek (self-normalized) IPW estimator at each requested time
#' point, with standard errors and confidence intervals.
#'
#' \strong{Competing risks:} IPW is not supported when a \code{competing}
#' column is specified. The Hajek estimator cannot correctly attribute
#' cause-specific outcomes because subjects absorbed by a competing event may
#' drop out of the regime-follower risk set, losing their known Y=0
#' contribution. Use \code{\link{estimate_gcomp}} or
#' \code{\link{estimate_tmle}} instead — both handle competing risks
#' explicitly through the outcome model.
#'
#' Weight truncation is controlled via \code{truncation} and
#' \code{truncation_quantile} in \code{\link{compute_weights}()}. For TMLE,
#' use \code{g_bounds} in \code{\link{estimate_tmle}()} instead — it bounds
#' the unstabilized cumulative propensity score used in the clever covariate.
#'
#' @param obj A `longy_data` object with treatment (and optionally
#'   censoring/observation) models already fit. Weights are computed
#'   automatically if not already present.
#' @param regime Character. Name of the regime.
#' @param times Numeric vector. Time points at which to estimate. If NULL,
#'   estimates at all time points with data.
#' @param inference Character. Inference method: `"ic"` (influence curve),
#'   `"bootstrap"`, `"sandwich"`, or `"none"` (point estimates only).
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param n_boot Integer. Number of bootstrap replicates (only for `"bootstrap"`).
#' @param cluster Character. Column name for clustered standard errors. If
#'   NULL (default), uses \code{nodes$cluster} from the \code{longy_data}
#'   object when available.
#' @param stabilized Logical. Use stabilized weights? Only applies when
#'   weights are auto-computed (ignored if weights already exist). Default TRUE.
#' @param stabilization Character. Controls the numerator of stabilized IPW
#'   weights. \code{"marginal"} (default) uses unconditional marginal rates.
#'   \code{"baseline"} fits models using only baseline covariates. Only applies
#'   when weights are auto-computed.
#' @param numerator_learners Character vector. SuperLearner library for
#'   baseline numerator models. Default \code{NULL} uses \code{"SL.glm"}.
#'   Only applies when \code{stabilization = "baseline"} and weights are
#'   auto-computed.
#' @param bounds Numeric(2) or NULL. Bounds for the cumulative AC product
#'   (cumprod of g_a * g_c), passed to \code{\link{compute_weights}()}.
#'   Default \code{c(0.01, 1)}. Only applies when weights are auto-computed.
#' @param truncation Numeric or NULL. Hard upper bound for weight truncation,
#'   passed to \code{\link{compute_weights}()}. Only applies when weights are
#'   auto-computed.
#' @param truncation_quantile Numeric in (0,1) or NULL. Quantile upper bound for weight
#'   truncation, passed to \code{\link{compute_weights}()}. Only applies when
#'   weights are auto-computed.
#'
#' @return Modified \code{longy_data} object with IPW results stored in
#'   \code{obj$results}.
#'
#' @export
estimate_ipw <- function(obj, regime = NULL, times = NULL, inference = "ic",
                         ci_level = 0.95, n_boot = 200L, cluster = NULL,
                         stabilized = TRUE,
                         stabilization = c("marginal", "baseline"),
                         numerator_learners = NULL,
                         bounds = c(0.01, 1),
                         truncation = NULL,
                         truncation_quantile = NULL) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)
  inference <- match.arg(inference, c("ic", "bootstrap", "sandwich", "none"))

  if (ci_level <= 0 || ci_level >= 1)
    stop("ci_level must be between 0 and 1.", call. = FALSE)

  # IPW is not valid with competing risks
  if (!is.null(obj$nodes$competing))
    stop("IPW is not supported with competing risks. ",
         "Subjects absorbed by a competing event may drop out of the ",
         "regime-follower risk set, losing their known Y=0 contribution. ",
         "Use estimate_gcomp() or estimate_tmle() instead.",
         call. = FALSE)

  # Default cluster from nodes if not explicitly provided
  if (is.null(cluster) && !is.null(obj$nodes$cluster)) {
    cluster <- obj$nodes$cluster
  }

  # Validate cluster column exists in data
  if (!is.null(cluster) && !cluster %in% names(obj$data))
    stop(sprintf("Cluster column '%s' not found in data.", cluster),
         call. = FALSE)

  # Preserve user's original times so each regime gets the same input
  user_times <- times

  for (rname in regime) {

  # Reset times for each regime to avoid cross-contamination
  times <- user_times

  # Auto-compute weights if not yet computed
  if (is.null(obj$weights[[rname]]) || length(obj$weights[[rname]]) == 0) {
    if (is.null(obj$fits$treatment[[rname]])) {
      stop(sprintf(
        "No treatment model fitted for regime '%s'. Run fit_treatment() first.",
        rname), call. = FALSE)
    }
    obj <- compute_weights(obj, regime = rname, stabilized = stabilized,
                            stabilization = match.arg(stabilization),
                            numerator_learners = numerator_learners,
                            bounds = bounds,
                            truncation = truncation,
                            truncation_quantile = truncation_quantile,
                            recompute = TRUE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights[[rname]]$weights_dt

  # Determine time points
  available_times <- sort(unique(w_dt$.time))
  if (is.null(times)) {
    times <- available_times
  } else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s have no weight data, removing.",
                      paste(bad_times, collapse = ", ")))
      times <- intersect(times, available_times)
    }
    if (length(times) == 0)
      stop(sprintf("No valid time points for regime '%s'. Available: %s",
                   rname, paste(available_times, collapse = ", ")),
           call. = FALSE)
  }

  # Columns to merge from data (include cluster if specified)
  merge_cols <- unique(c(id_col, nodes$outcome))
  if (!is.null(cluster)) merge_cols <- unique(c(merge_cols, cluster))

  # Full baseline population size (IC denominator)
  N <- obj$meta$n_subjects
  all_ids <- unique(obj$data[[id_col]])

  # Point estimates (also collect merged data for IC inference)
  est_list <- vector("list", length(times))
  merged_list <- vector("list", length(times))
  ic_list <- vector("list", length(times))
  for (k in seq_along(times)) {
    tt <- times[k]
    w_t <- w_dt[list(tt), on = ".time", nomatch = NULL]
    dt_t <- obj$data[list(tt), on = nodes$time, nomatch = NULL]
    merged <- merge(w_t, dt_t[, merge_cols, with = FALSE], by = id_col)

    # Handle NA outcomes
    na_idx <- is.na(merged[[nodes$outcome]])
    if (any(na_idx)) {
      if (!is.null(obj$fits$observation[[rname]])) {
        # Observation model fitted → weight table only has R=1 subjects.
        # Y=NA among observed subjects is a data error.
        stop(sprintf(
          "Time %s: %d subjects have NA outcomes despite observation model indicating R=1. ",
          tt, sum(na_idx)),
          "This suggests a data integrity issue — observed subjects should have non-NA outcomes.",
          call. = FALSE)
      }
      # No observation model → unobserved subjects may remain; filter (complete case)
      merged <- merged[!na_idx, ]
    }
    merged_list[[k]] <- merged

    yi <- merged[[nodes$outcome]]
    wi <- merged$.final_weight
    n_at_risk <- nrow(merged)

    # Guard against empty data or all-zero weights (produces NaN)
    psi_hat <- if (n_at_risk > 0 && sum(wi) > 0) {
      stats::weighted.mean(yi, wi)
    } else {
      NA_real_
    }

    n_eff <- if (n_at_risk > 0 && sum(wi) > 0) .ess(wi) else 0

    # Warn if effective sample size is very low at this time point
    if (n_eff > 0 && n_at_risk >= 5 && n_eff / n_at_risk < 0.05) {
      warning(sprintf(
        "Time %s: ESS=%.1f is <5%% of n_at_risk=%d. Estimate dominated by few subjects.",
        tt, n_eff, n_at_risk), call. = FALSE)
    }

    est_list[[k]] <- data.table::data.table(
      time = tt,
      estimate = psi_hat,
      n_effective = n_eff,
      n_at_risk = n_at_risk
    )

    # Compute per-subject ICs for contrast inference
    # IC_i = N * w_i * (Y_i - psi) / sum(w) for at-risk; 0 otherwise
    if (n_at_risk >= 2 && sum(wi) > 0 && !is.na(psi_hat)) {
      IC_at_risk <- N * wi * (yi - psi_hat) / sum(wi)
      # Build full-population IC vector (0 for not-at-risk)
      ic_vec <- rep(0, N)
      names(ic_vec) <- as.character(all_ids)
      ic_vec[as.character(merged[[id_col]])] <- IC_at_risk
      ic_list[[k]] <- data.table::data.table(
        .id = all_ids,
        .time = tt,
        .ic = ic_vec
      )
      data.table::setnames(ic_list[[k]], ".id", id_col)
    }
  }
  estimates <- data.table::rbindlist(est_list)
  ic_dt <- if (length(Filter(Negate(is.null), ic_list)) > 0) {
    data.table::rbindlist(Filter(Negate(is.null), ic_list))
  } else {
    NULL
  }

  # Inference (computed on raw/unsmoothed estimates)
  if (inference == "ic") {
    inf_dt <- .ic_inference(estimates, obj, regime = rname, ci_level = ci_level,
                            cluster = cluster, merged_list = merged_list)
    estimates <- cbind(estimates, inf_dt)
  } else if (inference == "bootstrap") {
    inf_dt <- .bootstrap_inference(obj, regime = rname, times = times,
                                   n_boot = n_boot, ci_level = ci_level)
    estimates <- cbind(estimates, inf_dt)
  } else if (inference == "sandwich") {
    inf_dt <- .sandwich_inference(obj, regime = rname, times = times,
                                  ci_level = ci_level, cluster = cluster)
    estimates <- cbind(estimates, inf_dt)
  }

  # Isotonic smoothing for survival outcomes (enforce monotone non-decreasing)
  # Applied after inference so SEs are computed from unsmoothed ICs.
  # CIs are shifted to match smoothed estimates (preserves CI width).
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    if ("ci_lower" %in% names(estimates)) {
      shift_iso <- estimates$estimate - raw_est
      estimates$ci_lower <- estimates$ci_lower + shift_iso
      estimates$ci_upper <- estimates$ci_upper + shift_iso
    }
  }

  # Clamp CIs to [0, 1] for binary/survival outcomes
  if (nodes$outcome_type %in% c("binary", "survival") &&
      "ci_lower" %in% names(estimates)) {
    estimates$ci_lower <- pmax(estimates$ci_lower, 0)
    estimates$ci_upper <- pmin(estimates$ci_upper, 1)
  }

  result <- list(
    estimates = estimates,
    regime = rname,
    estimator = "ipw",
    inference = inference,
    ci_level = ci_level,
    ic = ic_dt
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_ipw")]] <- result

  } # end for (rname in regime)

  obj
}
