#' Compute Unadjusted (Crude) Estimates
#'
#' Computes naive means of the outcome among uncensored, observed subjects
#' who naturally followed the specified regime through each time point.
#' No models, weights, or causal adjustment are applied. This serves as a
#' reference to show how much causal adjustment matters.
#'
#' Because no nuisance models are needed, this estimator can be called
#' directly after \code{\link{define_regime}} without any \code{fit_*} steps.
#'
#' @param obj A \code{longy_data} object with at least one regime defined.
#' @param regime Character. Name(s) of the regime(s). NULL = all defined.
#' @param times Numeric vector. Time points at which to estimate. NULL = all.
#' @param ci_level Numeric. Confidence level (default 0.95).
#'
#' @return Modified \code{longy_data} object with unadjusted results stored in
#'   \code{obj$results} (keyed as \code{{regime}_unadjusted}).
#'
#' @export
estimate_unadjusted <- function(obj, regime = NULL, times = NULL,
                                ci_level = 0.95) {
  obj <- .as_longy_data(obj)
  regime <- .resolve_regimes(obj, regime)

  if (ci_level <= 0 || ci_level >= 1)
    stop("ci_level must be between 0 and 1.", call. = FALSE)

  for (rname in regime) {

  nodes <- obj$nodes
  id_col <- nodes$id
  time_col <- nodes$time
  y_col <- nodes$outcome
  r_col <- nodes$observation

  # Add tracking columns (regime consistency + censoring status)
  .add_tracking_columns(obj$data, nodes, obj$regimes[[rname]])

  # Determine time points
  available_times <- sort(unique(obj$data[[time_col]]))
  eval_times <- if (is.null(times)) available_times else {
    bad_times <- setdiff(times, available_times)
    if (length(bad_times) > 0) {
      warning(sprintf("Time(s) %s not in data, removing.",
                      paste(bad_times, collapse = ", ")))
    }
    et <- intersect(times, available_times)
    if (length(et) == 0)
      stop(sprintf("No valid time points. Available: %s",
                   paste(available_times, collapse = ", ")),
           call. = FALSE)
    et
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  has_sw <- !is.null(nodes$sampling_weights)

  if (has_sw) {
    warning("Unadjusted estimates use sampling_weights as survey weights for the crude mean. ",
            "No causal adjustment is applied.", call. = FALSE)
  }

  est_list <- vector("list", length(eval_times))
  for (k in seq_along(eval_times)) {
    tt <- eval_times[k]
    dt_t <- obj$data[obj$data[[time_col]] == tt, ]

    # Subset: naturally followed regime AND uncensored through t
    keep <- dt_t$.longy_cum_consist == 1L & dt_t$.longy_cum_uncens == 1L
    keep[is.na(keep)] <- FALSE

    # Also require observed at t if observation node exists
    if (!is.null(r_col)) {
      obs_ok <- dt_t[[r_col]] == 1L
      obs_ok[is.na(obs_ok)] <- FALSE
      keep <- keep & obs_ok
    }

    # Require non-NA outcome
    yi_full <- dt_t[[y_col]][keep]
    not_na <- !is.na(yi_full)
    yi <- yi_full[not_na]
    n <- length(yi)

    # Extract aligned sampling weights
    sw_i <- NULL
    if (has_sw) {
      sw_i <- dt_t[[nodes$sampling_weights]][keep][not_na]
    }

    # Extract aligned cluster IDs for cluster-robust SEs
    cl_i <- NULL
    has_cl <- !is.null(nodes$cluster)
    if (has_cl) {
      cl_i <- dt_t[[nodes$cluster]][keep][not_na]
    }

    if (n > 0) {
      psi_hat <- if (!is.null(sw_i)) stats::weighted.mean(yi, sw_i) else mean(yi)
      # SE: use Kish's effective sample size when sampling weights are present
      if (!is.null(sw_i)) {
        n_eff <- sum(sw_i)^2 / sum(sw_i^2)
        if (nodes$outcome_type %in% c("binary", "survival")) {
          se <- sqrt(psi_hat * (1 - psi_hat) / n_eff)
        } else {
          w_var <- sum(sw_i * (yi - psi_hat)^2) / sum(sw_i)
          se <- sqrt(w_var / n_eff)
        }
      } else if (has_cl) {
        # Cluster-robust SE: compute cluster-level means, then variance of
        # cluster means around psi_hat
        n_eff <- n
        cl_means <- tapply(yi, cl_i, mean)
        all_clusters <- unique(obj$data[[nodes$cluster]])
        N_cl <- length(all_clusters)
        # Pad with psi_hat for clusters not in this risk set (contribute 0 IC)
        missing_cl <- setdiff(as.character(all_clusters), names(cl_means))
        cl_means_full <- c(as.numeric(cl_means), rep(psi_hat, length(missing_cl)))
        se <- sqrt(stats::var(cl_means_full) / N_cl)
      } else {
        n_eff <- n
        if (nodes$outcome_type %in% c("binary", "survival")) {
          se <- sqrt(psi_hat * (1 - psi_hat) / n)
        } else {
          se <- stats::sd(yi) / sqrt(n)
        }
      }
    } else {
      psi_hat <- NA_real_
      se <- NA_real_
      n_eff <- 0
    }

    est_list[[k]] <- data.table::data.table(
      time = tt,
      estimate = psi_hat,
      n_effective = as.numeric(n_eff),
      n_at_risk = n,
      se = se,
      ci_lower = psi_hat - z * se,
      ci_upper = psi_hat + z * se
    )
  }
  estimates <- data.table::rbindlist(est_list)

  # Isotonic smoothing for survival outcomes (enforce monotone non-decreasing)
  if (nodes$outcome_type == "survival" && nrow(estimates) > 1) {
    raw_est <- estimates$estimate
    iso <- stats::isoreg(raw_est)
    estimates$estimate <- iso$yf
    estimates$ci_lower <- estimates$estimate - z * estimates$se
    estimates$ci_upper <- estimates$estimate + z * estimates$se
    # Clamp CIs to [0, 1] for survival outcomes
    estimates$ci_lower <- pmax(estimates$ci_lower, 0)
    estimates$ci_upper <- pmin(estimates$ci_upper, 1)
  }

  # Clean up tracking columns
  .remove_tracking_columns(obj$data)

  result <- list(
    estimates = estimates,
    regime = rname,
    estimator = "unadjusted",
    inference = "analytic",
    ci_level = ci_level
  )
  class(result) <- "longy_result"
  obj$results[[paste0(rname, "_unadjusted")]] <- result

  } # end for (rname in regime)

  obj
}
