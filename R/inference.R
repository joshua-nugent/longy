#' Influence-Curve Based Inference for IPW Estimator
#'
#' Computes standard errors and confidence intervals using the influence curve.
#'
#' @param estimates data.table with columns: time, estimate, Y_weighted, weights, ids
#' @param obj longy_data object
#' @param ci_level Confidence level
#' @param cluster Character. Column name for cluster-robust SEs.
#' @return data.table with se, ci_lower, ci_upper added
#' @noRd
.ic_inference <- function(estimates, obj, ci_level = 0.95, cluster = NULL) {
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights$weights_dt

  results <- vector("list", nrow(estimates))

  for (j in seq_len(nrow(estimates))) {
    tt <- estimates$time[j]
    psi_hat <- estimates$estimate[j]

    # Get individual data for this time point
    w_t <- w_dt[w_dt$.time == tt, ]
    # Merge with outcome
    dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
    merged <- merge(w_t, dt_t[, c(id_col, nodes$outcome), with = FALSE],
                    by = id_col)

    wi <- merged$.final_weight
    yi <- merged[[nodes$outcome]]
    n_i <- nrow(merged)
    sum_w <- sum(wi)

    # Person-level IC (small scale: psi_hat - psi ≈ sum(ic))
    ic <- wi * (yi - psi_hat) / sum_w

    # Scale to ltmle convention (large scale: psi_hat - psi ≈ mean(IC))
    # so that var(IC)/n gives the right variance. Uses centered variance
    # (robust in finite samples). Follows ltmle's HouseholdIC approach.
    IC <- ic * n_i

    if (!is.null(cluster) && cluster %in% names(dt_t)) {
      # Cluster-robust: sum IC within clusters, rescale by n_cl/n_i,
      # then var/n_cl (ltmle's HouseholdIC pattern)
      merged_cl <- merge(merged, dt_t[, c(id_col, cluster), with = FALSE],
                         by = id_col)
      cl_IC <- tapply(IC, merged_cl[[cluster]], sum)
      n_cl <- length(cl_IC)
      cl_IC <- as.numeric(cl_IC) * n_cl / n_i
      se <- sqrt(stats::var(cl_IC) / n_cl)
    } else {
      se <- sqrt(stats::var(IC) / n_i)
    }

    results[[j]] <- data.table::data.table(
      se = se,
      ci_lower = psi_hat - z * se,
      ci_upper = psi_hat + z * se
    )
  }

  data.table::rbindlist(results)
}

#' Bootstrap Inference for IPW Estimator
#'
#' Resamples subjects (with replacement), re-computes weights and estimates.
#'
#' @param obj longy_data object (fully fitted)
#' @param regime Character regime name
#' @param times Numeric vector of time points
#' @param n_boot Integer number of bootstrap samples
#' @param ci_level Confidence level
#' @return data.table with se, ci_lower, ci_upper
#' @noRd
.bootstrap_inference <- function(obj, regime, times, n_boot = 200L,
                                 ci_level = 0.95) {
  nodes <- obj$nodes
  id_col <- nodes$id
  ids <- unique(obj$data[[id_col]])
  n <- length(ids)

  boot_estimates <- matrix(NA_real_, nrow = n_boot, ncol = length(times))

  for (b in seq_len(n_boot)) {
    boot_ids <- sample(ids, n, replace = TRUE)
    # Create bootstrap dataset
    boot_id_dt <- data.table::data.table(.boot_orig_id = boot_ids,
                                          .boot_new_id = seq_len(n))

    boot_rows <- lapply(seq_len(n), function(i) {
      rows <- obj$data[obj$data[[id_col]] == boot_ids[i], ]
      rows <- data.table::copy(rows)
      rows[[id_col]] <- i
      rows
    })
    boot_dt <- data.table::rbindlist(boot_rows)

    # Re-run pipeline on bootstrap data
    boot_obj <- tryCatch({
      b_obj <- longy_data(
        data = boot_dt,
        id = id_col, time = nodes$time,
        outcome = nodes$outcome, treatment = nodes$treatment,
        censoring = nodes$censoring, observation = nodes$observation,
        baseline = nodes$baseline, timevarying = nodes$timevarying,
        sampling_weights = nodes$sampling_weights,
        outcome_type = nodes$outcome_type, verbose = FALSE
      )
      b_obj$regimes <- obj$regimes
      b_obj <- fit_treatment(b_obj, regime = regime,
                             covariates = obj$fits$treatment$covariates,
                             learners = obj$fits$treatment$learners,
                             bounds = obj$fits$treatment$bounds,
                             verbose = FALSE)
      if (length(obj$fits$censoring) > 0) {
        cov_c <- obj$fits$censoring[[1]]$covariates
        lrn_c <- obj$fits$censoring[[1]]$learners
        bnd_c <- obj$fits$censoring[[1]]$bounds
        b_obj <- fit_censoring(b_obj, regime = regime,
                               covariates = cov_c, learners = lrn_c,
                               bounds = bnd_c, verbose = FALSE)
      }
      if (!is.null(obj$fits$observation)) {
        b_obj <- fit_observation(b_obj, regime = regime,
                                 covariates = obj$fits$observation$covariates,
                                 learners = obj$fits$observation$learners,
                                 bounds = obj$fits$observation$bounds,
                                 verbose = FALSE)
      }
      b_obj <- compute_weights(b_obj, regime = regime,
                               stabilized = obj$weights$stabilized,
                               truncation = obj$weights$truncation,
                               truncation_quantile = obj$weights$truncation_quantile)
      b_obj
    }, error = function(e) NULL)

    if (is.null(boot_obj)) next

    # Compute point estimates
    for (k in seq_along(times)) {
      tt <- times[k]
      est <- .hajek_estimate(boot_obj, tt)
      boot_estimates[b, k] <- est
    }
  }

  # Percentile CIs
  alpha <- 1 - ci_level
  results <- vector("list", length(times))
  for (k in seq_along(times)) {
    boot_k <- boot_estimates[, k]
    boot_k <- boot_k[!is.na(boot_k)]
    if (length(boot_k) < 3) {
      results[[k]] <- data.table::data.table(se = NA_real_,
                                              ci_lower = NA_real_,
                                              ci_upper = NA_real_)
    } else {
      results[[k]] <- data.table::data.table(
        se = stats::sd(boot_k),
        ci_lower = stats::quantile(boot_k, alpha / 2),
        ci_upper = stats::quantile(boot_k, 1 - alpha / 2)
      )
    }
  }

  data.table::rbindlist(results)
}

#' Sandwich Inference Using survey Package
#'
#' @param obj longy_data object
#' @param times Time points
#' @param ci_level Confidence level
#' @param cluster Cluster column
#' @return data.table with se, ci_lower, ci_upper
#' @noRd
.sandwich_inference <- function(obj, times, ci_level = 0.95, cluster = NULL) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for sandwich inference.", call. = FALSE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights$weights_dt

  results <- vector("list", length(times))

  for (k in seq_along(times)) {
    tt <- times[k]
    w_t <- w_dt[w_dt$.time == tt, ]
    dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
    merged <- merge(w_t, dt_t, by = id_col)

    merged_df <- as.data.frame(merged)

    if (!is.null(cluster) && cluster %in% names(merged_df)) {
      des <- survey::svydesign(ids = stats::reformulate(cluster),
                               weights = ~.final_weight,
                               data = merged_df)
    } else {
      des <- survey::svydesign(ids = ~1,
                               weights = ~.final_weight,
                               data = merged_df)
    }

    fml <- stats::reformulate(nodes$outcome)
    est <- survey::svymean(fml, des)

    ci <- stats::confint(est, level = ci_level)

    results[[k]] <- data.table::data.table(
      se = survey::SE(est)[1],
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  }

  data.table::rbindlist(results)
}

#' Hajek point estimate at a single time
#' @noRd
.hajek_estimate <- function(obj, tt) {
  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights$weights_dt
  w_t <- w_dt[w_dt$.time == tt, ]
  dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
  merged <- merge(w_t, dt_t[, c(id_col, nodes$outcome), with = FALSE],
                  by = id_col)
  yi <- merged[[nodes$outcome]]
  wi <- merged$.final_weight
  if (length(wi) == 0) return(NA_real_)
  stats::weighted.mean(yi, wi)
}
