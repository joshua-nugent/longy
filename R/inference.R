#' Dispatch bootstrap replicates via future.apply (if available) or lapply
#'
#' When \code{future.apply} is installed and a non-sequential
#' \code{future::plan()} is active, replicates run in parallel via
#' \code{future_lapply}. Otherwise falls back to sequential \code{lapply}.
#'
#' @param n_boot Integer number of replicates
#' @param one_boot_fn Function(b) returning a numeric vector of estimates
#' @param verbose Logical. Print a progress message.
#' @return List of length \code{n_boot} (each element a numeric vector or NULL)
#' @noRd
.run_bootstrap <- function(n_boot, one_boot_fn, verbose) {
  use_future <- requireNamespace("future.apply", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) &&
    !inherits(future::plan(), "sequential")
  if (verbose) {
    .vmsg("  Running %d bootstrap replicates (%s)...",
          n_boot, if (use_future) "parallel" else "sequential")
  }
  if (use_future) {
    # Bootstrap closures capture the full longy_data object + fitted models.
    # Temporarily raise future's global size limit so serialization doesn't fail.
    old_maxSize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = +Inf)
    on.exit(options(future.globals.maxSize = old_maxSize), add = TRUE)
    future.apply::future_lapply(seq_len(n_boot), one_boot_fn,
                                future.seed = TRUE)
  } else {
    lapply(seq_len(n_boot), one_boot_fn)
  }
}

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
.ic_inference <- function(estimates, obj, regime, ci_level = 0.95,
                          cluster = NULL, merged_list = NULL) {
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights[[regime]]$weights_dt

  # Full baseline population size — the IC denominator.
  # Subjects not at risk at time t contribute IC=0. The weights already

  # account for censoring/dropout, so dividing by N (not n_t) avoids
  # double-counting attrition. This is the Hajek (self-normalized) IC:
  # IC_i = N * w_i * (Y_i - psi) / sum(w). See Hernán & Robins Ch 12.
  N <- obj$meta$n_subjects

  results <- vector("list", nrow(estimates))

  for (j in seq_len(nrow(estimates))) {
    tt <- estimates$time[j]
    psi_hat <- estimates$estimate[j]

    # Use pre-computed merged data if available, otherwise compute
    if (!is.null(merged_list) && j <= length(merged_list)) {
      merged <- merged_list[[j]]
    } else {
      w_t <- w_dt[w_dt$.time == tt, ]
      dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
      merge_cols <- unique(c(id_col, nodes$outcome))
      if (!is.null(cluster)) merge_cols <- unique(c(merge_cols, cluster))
      merged <- merge(w_t, dt_t[, merge_cols, with = FALSE], by = id_col)
      # Filter NA outcomes
      merged <- merged[!is.na(merged[[nodes$outcome]]), ]
    }

    wi <- merged$.final_weight
    yi <- merged[[nodes$outcome]]
    n_t <- nrow(merged)
    sum_w <- sum(wi)

    if (n_t < 2 || sum_w == 0) {
      results[[j]] <- data.table::data.table(
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
      next
    }

    # Person-level IC for at-risk subjects, scaled by N (full population).
    # IC_i = N * w_i * (Y_i - psi) / sum(w) for at-risk subjects;
    # IC_i = 0 for subjects not at risk (absorbed by weights).
    # Variance: var(IC_padded) / N where IC_padded includes N - n_t zeros.
    IC <- N * wi * (yi - psi_hat) / sum_w

    if (!is.null(cluster) && cluster %in% names(merged)) {
      # Cluster-robust: sum IC within clusters, pad with zero-IC clusters,
      # then var/N_cl. Follows ltmle's HouseholdIC pattern.
      cl_IC <- tapply(IC, merged[[cluster]], sum)
      # Total clusters in population (not just at-risk clusters)
      all_clusters <- unique(obj$data[[cluster]])
      N_cl <- length(all_clusters)
      # Pad with zeros for clusters not at risk at this time
      n_cl_at_risk <- length(cl_IC)
      cl_IC_full <- c(as.numeric(cl_IC), rep(0, N_cl - n_cl_at_risk))
      se <- sqrt(stats::var(cl_IC_full) / N_cl)
    } else {
      # Pad IC with zeros for subjects not at risk, compute var/N
      IC_full <- c(IC, rep(0, N - n_t))
      se <- sqrt(stats::var(IC_full) / N)
    }

    # Degenerate IC (all outcomes identical) gives SE=0;
    # report NA since zero variance is uninformative, not certain
    if (!is.na(se) && se == 0) {
      se <- NA_real_
    }

    results[[j]] <- data.table::data.table(
      se = se,
      ci_lower = if (is.na(se)) NA_real_ else psi_hat - z * se,
      ci_upper = if (is.na(se)) NA_real_ else psi_hat + z * se
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
                                 ci_level = 0.95, verbose = TRUE) {
  nodes <- obj$nodes
  id_col <- nodes$id
  ids <- unique(obj$data[[id_col]])
  n <- length(ids)

  # Use original censoring column name (not internal binary col names)
  cens_col <- nodes$censoring_col  # NULL if no censoring

  one_boot <- function(b) {
    boot_ids <- sample(ids, n, replace = TRUE)

    # Efficient bootstrap: single data.table join instead of n copies + rbind
    boot_map <- data.table::data.table(
      .boot_orig = boot_ids,
      .boot_new = seq_len(n)
    )
    data.table::setnames(boot_map, ".boot_orig", id_col)
    boot_dt <- merge(boot_map, obj$data, by = id_col, allow.cartesian = TRUE)
    boot_dt[[id_col]] <- boot_dt$.boot_new
    boot_dt[, .boot_new := NULL]

    boot_obj <- tryCatch({
      b_obj <- longy_data(
        data = boot_dt,
        id = id_col, time = nodes$time,
        outcome = nodes$outcome, treatment = nodes$treatment,
        censoring = cens_col, observation = nodes$observation,
        baseline = nodes$baseline, timevarying = nodes$timevarying,
        cluster = nodes$cluster, sampling_weights = nodes$sampling_weights,
        outcome_type = nodes$outcome_type,
        competing = nodes$competing, verbose = FALSE
      )
      b_obj$regimes <- obj$regimes
      # Force sequential SL inside bootstrap to prevent nested parallelism
      trt_fit <- obj$fits$treatment[[regime]]
      b_obj <- fit_treatment(b_obj, regime = regime,
                             covariates = trt_fit$covariates,
                             learners = trt_fit$learners,
                             sl_control = if (!is.null(trt_fit$sl_control)) trt_fit$sl_control else list(),
                             adaptive_cv = if (!is.null(trt_fit$adaptive_cv)) trt_fit$adaptive_cv else TRUE,
                             min_obs = if (!is.null(trt_fit$min_obs)) trt_fit$min_obs else 50L,
                             min_events = if (!is.null(trt_fit$min_events)) trt_fit$min_events else 20L,
                             risk_set = if (!is.null(trt_fit$risk_set)) trt_fit$risk_set else "all",
                             use_ffSL = FALSE,
                             verbose = FALSE)
      cens_fits <- obj$fits$censoring[[regime]]
      if (length(cens_fits) > 0) {
        cov_c <- cens_fits[[1]]$covariates
        lrn_c <- cens_fits[[1]]$learners
        slc_c <- if (!is.null(cens_fits[[1]]$sl_control)) cens_fits[[1]]$sl_control else list()
        acv_c <- if (!is.null(cens_fits[[1]]$adaptive_cv)) cens_fits[[1]]$adaptive_cv else TRUE
        mob_c <- if (!is.null(cens_fits[[1]]$min_obs)) cens_fits[[1]]$min_obs else 50L
        mev_c <- if (!is.null(cens_fits[[1]]$min_events)) cens_fits[[1]]$min_events else 20L
        b_obj <- fit_censoring(b_obj, regime = regime,
                               covariates = cov_c, learners = lrn_c,
                               sl_control = slc_c,
                               adaptive_cv = acv_c,
                               min_obs = mob_c, min_events = mev_c,
                               use_ffSL = FALSE,
                               verbose = FALSE)
      }
      obs_fit <- obj$fits$observation[[regime]]
      if (!is.null(obs_fit)) {
        b_obj <- fit_observation(b_obj, regime = regime,
                                 covariates = obs_fit$covariates,
                                 learners = obs_fit$learners,
                                 sl_control = if (!is.null(obs_fit$sl_control)) obs_fit$sl_control else list(),
                                 adaptive_cv = if (!is.null(obs_fit$adaptive_cv)) obs_fit$adaptive_cv else TRUE,
                                 min_obs = if (!is.null(obs_fit$min_obs)) obs_fit$min_obs else 50L,
                                 min_events = if (!is.null(obs_fit$min_events)) obs_fit$min_events else 20L,
                                 use_ffSL = FALSE,
                                 verbose = FALSE)
      }
      b_obj <- compute_weights(b_obj, regime = regime,
                               stabilized = obj$weights[[regime]]$stabilized,
                               stabilization = obj$weights[[regime]]$stabilization %||% "marginal",
                               numerator_learners = obj$weights[[regime]]$numerator_learners,
                               bounds = obj$weights[[regime]]$bounds,
                               truncation = obj$weights[[regime]]$truncation,
                               truncation_quantile = obj$weights[[regime]]$truncation_quantile,
                               verbose = FALSE)
      b_obj
    }, error = function(e) NULL)

    if (is.null(boot_obj)) return(rep(NA_real_, length(times)))

    vapply(times, function(tt) .hajek_estimate(boot_obj, regime, tt), numeric(1))
  }

  boot_list <- .run_bootstrap(n_boot, one_boot, verbose)
  boot_estimates <- do.call(rbind, boot_list)

  # Warn if many bootstrap replicates failed
  n_failed <- sum(vapply(boot_list, function(x) all(is.na(x)), logical(1)))
  if (n_failed == n_boot) {
    warning("All IPW bootstrap replicates failed. SEs and CIs will be NA.",
            call. = FALSE)
  } else if (n_failed > n_boot * 0.5) {
    warning(sprintf("IPW bootstrap: %d/%d replicates failed. SEs may be unreliable.",
                    n_failed, n_boot), call. = FALSE)
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
.sandwich_inference <- function(obj, regime, times, ci_level = 0.95, cluster = NULL) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for sandwich inference. Install with install.packages('survey').",
         call. = FALSE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights[[regime]]$weights_dt

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

#' Bootstrap Inference for G-Computation Estimator
#'
#' Resamples subjects (with replacement), re-fits outcome models, and
#' re-computes G-comp point estimates. Returns percentile CIs and bootstrap SEs.
#'
#' @param obj longy_data object (with outcome model fitted)
#' @param regime Character regime name
#' @param times Numeric vector of time points
#' @param n_boot Integer number of bootstrap samples
#' @param ci_level Confidence level
#' @param verbose Logical. Print progress.
#' @return data.table with se, ci_lower, ci_upper
#' @noRd
.bootstrap_gcomp_inference <- function(obj, regime, times, n_boot = 200L,
                                       ci_level = 0.95, verbose = TRUE) {
  nodes <- obj$nodes
  id_col <- nodes$id
  ids <- unique(obj$data[[id_col]])
  n <- length(ids)
  outcome_fit <- obj$fits$outcome[[regime]]

  # Use original censoring column name (not internal binary col names)
  cens_col <- nodes$censoring_col  # NULL if no censoring

  one_boot <- function(b) {
    boot_ids <- sample(ids, n, replace = TRUE)

    # Efficient bootstrap: single data.table join instead of n copies + rbind
    boot_map <- data.table::data.table(
      .boot_orig = boot_ids,
      .boot_new = seq_len(n)
    )
    data.table::setnames(boot_map, ".boot_orig", id_col)
    boot_dt <- merge(boot_map, obj$data, by = id_col, allow.cartesian = TRUE)
    boot_dt[[id_col]] <- boot_dt$.boot_new
    boot_dt[, .boot_new := NULL]

    boot_obj <- tryCatch({
      b_obj <- longy_data(
        data = boot_dt,
        id = id_col, time = nodes$time,
        outcome = nodes$outcome, treatment = nodes$treatment,
        censoring = cens_col, observation = nodes$observation,
        baseline = nodes$baseline, timevarying = nodes$timevarying,
        cluster = nodes$cluster, sampling_weights = nodes$sampling_weights,
        outcome_type = nodes$outcome_type,
        competing = nodes$competing, verbose = FALSE
      )
      b_obj$regimes <- obj$regimes

      # Force sequential SL inside bootstrap to prevent nested parallelism
      # (ffSL spawns futures inside future workers → resource explosion)
      b_obj <- fit_outcome(b_obj, regime = regime,
                           covariates = outcome_fit$covariates,
                           learners = outcome_fit$learners,
                           sl_control = if (!is.null(outcome_fit$sl_control)) outcome_fit$sl_control else list(),
                           adaptive_cv = if (!is.null(outcome_fit$adaptive_cv)) outcome_fit$adaptive_cv else TRUE,
                           min_obs = if (!is.null(outcome_fit$min_obs)) outcome_fit$min_obs else 50L,
                           bounds = outcome_fit$bounds,
                           risk_set = if (!is.null(outcome_fit$risk_set)) outcome_fit$risk_set else "all",
                           times = times,
                           use_ffSL = FALSE,
                           verbose = FALSE)
      b_obj
    }, error = function(e) NULL)

    if (is.null(boot_obj)) return(rep(NA_real_, length(times)))

    pred_dt <- boot_obj$fits$outcome[[regime]]$predictions
    if (is.null(pred_dt)) return(rep(NA_real_, length(times)))
    # Extract per-subject sampling weights from bootstrap data
    boot_sw_dt <- NULL
    boot_id_col <- nodes$id
    if (!is.null(nodes$sampling_weights)) {
      boot_sw_dt <- unique(boot_obj$data[, c(boot_id_col, nodes$sampling_weights), with = FALSE])
    }
    vapply(times, function(tt) {
      q_t <- pred_dt[pred_dt$.target_time == tt, ]
      if (nrow(q_t) > 0) {
        if (!is.null(boot_sw_dt)) {
          q_sw <- merge(q_t, boot_sw_dt, by = boot_id_col, all.x = TRUE)
          sw_vals <- q_sw[[nodes$sampling_weights]]
          sw_vals[is.na(sw_vals)] <- 1
          stats::weighted.mean(q_t$.Q_hat, sw_vals)
        } else {
          mean(q_t$.Q_hat)
        }
      } else NA_real_
    }, numeric(1))
  }

  boot_list <- .run_bootstrap(n_boot, one_boot, verbose)
  boot_estimates <- do.call(rbind, boot_list)

  # Warn if many bootstrap replicates failed
  n_failed <- sum(vapply(boot_list, function(x) all(is.na(x)), logical(1)))
  if (n_failed == n_boot) {
    warning("All G-comp bootstrap replicates failed. SEs and CIs will be NA.",
            call. = FALSE)
  } else if (n_failed > n_boot * 0.5) {
    warning(sprintf("G-comp bootstrap: %d/%d replicates failed. SEs may be unreliable.",
                    n_failed, n_boot), call. = FALSE)
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

#' Bootstrap Inference for TMLE Estimator
#'
#' Resamples subjects (with replacement), re-fits treatment + censoring +
#' outcome models, and re-runs TMLE point estimates. Returns percentile CIs
#' and bootstrap SEs.
#'
#' @param obj longy_data object (fully fitted)
#' @param regime Character regime name
#' @param times Numeric vector of time points
#' @param n_boot Integer number of bootstrap samples
#' @param ci_level Confidence level
#' @param g_bounds Numeric(2). Bounds for cumulative g.
#' @param outcome_range Numeric(2) or NULL. Range for continuous outcome scaling.
#' @param verbose Logical. Print progress.
#' @return data.table with se, ci_lower, ci_upper
#' @noRd
.bootstrap_tmle_inference <- function(obj, regime, times, n_boot = 200L,
                                      ci_level = 0.95, g_bounds = c(0.01, 1),
                                      outcome_range = NULL, verbose = TRUE) {
  nodes <- obj$nodes
  id_col <- nodes$id
  ids <- unique(obj$data[[id_col]])
  n <- length(ids)
  outcome_fit <- obj$fits$outcome[[regime]]

  # Use original censoring column name (not internal binary col names)
  cens_col <- nodes$censoring_col  # NULL if no censoring

  one_boot <- function(b) {
    boot_ids <- sample(ids, n, replace = TRUE)

    # Efficient bootstrap: single data.table join instead of n copies + rbind
    boot_map <- data.table::data.table(
      .boot_orig = boot_ids,
      .boot_new = seq_len(n)
    )
    data.table::setnames(boot_map, ".boot_orig", id_col)
    boot_dt <- merge(boot_map, obj$data, by = id_col, allow.cartesian = TRUE)
    boot_dt[[id_col]] <- boot_dt$.boot_new
    boot_dt[, .boot_new := NULL]

    boot_est <- tryCatch({
      b_obj <- longy_data(
        data = boot_dt,
        id = id_col, time = nodes$time,
        outcome = nodes$outcome, treatment = nodes$treatment,
        censoring = cens_col, observation = nodes$observation,
        baseline = nodes$baseline, timevarying = nodes$timevarying,
        cluster = nodes$cluster, sampling_weights = nodes$sampling_weights,
        outcome_type = nodes$outcome_type,
        competing = nodes$competing, verbose = FALSE
      )
      b_obj$regimes <- obj$regimes

      # Force sequential SL inside bootstrap to prevent nested parallelism
      # (ffSL spawns futures inside future workers → resource explosion)
      # Replay all fit parameters to match original estimation
      trt_fit <- obj$fits$treatment[[regime]]
      b_obj <- fit_treatment(b_obj, regime = regime,
                             covariates = trt_fit$covariates,
                             learners = trt_fit$learners,
                             sl_control = if (!is.null(trt_fit$sl_control)) trt_fit$sl_control else list(),
                             adaptive_cv = if (!is.null(trt_fit$adaptive_cv)) trt_fit$adaptive_cv else TRUE,
                             min_obs = if (!is.null(trt_fit$min_obs)) trt_fit$min_obs else 50L,
                             min_events = if (!is.null(trt_fit$min_events)) trt_fit$min_events else 20L,
                             risk_set = if (!is.null(trt_fit$risk_set)) trt_fit$risk_set else "all",
                             use_ffSL = FALSE,
                             verbose = FALSE)

      # Fit censoring model
      cens_fits_tmle <- obj$fits$censoring[[regime]]
      if (length(cens_fits_tmle) > 0) {
        cov_c <- cens_fits_tmle[[1]]$covariates
        lrn_c <- cens_fits_tmle[[1]]$learners
        slc_c <- if (!is.null(cens_fits_tmle[[1]]$sl_control)) cens_fits_tmle[[1]]$sl_control else list()
        acv_c <- if (!is.null(cens_fits_tmle[[1]]$adaptive_cv)) cens_fits_tmle[[1]]$adaptive_cv else TRUE
        mob_c <- if (!is.null(cens_fits_tmle[[1]]$min_obs)) cens_fits_tmle[[1]]$min_obs else 50L
        mev_c <- if (!is.null(cens_fits_tmle[[1]]$min_events)) cens_fits_tmle[[1]]$min_events else 20L
        b_obj <- fit_censoring(b_obj, regime = regime,
                               covariates = cov_c, learners = lrn_c,
                               sl_control = slc_c,
                               adaptive_cv = acv_c,
                               min_obs = mob_c, min_events = mev_c,
                               use_ffSL = FALSE,
                               verbose = FALSE)
      }

      # Fit observation model (if present in original)
      obs_fit_tmle <- obj$fits$observation[[regime]]
      if (!is.null(obs_fit_tmle)) {
        b_obj <- fit_observation(b_obj, regime = regime,
                                 covariates = obs_fit_tmle$covariates,
                                 learners = obs_fit_tmle$learners,
                                 sl_control = if (!is.null(obs_fit_tmle$sl_control)) obs_fit_tmle$sl_control else list(),
                                 adaptive_cv = if (!is.null(obs_fit_tmle$adaptive_cv)) obs_fit_tmle$adaptive_cv else TRUE,
                                 min_obs = if (!is.null(obs_fit_tmle$min_obs)) obs_fit_tmle$min_obs else 50L,
                                 min_events = if (!is.null(obs_fit_tmle$min_events)) obs_fit_tmle$min_events else 20L,
                                 use_ffSL = FALSE,
                                 verbose = FALSE)
      }

      # Fit outcome model
      b_obj <- fit_outcome(b_obj, regime = regime,
                           covariates = outcome_fit$covariates,
                           learners = outcome_fit$learners,
                           sl_control = if (!is.null(outcome_fit$sl_control)) outcome_fit$sl_control else list(),
                           adaptive_cv = if (!is.null(outcome_fit$adaptive_cv)) outcome_fit$adaptive_cv else TRUE,
                           min_obs = if (!is.null(outcome_fit$min_obs)) outcome_fit$min_obs else 50L,
                           bounds = outcome_fit$bounds,
                           risk_set = if (!is.null(outcome_fit$risk_set)) outcome_fit$risk_set else "all",
                           times = times,
                           use_ffSL = FALSE,
                           verbose = FALSE)

      # Run TMLE (point estimates only)
      tmle_result <- estimate_tmle(b_obj, regime = regime, times = times,
                                   inference = "none", g_bounds = g_bounds,
                                   outcome_range = outcome_range,
                                   verbose = FALSE)

      tmle_key <- paste0(regime, "_tmle")
      tmle_result$results[[tmle_key]]$estimates$estimate
    }, error = function(e) rep(NA_real_, length(times)))

    boot_est
  }

  boot_list <- .run_bootstrap(n_boot, one_boot, verbose)
  boot_estimates <- do.call(rbind, boot_list)

  # Warn if many bootstrap replicates failed
  n_failed <- sum(vapply(boot_list, function(x) all(is.na(x)), logical(1)))
  if (n_failed == n_boot) {
    warning("All TMLE bootstrap replicates failed. SEs and CIs will be NA.",
            call. = FALSE)
  } else if (n_failed > n_boot * 0.5) {
    warning(sprintf("TMLE bootstrap: %d/%d replicates failed. SEs may be unreliable.",
                    n_failed, n_boot), call. = FALSE)
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

#' Hajek point estimate at a single time
#' @noRd
.hajek_estimate <- function(obj, regime, tt) {
  nodes <- obj$nodes
  id_col <- nodes$id
  w_dt <- obj$weights[[regime]]$weights_dt
  w_t <- w_dt[w_dt$.time == tt, ]
  dt_t <- obj$data[obj$data[[nodes$time]] == tt, ]
  merged <- merge(w_t, dt_t[, c(id_col, nodes$outcome), with = FALSE],
                  by = id_col)
  # Filter NA outcomes
  merged <- merged[!is.na(merged[[nodes$outcome]]), ]
  yi <- merged[[nodes$outcome]]
  wi <- merged$.final_weight
  if (length(wi) == 0 || sum(wi) == 0) return(NA_real_)
  stats::weighted.mean(yi, wi)
}
