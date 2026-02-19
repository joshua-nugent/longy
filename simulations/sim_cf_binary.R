# ===========================================================================
# Simulation Study: Cross-Fitting vs Non-CF TMLE — Binary Outcome
# ===========================================================================
#
# Complex DGP with genuine nonlinearities to demonstrate cross-fitting's
# value when SuperLearner with aggressive learners is used. GLM is
# meaningfully misspecified due to:
#   - W1^2 in L1, A, C, outcome
#   - W1*W3 interaction in L1
#   - A*L1 interaction in L1, outcome
#   - sin(W3*pi) in L2, outcome
#   - I(L1>0) threshold in L2
#   - L1^2 in treatment model
#
# 4-arm comparison:
#   1. GLM (baseline — misspecified)
#   2. GLM + 5-fold CF (sample-splitting benefit alone)
#   3. SL (flexible learners without protection)
#   4. SL + 5-fold CF (recommended approach)
#
# Usage:
#   source("longy/simulations/sim_cf_binary.R")
#   run_all()                    # full simulation (both scenarios)
#   run_all(n_sim = 10)          # quick smoke test
# ===========================================================================

library(longy)

# === SECTION 1: DGP ===========================================================

generate_data <- function(n, seed, scenario = "alternative") {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 0.8

  # Baseline covariates
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.3)
  W3 <- runif(n, -1, 1)

  L1_prev <- rep(0, n)
  L2_prev <- rep(0L, n)
  A_prev  <- rep(0L, n)
  alive   <- rep(TRUE, n)

  time_data <- vector("list", 5)

  for (t in 0:4) {
    idx <- which(alive)
    ni  <- length(idx)
    if (ni == 0) break

    # --- L1 (continuous, AR-1 with nonlinear effects) ---
    L1 <- numeric(n)
    if (t == 0) {
      L1[idx] <- 0.5 * W1[idx] + 0.3 * W1[idx]^2 - 0.5 * W2[idx] +
        rnorm(ni, 0, 0.5)
    } else {
      L1[idx] <- 0.3 * W1[idx] + 0.2 * W1[idx] * W3[idx] +
        0.4 * L1_prev[idx] + beta_LA * A_prev[idx] +
        0.2 * beta_LA * A_prev[idx] * L1_prev[idx] +
        rnorm(ni, 0, 0.5)
    }

    # --- L2 (binary, with threshold and trigonometric terms) ---
    L2 <- integer(n)
    if (t == 0) {
      L2[idx] <- rbinom(ni, 1, plogis(-0.5 + 0.3 * W1[idx] +
                                        0.5 * W2[idx] +
                                        0.3 * sin(W3[idx] * pi)))
    } else {
      L2[idx] <- rbinom(ni, 1, plogis(-0.3 + 0.3 * W1[idx] +
                                        0.4 * L1[idx] +
                                        0.3 * as.integer(L1[idx] > 0) +
                                        0.2 * A_prev[idx]))
    }

    # --- Treatment ---
    A <- integer(n)
    A[idx] <- rbinom(ni, 1, plogis(-0.5 + 0.2 * W1[idx] +
                                    0.3 * W2[idx] +
                                    0.4 * L1[idx] +
                                    0.2 * L1[idx]^2 +
                                    0.3 * L2[idx] +
                                    0.2 * W3[idx] * L1[idx]))

    # --- Censoring (absorbing, ~3-5% per time) ---
    C <- integer(n)
    C[idx] <- rbinom(ni, 1, plogis(-3.5 + 0.2 * L1[idx] +
                                    0.3 * L2[idx] +
                                    0.1 * W1[idx]^2))

    # --- Outcome (binary) ---
    Y <- rep(NA_integer_, n)
    unc <- idx[C[idx] == 0L]
    if (length(unc) > 0) {
      lp_Y <- -1.5 + 0.3 * W1[unc] + 0.2 * W1[unc]^2 +
        0.4 * L1[unc] + 0.3 * L2[unc] +
        beta_YA * A[unc] +
        0.3 * sin(W3[unc] * pi)
      if (scenario != "null") {
        lp_Y <- lp_Y + 0.2 * A[unc] * L1[unc]
      }
      Y[unc] <- rbinom(length(unc), 1, plogis(lp_Y))
    }

    # --- Observation (intermittent, ~80%) ---
    R <- integer(n)
    Y_obs <- rep(NA_integer_, n)
    if (length(unc) > 0) {
      R[unc] <- rbinom(length(unc), 1, plogis(1.5 - 0.3 * L1[unc] -
                                                0.2 * L2[unc] +
                                                0.1 * W2[unc]))
      Y_obs[unc] <- Y[unc]
      Y_obs[unc[R[unc] == 0L]] <- NA_integer_
    }

    time_data[[t + 1]] <- list(
      idx = idx, L1 = L1, L2 = L2, A = A, C = C, R = R,
      Y_obs = Y_obs, Y_true = Y
    )

    L1_prev <- L1
    L2_prev <- L2
    A_prev  <- A
    alive[idx[C[idx] == 1L]] <- FALSE
  }

  # Assemble long-format data.frame
  all_rows <- vector("list", 5)
  for (t in 0:4) {
    td <- time_data[[t + 1]]
    if (is.null(td)) next
    idx <- td$idx
    all_rows[[t + 1]] <- data.frame(
      id = idx, time = t,
      W1 = W1[idx], W2 = W2[idx], W3 = W3[idx],
      L1 = td$L1[idx], L2 = td$L2[idx],
      A = td$A[idx], C = td$C[idx], R = td$R[idx], Y = td$Y_obs[idx],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, all_rows)
}


# === SECTION 2: TRUE VALUES ====================================================

true_values <- function(scenario = "alternative", n_mc = 1e6, seed = 42) {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 0.8

  W1 <- rnorm(n_mc)
  W2 <- rbinom(n_mc, 1, 0.3)
  W3 <- runif(n_mc, -1, 1)

  out <- list()
  for (regime in c("always", "never")) {
    a_val <- if (regime == "always") 1 else 0

    L1_prev <- rep(0, n_mc)
    A_prev  <- rep(0, n_mc)

    means <- numeric(5)
    for (t in 0:4) {
      # L1
      if (t == 0) {
        L1 <- 0.5 * W1 + 0.3 * W1^2 - 0.5 * W2 + rnorm(n_mc, 0, 0.5)
      } else {
        L1 <- 0.3 * W1 + 0.2 * W1 * W3 + 0.4 * L1_prev +
          beta_LA * A_prev + 0.2 * beta_LA * A_prev * L1_prev +
          rnorm(n_mc, 0, 0.5)
      }

      # L2 (draw for feedback but not needed for Y)
      # (included for completeness — L2 is binary so we draw it)
      if (t == 0) {
        L2 <- rbinom(n_mc, 1, plogis(-0.5 + 0.3 * W1 + 0.5 * W2 +
                                       0.3 * sin(W3 * pi)))
      } else {
        L2 <- rbinom(n_mc, 1, plogis(-0.3 + 0.3 * W1 + 0.4 * L1 +
                                       0.3 * as.integer(L1 > 0) +
                                       0.2 * A_prev))
      }

      # E[Y(t)] under regime — use plogis for lower MC variance
      lp_Y <- -1.5 + 0.3 * W1 + 0.2 * W1^2 + 0.4 * L1 + 0.3 * L2 +
        beta_YA * a_val + 0.3 * sin(W3 * pi)
      if (scenario != "null") {
        lp_Y <- lp_Y + 0.2 * a_val * L1
      }
      means[t + 1] <- mean(plogis(lp_Y))

      L1_prev <- L1
      A_prev  <- rep(a_val, n_mc)
    }

    out[[regime]] <- data.frame(regime = regime, time = 0:4,
                                true_mean = means,
                                stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

mc_verify_true_values <- function(scenario = "alternative", n_mc = 1e6) {
  tv1 <- true_values(scenario, n_mc = n_mc, seed = 42)
  tv2 <- true_values(scenario, n_mc = n_mc, seed = 123)
  max_diff <- max(abs(tv1$true_mean - tv2$true_mean))
  cat(sprintf("MC verification (%s): max |diff| across seeds = %.6f\n",
              scenario, max_diff))
  if (max_diff > 0.005) {
    warning("MC true values may not be stable — consider increasing n_mc")
  }
  invisible(data.frame(time = tv1$time, regime = tv1$regime,
                       seed42 = tv1$true_mean, seed123 = tv2$true_mean,
                       diff = tv1$true_mean - tv2$true_mean))
}


# === SECTION 3: SIMULATION RUNNER ==============================================

run_simulation <- function(n = 1000, n_sim = 500, scenario = "alternative",
                           seed_start = 1, n_cores = 10) {
  cat(sprintf("\nRunning %d CF-comparison TMLE simulations (binary, n=%d, %s, %d cores)\n",
              n_sim, n, scenario, n_cores))
  t0 <- proc.time()

  sl_library <- c("SL.glm", "SL.gam", "SL.earth", "SL.nnet")

  arms <- list(
    list(label = "glm",    learners = NULL,       cross_fit = NULL),
    list(label = "glm_cf", learners = NULL,       cross_fit = 5L),
    list(label = "sl",     learners = sl_library, cross_fit = NULL),
    list(label = "sl_cf",  learners = sl_library, cross_fit = 5L)
  )

  run_one_sim <- function(sim_id) {
    d <- generate_data(n = n, seed = seed_start + sim_id - 1,
                       scenario = scenario)

    sim_results <- vector("list", length(arms))
    for (a in seq_along(arms)) {
      arm <- arms[[a]]

      res <- tryCatch(
        suppressWarnings(longy(
          data = d,
          id = "id", time = "time", outcome = "Y",
          treatment = "A", censoring = "C", observation = "R",
          baseline = c("W1", "W2", "W3"),
          timevarying = c("L1", "L2"),
          outcome_type = "binary",
          regimes = list(always = 1L, never = 0L),
          estimator = "tmle",
          learners = arm$learners,
          cross_fit = arm$cross_fit,
          sl_fn = "SuperLearner",
          n_boot = 0L,
          verbose = FALSE
        )),
        error = function(e) {
          warning(sprintf("Sim %d arm %s failed: %s",
                          sim_id, arm$label, e$message))
          NULL
        }
      )

      if (is.null(res)) next

      sim_rows <- list()
      for (regime in c("always", "never")) {
        est <- res[[regime]]$estimates
        sim_rows[[regime]] <- data.frame(
          sim      = sim_id,
          arm      = arm$label,
          regime   = regime,
          time     = est$time,
          estimate = est$estimate,
          se       = if ("se" %in% names(est)) est$se else NA_real_,
          ci_lower = if ("ci_lower" %in% names(est)) est$ci_lower else NA_real_,
          ci_upper = if ("ci_upper" %in% names(est)) est$ci_upper else NA_real_,
          stringsAsFactors = FALSE
        )
      }
      sim_results[[a]] <- do.call(rbind, sim_rows)
    }

    do.call(rbind, sim_results)
  }

  results <- parallel::mclapply(seq_len(n_sim), run_one_sim,
                                mc.cores = n_cores)

  elapsed <- (proc.time() - t0)[3]
  n_ok <- sum(vapply(results, function(x) !is.null(x), logical(1)))
  cat(sprintf("  Completed %d/%d in %.1f seconds (%.2f s/sim)\n",
              n_ok, n_sim, elapsed, elapsed / n_sim))

  do.call(rbind, results)
}


# === SECTION 4: SUMMARIZE =====================================================

summarize_results <- function(sim_results, truth) {
  merged <- merge(sim_results, truth, by = c("regime", "time"))

  summ <- do.call(rbind, lapply(
    split(merged, list(merged$arm, merged$regime, merged$time), drop = TRUE),
    function(x) {
      covers <- !is.na(x$ci_lower) & !is.na(x$ci_upper) &
        x$ci_lower <= x$true_mean & x$ci_upper >= x$true_mean
      emp_se <- sd(x$estimate)
      mean_se <- mean(x$se, na.rm = TRUE)
      data.frame(
        arm       = x$arm[1],
        regime    = x$regime[1],
        time      = x$time[1],
        true_mean = x$true_mean[1],
        mean_est  = mean(x$estimate),
        bias      = mean(x$estimate) - x$true_mean[1],
        rmse      = sqrt(mean((x$estimate - x$true_mean[1])^2)),
        emp_se    = emp_se,
        mean_se   = mean_se,
        se_ratio  = if (emp_se > 0) mean_se / emp_se else NA_real_,
        coverage  = mean(covers),
        n_sims    = nrow(x),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(summ) <- NULL
  summ[order(summ$arm, summ$regime, summ$time), ]
}

summarize_ate <- function(sim_results, truth) {
  a_res <- sim_results[sim_results$regime == "always",
                       c("sim", "arm", "time", "estimate", "se",
                         "ci_lower", "ci_upper")]
  n_res <- sim_results[sim_results$regime == "never",
                       c("sim", "arm", "time", "estimate", "se",
                         "ci_lower", "ci_upper")]
  names(a_res)[4:7] <- paste0(names(a_res)[4:7], "_a")
  names(n_res)[4:7] <- paste0(names(n_res)[4:7], "_n")

  ate <- merge(a_res, n_res, by = c("sim", "arm", "time"))
  ate$ate_est <- ate$estimate_a - ate$estimate_n
  ate$ate_se  <- sqrt(ate$se_a^2 + ate$se_n^2)
  z <- qnorm(0.975)
  ate$ate_ci_lo <- ate$ate_est - z * ate$ate_se
  ate$ate_ci_hi <- ate$ate_est + z * ate$ate_se

  ta <- truth$true_mean[truth$regime == "always"]
  tn <- truth$true_mean[truth$regime == "never"]
  true_ate <- data.frame(time = 0:4, true_ate = ta - tn)
  ate <- merge(ate, true_ate, by = "time")

  do.call(rbind, lapply(
    split(ate, list(ate$arm, ate$time), drop = TRUE),
    function(x) {
      covers <- !is.na(x$ate_ci_lo) & !is.na(x$ate_ci_hi) &
        x$ate_ci_lo <= x$true_ate & x$ate_ci_hi >= x$true_ate
      emp_se <- sd(x$ate_est)
      mean_se <- mean(x$ate_se, na.rm = TRUE)
      data.frame(
        arm      = x$arm[1],
        time     = x$time[1],
        true_ate = x$true_ate[1],
        mean_ate = mean(x$ate_est),
        bias     = mean(x$ate_est) - x$true_ate[1],
        rmse     = sqrt(mean((x$ate_est - x$true_ate[1])^2)),
        emp_se   = emp_se,
        mean_se  = mean_se,
        se_ratio = if (emp_se > 0) mean_se / emp_se else NA_real_,
        coverage = mean(covers),
        n_sims   = nrow(x),
        stringsAsFactors = FALSE
      )
    }
  ))
}

print_summary <- function(scenario_label, regime_summ, ate_summ) {
  cat(sprintf("\n%s\n%s\n", scenario_label, strrep("=", nchar(scenario_label))))

  arm_labels <- c("glm", "glm_cf", "sl", "sl_cf")
  arm_display <- c("GLM", "GLM+CF", "SL", "SL+CF")

  cat("\nATE (always - never) by arm:\n")
  cat(sprintf("  %s\n", strrep("-", 100)))
  cat(sprintf("  %-8s %-6s", "Arm", "Time"))
  cat(sprintf("  %-8s %-+9s %-8s %-8s %-8s %-7s %-7s\n",
              "Truth", "Bias", "RMSE", "EmpSE", "MeanSE", "SE_rat", "Cover"))
  cat(sprintf("  %s\n", strrep("-", 100)))

  for (a in seq_along(arm_labels)) {
    sub <- ate_summ[ate_summ$arm == arm_labels[a], ]
    if (nrow(sub) == 0) next
    sub <- sub[order(sub$time), ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      cat(sprintf("  %-8s %-6d  %-8.3f %-+9.4f %-8.4f %-8.4f %-8.4f %-7.2f %-7.3f\n",
                  arm_display[a], r$time, r$true_ate, r$bias, r$rmse,
                  r$emp_se, r$mean_se, r$se_ratio, r$coverage))
    }
    cat(sprintf("  %s\n", strrep("-", 100)))
  }
}

print_comparison <- function(ate_summ) {
  cat("\n")
  cat("==================================================================\n")
  cat("  CROSS-ARM COMPARISON (ATE)\n")
  cat("==================================================================\n")
  cat(sprintf("  %-6s", "Time"))
  for (arm in c("GLM", "GLM+CF", "SL", "SL+CF")) {
    cat(sprintf("  %-+7s %-6s", paste0(arm, ":B"), "Cov"))
  }
  cat("\n")
  cat(sprintf("  %s\n", strrep("-", 70)))

  arm_labels <- c("glm", "glm_cf", "sl", "sl_cf")
  for (t in 0:4) {
    cat(sprintf("  %-6d", t))
    for (a in arm_labels) {
      sub <- ate_summ[ate_summ$arm == a & ate_summ$time == t, ]
      if (nrow(sub) == 0) {
        cat(sprintf("  %7s %6s", "NA", "NA"))
      } else {
        cat(sprintf("  %-+7.3f %-6.3f", sub$bias, sub$coverage))
      }
    }
    cat("\n")
  }
  cat(sprintf("  %s\n", strrep("-", 70)))
}


# === SECTION 5: RUN ============================================================

run_all <- function(n = 1000, n_sim = 500, seed_start = 1, n_cores = 10) {
  all_out <- list()

  # --- MC verification ---
  cat("Verifying MC true values stability...\n")
  mc_verify_true_values("alternative")
  mc_verify_true_values("null")
  cat("\n")

  # --- ALTERNATIVE SCENARIO ---
  truth_alt <- true_values("alternative")
  sim_alt   <- run_simulation(n, n_sim, "alternative", seed_start, n_cores)
  regime_alt <- summarize_results(sim_alt, truth_alt)
  ate_alt    <- summarize_ate(sim_alt, truth_alt)

  ate_vals <- truth_alt$true_mean[truth_alt$regime == "always"] -
    truth_alt$true_mean[truth_alt$regime == "never"]
  ate_label <- paste(sprintf("%.3f", ate_vals), collapse = ", ")
  print_summary(
    sprintf("ALTERNATIVE (binary) -- ATE ~ %s", ate_label),
    regime_alt, ate_alt
  )
  print_comparison(ate_alt)

  # --- NULL SCENARIO ---
  truth_null <- true_values("null")
  sim_null   <- run_simulation(n, n_sim, "null", seed_start + n_sim, n_cores)
  regime_null <- summarize_results(sim_null, truth_null)
  ate_null    <- summarize_ate(sim_null, truth_null)
  print_summary(
    "NULL (binary) -- ATE = 0 at all times",
    regime_null, ate_null
  )
  print_comparison(ate_null)

  # --- DIAGNOSTICS ---
  cat("\n--- Diagnostics ---\n")
  for (arm in c("glm", "glm_cf", "sl", "sl_cf")) {
    sub_alt  <- ate_alt[ate_alt$arm == arm, ]
    sub_null <- ate_null[ate_null$arm == arm, ]
    if (nrow(sub_alt) > 0 && nrow(sub_null) > 0) {
      cat(sprintf("  %s: Alt max|bias|=%.4f, coverage=[%s] | Null type-I=[%s]\n",
                  toupper(arm),
                  max(abs(sub_alt$bias)),
                  paste(sprintf("%.3f", sub_alt$coverage), collapse=","),
                  paste(sprintf("%.3f", 1 - sub_null$coverage), collapse=",")))
    }
  }

  all_out <- list(
    alternative = list(sim = sim_alt, regime = regime_alt,
                       ate = ate_alt, truth = truth_alt),
    null        = list(sim = sim_null, regime = regime_null,
                       ate = ate_null, truth = truth_null)
  )

  invisible(all_out)
}
