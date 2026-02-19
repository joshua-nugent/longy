# ===========================================================================
# Simulation Study: G-comp with Continuous Outcome (3 Time Points)
# ===========================================================================
#
# Verifies longy's G-computation estimator (backward ICE) for continuous
# outcomes. Point estimates only (bootstrap too slow for simulation study).
# Reports bias and RMSE.
#
# Two DGP variants toggleable via use_ylag (see sim_ipw_continuous.R for
# full DGP documentation):
#
#   use_ylag = FALSE: L depends on W, L(t-1), A(t-1). No outcome feedback.
#   use_ylag = TRUE:  Y(t-1) feeds back into L, A, Y (coefficients 0.1, 0.1, 0.2).
#
# True ATE without Y_lag: 1.0, 1.2, 1.3
# True ATE with Y_lag:    1.0, 1.45, 1.69
#
# Usage:
#   source("longy/simulations/sim_gcomp_continuous.R")
#   run_all()                    # full simulation (~60 sec, both variants)
#   run_all(n_sim = 10)          # quick smoke test
# ===========================================================================

library(longy)

# === SECTION 1: DGP ===========================================================

generate_data <- function(n, seed, scenario = "alternative",
                          use_ylag = FALSE) {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0
  beta_LY <- if (use_ylag) 0.1 else 0
  beta_AY <- if (use_ylag) 0.1 else 0
  beta_YY <- if (use_ylag) 0.2 else 0

  W <- rnorm(n)
  L_prev <- rep(0, n)
  A_prev <- rep(0L, n)
  Y_true <- rep(0, n)
  alive  <- rep(TRUE, n)

  time_data <- vector("list", 3)

  for (t in 0:2) {
    idx <- which(alive)
    ni  <- length(idx)
    if (ni == 0) break

    L <- numeric(n)
    L[idx] <- 0.5 * W[idx] + 0.5 * L_prev[idx] + beta_LA * A_prev[idx] +
      beta_LY * Y_true[idx] + rnorm(ni, 0, 0.5)

    A <- integer(n)
    A[idx] <- rbinom(ni, 1, plogis(-0.3 + 0.3 * W[idx] + 0.5 * L[idx] +
                                     beta_AY * Y_true[idx]))

    C <- integer(n)
    C[idx] <- rbinom(ni, 1, plogis(-3.5 + 0.3 * L[idx]))

    Y <- rep(NA_real_, n)
    unc <- idx[C[idx] == 0L]
    if (length(unc) > 0) {
      Y[unc] <- 2 + 0.3 * W[unc] + 0.5 * L[unc] +
        beta_YA * A[unc] + beta_YY * Y_true[unc] + rnorm(length(unc))
    }

    R <- integer(n)
    Y_obs <- rep(NA_real_, n)
    if (length(unc) > 0) {
      R[unc] <- rbinom(length(unc), 1, plogis(1.5 - 0.3 * L[unc]))
      Y_obs[unc] <- Y[unc]
      Y_obs[unc[R[unc] == 0L]] <- NA_real_
    }

    time_data[[t + 1]] <- list(
      idx = idx, L = L, A = A, C = C, R = R,
      Y_obs = Y_obs, Y_true = Y
    )

    L_prev <- L
    A_prev <- A
    if (length(unc) > 0) Y_true[unc] <- Y[unc]
    alive[idx[C[idx] == 1L]] <- FALSE
  }

  all_rows <- vector("list", 3)
  last_obs_Y <- rep(0, n)

  for (t in 0:2) {
    td <- time_data[[t + 1]]
    if (is.null(td)) next
    idx <- td$idx

    if (use_ylag) {
      Y_lag <- rep(0, n)
      Y_lag_carried <- rep(0L, n)
      if (t == 0) {
        Y_lag[idx] <- 0
        Y_lag_carried[idx] <- 0L
      } else {
        Y_lag[idx] <- last_obs_Y[idx]
        td_prev <- time_data[[t]]
        prev_idx <- td_prev$idx
        prev_R   <- td_prev$R
        for (i in idx) {
          if (i %in% prev_idx && prev_R[i] == 1L) {
            Y_lag_carried[i] <- 0L
          } else {
            Y_lag_carried[i] <- 1L
          }
        }
      }
      all_rows[[t + 1]] <- data.frame(
        id = idx, time = t,
        W = W[idx], L = td$L[idx],
        Y_lag = Y_lag[idx], Y_lag_carried = Y_lag_carried[idx],
        A = td$A[idx], C = td$C[idx], R = td$R[idx], Y = td$Y_obs[idx],
        stringsAsFactors = FALSE
      )
    } else {
      all_rows[[t + 1]] <- data.frame(
        id = idx, time = t,
        W = W[idx], L = td$L[idx],
        A = td$A[idx], C = td$C[idx], R = td$R[idx], Y = td$Y_obs[idx],
        stringsAsFactors = FALSE
      )
    }

    obs_idx <- idx[td$R[idx] == 1L & td$C[idx] == 0L]
    if (length(obs_idx) > 0) last_obs_Y[obs_idx] <- td$Y_true[obs_idx]
  }

  do.call(rbind, all_rows)
}


# === SECTION 2: TRUE VALUES ====================================================

true_values <- function(scenario = "alternative", use_ylag = FALSE) {
  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0
  beta_LY <- if (use_ylag) 0.1 else 0
  beta_YY <- if (use_ylag) 0.2 else 0

  out <- list()
  for (regime in c("always", "never")) {
    a_val <- if (regime == "always") 1 else 0
    EL <- numeric(3)
    EY <- numeric(3)

    EL[1] <- 0
    EY[1] <- 2 + 0.5 * EL[1] + beta_YA * a_val

    for (k in 2:3) {
      EL[k] <- 0.5 * EL[k - 1] + beta_LA * a_val + beta_LY * EY[k - 1]
      EY[k] <- 2 + 0.5 * EL[k] + beta_YA * a_val + beta_YY * EY[k - 1]
    }

    out[[regime]] <- data.frame(regime = regime, time = 0:2, true_mean = EY)
  }
  do.call(rbind, out)
}


# === SECTION 3: SIMULATION RUNNER ==============================================

run_simulation <- function(n = 500, n_sim = 500, scenario = "alternative",
                           seed_start = 1, use_ylag = FALSE) {
  ylag_label <- if (use_ylag) "with Y_lag" else "without Y_lag"
  cat(sprintf("\nRunning %d G-comp simulations (n=%d, %s, %s)\n",
              n_sim, n, scenario, ylag_label))
  t0 <- proc.time()

  tv <- if (use_ylag) c("L", "Y_lag", "Y_lag_carried") else "L"
  all_results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    d <- generate_data(n = n, seed = seed_start + i - 1,
                       scenario = scenario, use_ylag = use_ylag)

    res <- tryCatch(
      suppressWarnings(longy(
        data = d,
        id = "id", time = "time", outcome = "Y",
        treatment = "A", censoring = "C", observation = "R",
        baseline = "W",
        timevarying = tv,
        outcome_type = "continuous",
        regimes = list(always = 1L, never = 0L),
        estimator = "gcomp",
        learners = NULL,
        n_boot = 0L,
        verbose = FALSE
      )),
      error = function(e) {
        warning(sprintf("Sim %d failed: %s", i, e$message))
        NULL
      }
    )

    if (is.null(res)) next

    sim_rows <- list()
    for (regime in c("always", "never")) {
      est <- res[[regime]]$estimates
      sim_rows[[regime]] <- data.frame(
        sim      = i,
        regime   = regime,
        time     = est$time,
        estimate = est$estimate,
        stringsAsFactors = FALSE
      )
    }
    all_results[[i]] <- do.call(rbind, sim_rows)

    if (i %% 50 == 0) {
      elapsed <- (proc.time() - t0)[3]
      rate    <- elapsed / i
      cat(sprintf("  %d/%d done (%.1f s elapsed, ~%.0f s remaining)\n",
                  i, n_sim, elapsed, rate * (n_sim - i)))
    }
  }

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  Completed in %.1f seconds (%.2f s/sim)\n", elapsed,
              elapsed / n_sim))
  do.call(rbind, all_results)
}


# === SECTION 4: SUMMARIZE =====================================================

summarize_results <- function(sim_results, truth) {
  merged <- merge(sim_results, truth, by = c("regime", "time"))

  summ <- do.call(rbind, lapply(
    split(merged, list(merged$regime, merged$time), drop = TRUE),
    function(x) {
      data.frame(
        regime    = x$regime[1],
        time      = x$time[1],
        true_mean = x$true_mean[1],
        mean_est  = mean(x$estimate),
        bias      = mean(x$estimate) - x$true_mean[1],
        rmse      = sqrt(mean((x$estimate - x$true_mean[1])^2)),
        emp_se    = sd(x$estimate),
        n_sims    = nrow(x),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(summ) <- NULL
  summ[order(summ$regime, summ$time), ]
}

summarize_ate <- function(sim_results, truth) {
  a <- sim_results[sim_results$regime == "always", c("sim", "time", "estimate")]
  n <- sim_results[sim_results$regime == "never",  c("sim", "time", "estimate")]
  names(a)[3] <- "estimate_a"
  names(n)[3] <- "estimate_n"

  ate <- merge(a, n, by = c("sim", "time"))
  ate$ate_est <- ate$estimate_a - ate$estimate_n

  ta <- truth$true_mean[truth$regime == "always"]
  tn <- truth$true_mean[truth$regime == "never"]
  true_ate <- data.frame(time = 0:2, true_ate = ta - tn)
  ate <- merge(ate, true_ate, by = "time")

  do.call(rbind, lapply(split(ate, ate$time), function(x) {
    data.frame(
      time     = x$time[1],
      true_ate = x$true_ate[1],
      mean_ate = mean(x$ate_est),
      bias     = mean(x$ate_est) - x$true_ate[1],
      rmse     = sqrt(mean((x$ate_est - x$true_ate[1])^2)),
      emp_se   = sd(x$ate_est),
      n_sims   = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
}

print_summary <- function(scenario_label, regime_summ, ate_summ) {
  cat(sprintf("\n%s\n%s\n", scenario_label, strrep("=", nchar(scenario_label))))

  cat("\nPer-regime estimates:\n")
  cat(sprintf("  %-8s %-6s %-8s %-8s %-+8s %-8s %-8s\n",
              "Regime", "Time", "Truth", "Mean", "Bias", "RMSE", "EmpSE"))
  cat(sprintf("  %s\n", strrep("-", 56)))
  for (i in seq_len(nrow(regime_summ))) {
    r <- regime_summ[i, ]
    cat(sprintf("  %-8s %-6d %-8.3f %-8.3f %-+8.4f %-8.4f %-8.4f\n",
                r$regime, r$time, r$true_mean, r$mean_est, r$bias,
                r$rmse, r$emp_se))
  }

  cat("\nATE (always - never):\n")
  cat(sprintf("  %-6s %-8s %-8s %-+8s %-8s %-8s\n",
              "Time", "Truth", "Mean", "Bias", "RMSE", "EmpSE"))
  cat(sprintf("  %s\n", strrep("-", 48)))
  for (i in seq_len(nrow(ate_summ))) {
    r <- ate_summ[i, ]
    cat(sprintf("  %-6d %-8.3f %-8.3f %-+8.4f %-8.4f %-8.4f\n",
                r$time, r$true_ate, r$mean_ate, r$bias, r$rmse, r$emp_se))
  }
}


# === SECTION 5: RUN ============================================================

run_all <- function(n = 500, n_sim = 500, seed_start = 1) {
  all_out <- list()

  for (use_ylag in c(FALSE, TRUE)) {
    ylag_tag <- if (use_ylag) "ylag" else "no_ylag"
    ylag_label <- if (use_ylag) "WITH Y_lag" else "WITHOUT Y_lag"

    truth_alt <- true_values("alternative", use_ylag = use_ylag)
    sim_alt   <- run_simulation(n, n_sim, "alternative", seed_start,
                                use_ylag = use_ylag)
    regime_alt <- summarize_results(sim_alt, truth_alt)
    ate_alt    <- summarize_ate(sim_alt, truth_alt)

    ate_vals <- truth_alt$true_mean[truth_alt$regime == "always"] -
      truth_alt$true_mean[truth_alt$regime == "never"]
    ate_label <- paste(sprintf("%.2f", ate_vals), collapse = ", ")
    print_summary(
      sprintf("ALTERNATIVE — %s (ATE = %s)", ylag_label, ate_label),
      regime_alt, ate_alt
    )

    truth_null <- true_values("null", use_ylag = use_ylag)
    sim_null   <- run_simulation(n, n_sim, "null",
                                 seed_start + n_sim,
                                 use_ylag = use_ylag)
    regime_null <- summarize_results(sim_null, truth_null)
    ate_null    <- summarize_ate(sim_null, truth_null)
    print_summary(
      sprintf("NULL — %s (ATE = 0 at all times)", ylag_label),
      regime_null, ate_null
    )

    cat(sprintf("\n--- Diagnostics [%s] ---\n", ylag_label))
    cat(sprintf("Alt: max |bias| = %.4f, max RMSE = %.4f\n",
                max(abs(regime_alt$bias)), max(regime_alt$rmse)))
    cat(sprintf("Null: max |ATE bias| = %.4f\n",
                max(abs(ate_null$bias))))

    all_out[[ylag_tag]] <- list(
      alternative = list(sim = sim_alt, regime = regime_alt,
                         ate = ate_alt, truth = truth_alt),
      null        = list(sim = sim_null, regime = regime_null,
                         ate = ate_null, truth = truth_null)
    )

    seed_start <- seed_start + 2 * n_sim
  }

  # Cross-variant comparison
  cat("\n\n")
  cat("==================================================\n")
  cat("  COMPARISON: Y_lag effect on bias and RMSE\n")
  cat("==================================================\n")
  cat(sprintf("  %-10s %-6s  %-10s %-10s  %-10s %-10s\n",
              "", "", "No Y_lag", "", "With Y_lag", ""))
  cat(sprintf("  %-10s %-6s  %-10s %-10s  %-10s %-10s\n",
              "Scenario", "Time", "Bias", "RMSE", "Bias", "RMSE"))
  cat(sprintf("  %s\n", strrep("-", 62)))
  for (scenario in c("alternative", "null")) {
    an <- all_out$no_ylag[[scenario]]$ate
    ay <- all_out$ylag[[scenario]]$ate
    label <- if (scenario == "alternative") "Alt ATE" else "Null ATE"
    for (t in 0:2) {
      bn <- an$bias[an$time == t]; rn <- an$rmse[an$time == t]
      by <- ay$bias[ay$time == t]; ry <- ay$rmse[ay$time == t]
      cat(sprintf("  %-10s %-6d  %-+10.4f %-10.4f  %-+10.4f %-10.4f\n",
                  label, t, bn, rn, by, ry))
    }
    cat(sprintf("  %s\n", strrep("-", 62)))
  }

  invisible(all_out)
}
