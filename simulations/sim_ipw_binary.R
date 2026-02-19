# ===========================================================================
# Simulation Study: IPW with Binary Outcome (3 Time Points)
# ===========================================================================
#
# Verifies longy's IPW estimator under a DGP with time-varying confounding
# and a binary outcome. Two DGP variants are toggleable via use_ylag:
#
#   use_ylag = FALSE: Simple DGP. L depends on W, L(t-1), A(t-1). A depends
#     on W, L(t). Y depends on W, L(t), A(t). No outcome feedback.
#
#   use_ylag = TRUE: Y(t-1) feeds back into L(t), A(t), and Y(t). When R
#     masks Y, LOCF carries the last observed Y forward, with a missingness
#     indicator. Creates richer confounding that stresses the IPW estimator.
#
# Structural equations:
#   W ~ N(0, 1)
#   L(t) = 0.5*W + 0.5*L(t-1) + 0.4*A(t-1) [+ 0.3*Y(t-1)] + N(0, 0.25)
#   A(t) ~ Bern(expit(-0.3 + 0.3*W + 0.5*L(t) [+ 0.3*Y(t-1)]))
#   C(t) ~ Bern(expit(-3.5 + 0.3*L(t)))         [~3-5% per time]
#   Y(t) ~ Bern(expit(-1 + 0.3*W + 0.5*L(t) + 1.0*A(t) [+ 0.5*Y(t-1)]))
#   R(t) ~ Bern(expit(1.5 - 0.3*L(t)))           [~82% observation rate]
#
# Y(-1) = 0, A(-1) = 0, L(-1) = 0 by convention.
#
# Null scenario: beta_YA = 0, beta_LA = 0 (treatment has no effect).
# Y_lag feedback coefficients (if enabled) remain in the null.
#
# True counterfactual means are computed by Monte Carlo (1e6 subjects,
# intervene on A, no censoring or missingness, use plogis for precision).
# Without Y_lag, no Y draws are needed (simpler MC, less variance).
#
# Usage:
#   source("longy/simulations/sim_ipw_binary.R")
#   run_all()                    # full simulation (~60 sec, both variants)
#   run_all(n_sim = 10)          # quick smoke test
#   mc_verify_true_values()      # verify MC truth is stable
# ===========================================================================

library(longy)

# === SECTION 1: DGP ===========================================================

#' Generate one dataset from the binary outcome DGP
#'
#' @param n Integer. Number of subjects.
#' @param seed Integer. Random seed.
#' @param scenario Character. "alternative" or "null".
#' @param use_ylag Logical. Include Y_lag feedback in L, A, and Y equations?
#' @return data.frame in longy long format.
generate_data <- function(n, seed, scenario = "alternative",
                          use_ylag = FALSE) {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0
  beta_LY <- if (use_ylag) 0.3 else 0   # Y(t-1) -> L(t)
  beta_AY <- if (use_ylag) 0.3 else 0   # Y(t-1) -> A(t)
  beta_YY <- if (use_ylag) 0.5 else 0   # Y(t-1) -> Y(t)

  W <- rnorm(n)
  L_prev   <- rep(0, n)       # L(-1) = 0
  A_prev   <- rep(0L, n)      # A(-1) = 0
  Y_true   <- rep(0L, n)      # Y(-1) = 0
  alive    <- rep(TRUE, n)

  time_data <- vector("list", 3)

  for (t in 0:2) {
    idx <- which(alive)
    ni  <- length(idx)
    if (ni == 0) break

    # L(t) = 0.5*W + 0.5*L(t-1) + beta_LA*A(t-1) + beta_LY*Y(t-1) + eps_L
    L <- numeric(n)
    L[idx] <- 0.5 * W[idx] + 0.5 * L_prev[idx] + beta_LA * A_prev[idx] +
      beta_LY * Y_true[idx] + rnorm(ni, 0, 0.5)

    # A(t) ~ Bern(expit(-0.3 + 0.3*W + 0.5*L(t) + beta_AY*Y(t-1)))
    A <- integer(n)
    A[idx] <- rbinom(ni, 1, plogis(-0.3 + 0.3 * W[idx] + 0.5 * L[idx] +
                                     beta_AY * Y_true[idx]))

    # C(t) ~ Bern(expit(-3.5 + 0.3*L(t)))
    C <- integer(n)
    C[idx] <- rbinom(ni, 1, plogis(-3.5 + 0.3 * L[idx]))

    # Y(t) ~ Bern(expit(-1 + 0.3*W + 0.5*L(t) + beta_YA*A(t) + beta_YY*Y(t-1)))
    Y <- rep(NA_integer_, n)
    unc <- idx[C[idx] == 0L]
    if (length(unc) > 0) {
      Y[unc] <- rbinom(length(unc), 1,
                        plogis(-1 + 0.3 * W[unc] + 0.5 * L[unc] +
                                 beta_YA * A[unc] + beta_YY * Y_true[unc]))
    }

    # R(t) ~ Bern(expit(1.5 - 0.3*L(t)))
    R <- integer(n)
    Y_obs <- rep(NA_integer_, n)
    if (length(unc) > 0) {
      R[unc] <- rbinom(length(unc), 1, plogis(1.5 - 0.3 * L[unc]))
      Y_obs[unc] <- Y[unc]
      Y_obs[unc[R[unc] == 0L]] <- NA_integer_
    }

    time_data[[t + 1]] <- list(
      idx = idx, L = L, A = A, C = C, R = R,
      Y_obs = Y_obs, Y_true = Y
    )

    # Update state (use TRUE Y, not observed)
    L_prev <- L
    A_prev <- A
    if (length(unc) > 0) Y_true[unc] <- Y[unc]
    alive[idx[C[idx] == 1L]] <- FALSE
  }

  # --- Build output (with optional Y_lag via LOCF) ---
  all_rows <- vector("list", 3)
  last_obs_Y <- rep(0L, n)  # Y(-1) = 0

  for (t in 0:2) {
    td <- time_data[[t + 1]]
    if (is.null(td)) next
    idx <- td$idx

    if (use_ylag) {
      Y_lag <- rep(0L, n)
      Y_lag_carried <- rep(0L, n)

      if (t == 0) {
        Y_lag[idx] <- 0L
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

    # Update LOCF tracker
    obs_idx <- idx[td$R[idx] == 1L & td$C[idx] == 0L]
    if (length(obs_idx) > 0) {
      last_obs_Y[obs_idx] <- td$Y_true[obs_idx]
    }
  }

  do.call(rbind, all_rows)
}


# === SECTION 2: TRUE VALUES ====================================================

#' Monte Carlo true counterfactual means for binary outcome
#'
#' With nonlinear (logistic) structural equations, there is no closed-form
#' expression for the marginal mean. We use a large Monte Carlo sample under
#' each intervention, with no censoring or missingness. For precision, we
#' average plogis(LP) instead of drawing binary outcomes.
#'
#' Without Y_lag, no binary Y draws are needed at all (L doesn't depend on Y),
#' so the MC has lower variance. With Y_lag, binary draws are needed for
#' forward simulation of Y feedback through L, A, and Y.
#'
#' @param scenario "alternative" or "null"
#' @param use_ylag Logical.
#' @param n_mc Integer. Monte Carlo sample size (default 1e6).
#' @param seed Random seed.
#' @return data.frame with columns: regime, time, true_mean
true_values <- function(scenario = "alternative", use_ylag = FALSE,
                        n_mc = 1e6, seed = 42) {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0
  beta_LY <- if (use_ylag) 0.3 else 0
  beta_AY <- if (use_ylag) 0.3 else 0
  beta_YY <- if (use_ylag) 0.5 else 0

  W <- rnorm(n_mc)

  out <- list()
  for (regime in c("always", "never")) {
    a_val <- if (regime == "always") 1 else 0

    L_prev <- rep(0, n_mc)
    A_prev <- rep(0, n_mc)
    Y_prev <- rep(0, n_mc)   # Y(-1) = 0

    means <- numeric(3)
    for (t in 0:2) {
      L <- 0.5 * W + 0.5 * L_prev + beta_LA * A_prev +
        beta_LY * Y_prev + rnorm(n_mc, 0, 0.5)

      # Use plogis(LP) for precision
      lp_Y <- -1 + 0.3 * W + 0.5 * L + beta_YA * a_val + beta_YY * Y_prev
      p_Y  <- plogis(lp_Y)
      means[t + 1] <- mean(p_Y)

      # Draw binary Y for forward simulation (only matters when use_ylag=TRUE;
      # when FALSE, Y_prev is never used so draws are harmless but inert)
      Y <- rbinom(n_mc, 1, p_Y)

      L_prev <- L
      A_prev <- rep(a_val, n_mc)
      Y_prev <- Y
    }

    out[[regime]] <- data.frame(regime = regime, time = 0:2,
                                true_mean = means,
                                stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}


#' Monte Carlo verification — compare two seeds for stability
#'
#' @param n_mc Integer. Monte Carlo sample size (default 1e6).
#' @param scenario "alternative" or "null"
#' @param use_ylag Logical.
mc_verify_true_values <- function(n_mc = 1e6, scenario = "alternative",
                                  use_ylag = FALSE) {
  truth1 <- true_values(scenario, use_ylag, n_mc = n_mc, seed = 42)
  truth2 <- true_values(scenario, use_ylag, n_mc = n_mc, seed = 123)

  ylag_label <- if (use_ylag) "with Y_lag" else "without Y_lag"
  cat(sprintf("\nMonte Carlo verification (n = %s, %s, %s)\n",
              format(n_mc, big.mark = ","), scenario, ylag_label))
  cat(sprintf("%-8s %-6s %-12s %-12s %-10s\n",
              "Regime", "Time", "Seed 42", "Seed 123", "Diff"))
  cat(strrep("-", 50), "\n")

  for (i in seq_len(nrow(truth1))) {
    cat(sprintf("%-8s %-6d %-12.4f %-12.4f %-10.4f\n",
                truth1$regime[i], truth1$time[i],
                truth1$true_mean[i], truth2$true_mean[i],
                truth1$true_mean[i] - truth2$true_mean[i]))
  }
  invisible(truth1)
}


# === SECTION 3: SIMULATION RUNNER ==============================================

#' Run the simulation study
#'
#' @param n Sample size per simulation.
#' @param n_sim Number of simulation replicates.
#' @param scenario "alternative" or "null".
#' @param seed_start First seed.
#' @param use_ylag Logical.
#' @return data.frame with one row per (sim, regime, time).
run_simulation <- function(n = 500, n_sim = 500, scenario = "alternative",
                           seed_start = 1, use_ylag = FALSE) {
  ylag_label <- if (use_ylag) "with Y_lag" else "without Y_lag"
  cat(sprintf("\nRunning %d simulations (n=%d, %s, %s)\n",
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
        outcome_type = "binary",
        regimes = list(always = 1L, never = 0L),
        estimator = "ipw",
        learners = NULL,
        inference = "ic",
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
        se       = if ("se" %in% names(est)) est$se else NA_real_,
        ci_lower = if ("ci_lower" %in% names(est)) est$ci_lower else NA_real_,
        ci_upper = if ("ci_upper" %in% names(est)) est$ci_upper else NA_real_,
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

#' Summarize simulation results: per-regime metrics
#'
#' @param sim_results data.frame from run_simulation().
#' @param truth data.frame from true_values().
#' @return data.frame with bias, rmse, coverage, mean_se, emp_se per (regime, time).
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
        mean_se   = mean(x$se, na.rm = TRUE),
        coverage  = mean(x$ci_lower <= x$true_mean[1] &
                         x$ci_upper >= x$true_mean[1], na.rm = TRUE),
        n_sims    = nrow(x),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(summ) <- NULL
  summ[order(summ$regime, summ$time), ]
}


#' Summarize ATE across regimes
#'
#' Computes ATE = E[Y^always] - E[Y^never] per simulation, then summarizes.
#' ATE standard error is approximated as sqrt(se_always^2 + se_never^2)
#' (conservative: ignores positive correlation from shared data).
#'
#' @param sim_results data.frame from run_simulation().
#' @param truth data.frame from true_values().
#' @return data.frame with one row per time: bias, rmse, coverage, type1_error.
summarize_ate <- function(sim_results, truth) {
  a <- sim_results[sim_results$regime == "always",
                   c("sim", "time", "estimate", "se", "ci_lower", "ci_upper")]
  n <- sim_results[sim_results$regime == "never",
                   c("sim", "time", "estimate", "se", "ci_lower", "ci_upper")]
  names(a)[3:6] <- paste0(names(a)[3:6], "_a")
  names(n)[3:6] <- paste0(names(n)[3:6], "_n")

  ate <- merge(a, n, by = c("sim", "time"))
  ate$ate_est <- ate$estimate_a - ate$estimate_n
  ate$ate_se  <- sqrt(ate$se_a^2 + ate$se_n^2)
  z <- qnorm(0.975)
  ate$ate_lo <- ate$ate_est - z * ate$ate_se
  ate$ate_hi <- ate$ate_est + z * ate$ate_se

  # True ATE
  ta <- truth$true_mean[truth$regime == "always"]
  tn <- truth$true_mean[truth$regime == "never"]
  true_ate <- data.frame(time = 0:2, true_ate = ta - tn)

  ate <- merge(ate, true_ate, by = "time")

  do.call(rbind, lapply(split(ate, ate$time), function(x) {
    data.frame(
      time        = x$time[1],
      true_ate    = x$true_ate[1],
      mean_ate    = mean(x$ate_est),
      bias        = mean(x$ate_est) - x$true_ate[1],
      rmse        = sqrt(mean((x$ate_est - x$true_ate[1])^2)),
      emp_se      = sd(x$ate_est),
      mean_se     = mean(x$ate_se, na.rm = TRUE),
      coverage    = mean(x$ate_lo <= x$true_ate[1] &
                         x$ate_hi >= x$true_ate[1], na.rm = TRUE),
      reject_rate = mean(x$ate_lo > 0 | x$ate_hi < 0, na.rm = TRUE),
      n_sims      = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
}


#' Print a formatted summary table
print_summary <- function(scenario_label, regime_summ, ate_summ) {
  cat(sprintf("\n%s\n%s\n", scenario_label, strrep("=", nchar(scenario_label))))

  cat("\nPer-regime estimates:\n")
  cat(sprintf("  %-8s %-6s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
              "Regime", "Time", "Truth", "Mean", "Bias", "RMSE",
              "EmpSE", "MeanSE", "Covg"))
  cat(sprintf("  %s\n", strrep("-", 72)))
  for (i in seq_len(nrow(regime_summ))) {
    r <- regime_summ[i, ]
    cat(sprintf("  %-8s %-6d %-8.3f %-8.3f %-+8.4f %-8.4f %-8.4f %-8.4f %-8.3f\n",
                r$regime, r$time, r$true_mean, r$mean_est, r$bias,
                r$rmse, r$emp_se, r$mean_se, r$coverage))
  }

  cat("\nATE (always - never):\n")
  cat(sprintf("  %-6s %-8s %-8s %-+8s %-8s %-8s %-8s %-8s %-8s\n",
              "Time", "Truth", "Mean", "Bias", "RMSE",
              "EmpSE", "MeanSE", "Covg", "Reject"))
  cat(sprintf("  %s\n", strrep("-", 72)))
  for (i in seq_len(nrow(ate_summ))) {
    r <- ate_summ[i, ]
    cat(sprintf("  %-6d %-8.3f %-8.3f %-+8.4f %-8.4f %-8.4f %-8.4f %-8.3f %-8.3f\n",
                r$time, r$true_ate, r$mean_ate, r$bias,
                r$rmse, r$emp_se, r$mean_se, r$coverage, r$reject_rate))
  }
}


# === SECTION 5: RUN ============================================================

#' Run the full simulation study (both Y_lag variants)
#'
#' @param n Sample size per dataset.
#' @param n_sim Number of replicates.
#' @param seed_start First seed.
run_all <- function(n = 500, n_sim = 500, seed_start = 1) {
  all_out <- list()

  for (use_ylag in c(FALSE, TRUE)) {
    ylag_tag <- if (use_ylag) "ylag" else "no_ylag"
    ylag_label <- if (use_ylag) "WITH Y_lag" else "WITHOUT Y_lag"

    # --- Alternative scenario ---
    truth_alt <- true_values("alternative", use_ylag = use_ylag)
    sim_alt   <- run_simulation(n, n_sim, "alternative", seed_start,
                                use_ylag = use_ylag)
    regime_alt <- summarize_results(sim_alt, truth_alt)
    ate_alt    <- summarize_ate(sim_alt, truth_alt)

    ate_vals <- truth_alt$true_mean[truth_alt$regime == "always"] -
      truth_alt$true_mean[truth_alt$regime == "never"]
    ate_label <- paste(sprintf("%.2f", ate_vals), collapse = ", ")
    print_summary(
      sprintf("ALTERNATIVE — %s (ATE ~ %s)", ylag_label, ate_label),
      regime_alt, ate_alt
    )

    # --- Null scenario ---
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
    cat(sprintf("Alt: max |bias| = %.4f (target < 0.05)\n",
                max(abs(regime_alt$bias))))
    cat(sprintf("Alt: coverage range = [%.3f, %.3f] (target ~0.95)\n",
                min(regime_alt$coverage), max(regime_alt$coverage)))
    cat(sprintf("Null: type I error range = [%.3f, %.3f] (target ~0.05)\n",
                min(ate_null$reject_rate), max(ate_null$reject_rate)))

    all_out[[ylag_tag]] <- list(
      alternative = list(sim = sim_alt, regime = regime_alt,
                         ate = ate_alt, truth = truth_alt),
      null        = list(sim = sim_null, regime = regime_null,
                         ate = ate_null, truth = truth_null)
    )

    seed_start <- seed_start + 2 * n_sim
  }

  # --- Cross-variant comparison ---
  cat("\n\n")
  cat("==================================================\n")
  cat("  COMPARISON: Y_lag effect on bias and coverage\n")
  cat("==================================================\n")
  cat(sprintf("  %-10s %-6s  %-10s %-10s  %-10s %-10s\n",
              "", "", "No Y_lag", "", "With Y_lag", ""))
  cat(sprintf("  %-10s %-6s  %-10s %-10s  %-10s %-10s\n",
              "Scenario", "Time", "Bias", "Covg", "Bias", "Covg"))
  cat(sprintf("  %s\n", strrep("-", 62)))
  for (scenario in c("alternative", "null")) {
    an <- all_out$no_ylag[[scenario]]$ate
    ay <- all_out$ylag[[scenario]]$ate
    label <- if (scenario == "alternative") "Alt ATE" else "Null ATE"
    for (t in 0:2) {
      bn <- an$bias[an$time == t]
      cn <- an$coverage[an$time == t]
      by <- ay$bias[ay$time == t]
      cy <- ay$coverage[ay$time == t]
      cat(sprintf("  %-10s %-6d  %-+10.4f %-10.3f  %-+10.4f %-10.3f\n",
                  label, t, bn, cn, by, cy))
    }
    cat(sprintf("  %s\n", strrep("-", 62)))
  }

  invisible(all_out)
}
