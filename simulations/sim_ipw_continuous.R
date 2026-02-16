# ===========================================================================
# Simulation Study: IPW with Continuous Outcome (3 Time Points)
# ===========================================================================
#
# Verifies longy's IPW estimator under a DGP with time-varying confounding.
# Not a unit test — too slow for every build. Run manually at milestones.
#
# DGP: 3 time points, continuous Y, binary A, absorbing C, no intermittent
# missingness. Past treatment affects future confounders (L), which affect
# both future treatment and outcome. True counterfactual means are analytic.
#
# Structural equations:
#   L(t) = 0.5*W + 0.5*L(t-1) + 0.4*A(t-1) + N(0, 0.25)
#   A(t) ~ Bern(expit(-0.3 + 0.3*W + 0.5*L(t)))
#   C(t) ~ Bern(expit(-3.5 + 0.3*L(t)))
#   Y(t) = 2 + 0.3*W + 0.5*L(t) + 1.0*A(t) + N(0, 1)
#
# True ATE: 1.0 (t=0), 1.2 (t=1), 1.3 (t=2)
# The increasing ATE reflects indirect effects accumulating through L.
#
# Expected results (n=500, n_sim=500):
#   Times 0-1: bias < 0.01, coverage ~95%, type I error ~3-5%
#   Time 2: modest positive bias (~0.03-0.06), coverage ~85-93%.
#     This degradation reflects the known finite-sample behavior of IPW
#     with cumulative weights over 3 time points — weight variability
#     and Hajek estimator bias increase with the number of time steps.
#
# Usage:
#   source("longy/simulations/sim_ipw_continuous.R")
#   run_all()                    # full simulation (~30 sec)
#   run_all(n_sim = 10)          # quick smoke test
#   mc_verify_true_values()      # verify analytic truth via Monte Carlo
# ===========================================================================

library(longy)

# === SECTION 1: DGP ===========================================================

#' Generate one dataset from the simulation DGP
#'
#' @param n Integer. Number of subjects.
#' @param seed Integer. Random seed.
#' @param scenario Character. "alternative" (treatment effect present) or
#'   "null" (no treatment effect, no treatment -> confounder effect).
#' @return data.frame in longy long format (id, time, W, L, A, C, Y).
generate_data <- function(n, seed, scenario = "alternative") {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0

  W <- rnorm(n)
  L_prev <- rep(0, n)       # L(-1) = 0
  A_prev <- rep(0L, n)      # A(-1) = 0
  alive   <- rep(TRUE, n)

  all_rows <- vector("list", 3)

  for (t in 0:2) {
    idx <- which(alive)
    ni  <- length(idx)
    if (ni == 0) break

    # L(t) = 0.5*W + 0.5*L(t-1) + beta_LA*A(t-1) + eps_L, eps_L ~ N(0,0.25)
    L <- numeric(n)
    L[idx] <- 0.5 * W[idx] + 0.5 * L_prev[idx] + beta_LA * A_prev[idx] +
      rnorm(ni, 0, 0.5)

    # A(t) ~ Bernoulli(expit(-0.3 + 0.3*W + 0.5*L(t)))
    # Coefficients chosen to give ~43% treatment rate with bounded weights.
    # Stronger coefficients (e.g. 0.8*L) create near-positivity violations at
    # later time points, causing extreme cumulative weights and IPW instability.
    A <- integer(n)
    A[idx] <- rbinom(ni, 1, plogis(-0.3 + 0.3 * W[idx] + 0.5 * L[idx]))

    # C(t) ~ Bernoulli(expit(-3.5 + 0.3*L(t)))  [~3-5% per time]
    C <- integer(n)
    C[idx] <- rbinom(ni, 1, plogis(-3.5 + 0.3 * L[idx]))

    # Y(t) = 2 + 0.3*W + 0.5*L(t) + beta_YA*A(t) + eps_Y, eps_Y ~ N(0,1)
    # Only observed when uncensored
    Y <- rep(NA_real_, n)
    unc <- idx[C[idx] == 0L]
    if (length(unc) > 0) {
      Y[unc] <- 2 + 0.3 * W[unc] + 0.5 * L[unc] +
        beta_YA * A[unc] + rnorm(length(unc))
    }

    all_rows[[t + 1]] <- data.frame(
      id = idx, time = t,
      W = W[idx], L = L[idx],
      A = A[idx], C = C[idx], Y = Y[idx]
    )

    # Update state for next time point
    L_prev[idx] <- L[idx]
    A_prev[idx] <- A[idx]
    alive[idx[C[idx] == 1L]] <- FALSE
  }

  do.call(rbind, all_rows)
}


# === SECTION 2: TRUE VALUES ====================================================

#' Analytic true counterfactual means
#'
#' All structural equations are linear with E[W]=0, so the marginal means
#' simplify to functions of the coefficients and regime values only.
#'
#' @param scenario "alternative" or "null"
#' @return data.frame with columns: regime, time, true_mean
true_values <- function(scenario = "alternative") {
  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0

  out <- list()
  for (regime in c("always", "never")) {
    a_val <- if (regime == "always") 1 else 0
    EL <- numeric(3)
    # E[L(0)] = 0.5*E[W] + 0.5*0 + beta_LA*0 = 0  (A(-1)=0)
    EL[1] <- 0
    # E[L(t)] = 0.5*E[L(t-1)] + beta_LA * a_val  (E[W]=0 drops)
    for (k in 2:3) EL[k] <- 0.5 * EL[k - 1] + beta_LA * a_val
    # E[Y(t)] = 2 + 0.3*E[W] + 0.5*E[L(t)] + beta_YA*a_val = 2 + 0.5*E[L(t)] + beta_YA*a_val
    EY <- 2 + 0.5 * EL + beta_YA * a_val
    out[[regime]] <- data.frame(regime = regime, time = 0:2, true_mean = EY)
  }
  do.call(rbind, out)
}


#' Monte Carlo verification of analytic true values
#'
#' Generates a huge sample under each intervention (ignoring censoring),
#' computes sample means. Should match analytic values to ~2 decimal places.
#'
#' @param n_mc Integer. Monte Carlo sample size (default 1e6).
#' @param scenario "alternative" or "null"
#' @param seed Random seed.
mc_verify_true_values <- function(n_mc = 1e6, scenario = "alternative",
                                  seed = 99) {
  set.seed(seed)

  beta_LA <- if (scenario == "null") 0 else 0.4
  beta_YA <- if (scenario == "null") 0 else 1.0

  truth <- true_values(scenario)

  cat(sprintf("\nMonte Carlo verification (n = %s, scenario = %s)\n",
              format(n_mc, big.mark = ","), scenario))
  cat(sprintf("%-8s %-6s %-12s %-12s %-10s\n",
              "Regime", "Time", "MC Mean", "Analytic", "Diff"))
  cat(strrep("-", 50), "\n")

  W <- rnorm(n_mc)

  for (regime in c("always", "never")) {
    a_val <- if (regime == "always") 1 else 0
    L_prev <- rep(0, n_mc)
    A_prev <- rep(0, n_mc)  # A(-1) = 0

    for (t in 0:2) {
      L <- 0.5 * W + 0.5 * L_prev + beta_LA * A_prev + rnorm(n_mc, 0, 0.5)
      Y <- 2 + 0.3 * W + 0.5 * L + beta_YA * a_val + rnorm(n_mc)

      mc_mean  <- mean(Y)
      analytic <- truth$true_mean[truth$regime == regime & truth$time == t]
      cat(sprintf("%-8s %-6d %-12.4f %-12.4f %-10.4f\n",
                  regime, t, mc_mean, analytic, mc_mean - analytic))

      L_prev <- L
      A_prev <- rep(a_val, n_mc)
    }
  }
  invisible(truth)
}


# === SECTION 3: SIMULATION RUNNER ==============================================

#' Run the simulation study
#'
#' @param n Sample size per simulation.
#' @param n_sim Number of simulation replicates.
#' @param scenario "alternative" or "null".
#' @param seed_start First seed (seeds go from seed_start to seed_start+n_sim-1).
#' @return data.frame with one row per (sim, regime, time), columns:
#'   sim, regime, time, estimate, se, ci_lower, ci_upper.
run_simulation <- function(n = 500, n_sim = 500, scenario = "alternative",
                           seed_start = 1) {
  cat(sprintf("\nRunning %d simulations (n=%d, scenario=%s)\n",
              n_sim, n, scenario))
  t0 <- proc.time()

  all_results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    d <- generate_data(n = n, seed = seed_start + i - 1, scenario = scenario)

    res <- tryCatch(
      suppressWarnings(longy(
        data = d,
        id = "id", time = "time", outcome = "Y",
        treatment = "A", censoring = "C",
        baseline = "W", timevarying = "L",
        outcome_type = "continuous",
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

#' Run the full simulation study
#'
#' @param n Sample size per dataset.
#' @param n_sim Number of replicates.
#' @param seed_start First seed.
run_all <- function(n = 500, n_sim = 500, seed_start = 1) {
  # --- Alternative scenario (treatment effect present) ---
  truth_alt <- true_values("alternative")
  sim_alt   <- run_simulation(n, n_sim, "alternative", seed_start)
  regime_alt <- summarize_results(sim_alt, truth_alt)
  ate_alt    <- summarize_ate(sim_alt, truth_alt)
  print_summary("ALTERNATIVE SCENARIO (ATE = 1.0, 1.2, 1.3)", regime_alt, ate_alt)

  # --- Null scenario (no treatment effect) ---
  truth_null <- true_values("null")
  sim_null   <- run_simulation(n, n_sim, "null", seed_start + n_sim)
  regime_null <- summarize_results(sim_null, truth_null)
  ate_null    <- summarize_ate(sim_null, truth_null)
  print_summary("NULL SCENARIO (ATE = 0 at all times)", regime_null, ate_null)

  cat("\n--- Diagnostic checks ---\n")
  cat(sprintf("Alt: max |bias| = %.4f (target < 0.05)\n",
              max(abs(regime_alt$bias))))
  cat(sprintf("Alt: coverage range = [%.3f, %.3f] (target ~0.95)\n",
              min(regime_alt$coverage), max(regime_alt$coverage)))
  cat(sprintf("Null: type I error range = [%.3f, %.3f] (target ~0.05)\n",
              min(ate_null$reject_rate), max(ate_null$reject_rate)))

  invisible(list(
    alternative = list(sim = sim_alt, regime = regime_alt, ate = ate_alt, truth = truth_alt),
    null        = list(sim = sim_null, regime = regime_null, ate = ate_null, truth = truth_null)
  ))
}
