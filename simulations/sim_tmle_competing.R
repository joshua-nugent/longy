# Simulation: TMLE with competing risks
# Estimates cause-specific cumulative incidence of the primary event
# under always-treat and never-treat regimes.

library(longy)

# --- DGP: Competing risks with confounding and censoring ---
simulate_cr_dgp <- function(n, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    alive <- TRUE
    primary_occurred <- FALSE
    competing_occurred <- FALSE
    for (tt in 0:(K - 1)) {
      if (!alive) break
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.3 * L1))
      p_a <- plogis(-0.2 + 0.4 * L1 + 0.3 * W1 - 0.2 * L2)
      A <- rbinom(1, 1, p_a)
      p_c <- plogis(-3 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
      C <- rbinom(1, 1, p_c)
      if (C == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = 1L, R = 1L,
          Y = NA_integer_, D = NA_integer_
        )
        alive <- FALSE
        next
      }
      if (primary_occurred) {
        Y <- 1L; D <- 0L
      } else if (competing_occurred) {
        Y <- 0L; D <- 1L
      } else {
        p_primary <- plogis(-3 + 0.3 * L1 + 0.2 * A + 0.15 * tt)
        p_competing <- plogis(-3.5 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
        u <- runif(1)
        if (u < p_primary) {
          Y <- 1L; D <- 0L; primary_occurred <- TRUE
        } else if (u < p_primary + p_competing) {
          Y <- 0L; D <- 1L; competing_occurred <- TRUE
        } else {
          Y <- 0L; D <- 0L
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = 0L, R = 1L,
        Y = Y, D = D
      )
    }
  }
  do.call(rbind, rows)
}

# --- Monte Carlo true values ---
# Simulate large dataset under intervention (always treat) with no censoring
mc_true <- function(n_mc = 50000, K = 5, regime_a = 1L, seed = 999) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n_mc)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    primary_occurred <- FALSE
    competing_occurred <- FALSE
    for (tt in 0:(K - 1)) {
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.3 * L1))
      A <- regime_a  # intervened
      if (primary_occurred) {
        Y <- 1L; D <- 0L
      } else if (competing_occurred) {
        Y <- 0L; D <- 1L
      } else {
        p_primary <- plogis(-3 + 0.3 * L1 + 0.2 * A + 0.15 * tt)
        p_competing <- plogis(-3.5 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
        u <- runif(1)
        if (u < p_primary) {
          Y <- 1L; D <- 0L; primary_occurred <- TRUE
        } else if (u < p_primary + p_competing) {
          Y <- 0L; D <- 1L; competing_occurred <- TRUE
        } else {
          Y <- 0L; D <- 0L
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, Y = Y, D = D
      )
    }
  }
  mc_dt <- do.call(rbind, rows)
  # CIF at each time = P(Y=1 at time t)
  tapply(mc_dt$Y, mc_dt$time, mean)
}

cat("Computing Monte Carlo true values...\n")
true_always <- mc_true(regime_a = 1L, seed = 999)
true_never  <- mc_true(regime_a = 0L, seed = 1000)
cat("True CIF (always):", round(true_always, 4), "\n")
cat("True CIF (never):", round(true_never, 4), "\n")

# --- Simulation ---
n_sim <- 200
n_obs <- 500
K <- 5
results <- list()

cat(sprintf("\nRunning %d simulations (n=%d, K=%d)...\n", n_sim, n_obs, K))

for (sim in seq_len(n_sim)) {
  if (sim %% 10 == 0) cat(sprintf("  sim %d/%d\n", sim, n_sim))

  dat <- simulate_cr_dgp(n = n_obs, K = K, seed = sim)

  res <- tryCatch({
    longy(
      data = dat, id = "id", time = "time", outcome = "Y",
      treatment = "A", censoring = "C",
      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
      outcome_type = "survival", competing = "D",
      regimes = list(always = 1L, never = 0L),
      estimator = "tmle", learners = NULL,
      inference = "eif", verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  sim %d failed: %s\n", sim, e$message))
    NULL
  })

  if (!is.null(res)) {
    est_always <- res$always_tmle$estimates
    est_never  <- res$never_tmle$estimates
    results[[sim]] <- list(
      always = est_always$estimate,
      never = est_never$estimate,
      always_se = if ("se" %in% names(est_always)) est_always$se else rep(NA, nrow(est_always)),
      never_se = if ("se" %in% names(est_never)) est_never$se else rep(NA, nrow(est_never)),
      always_ci_lo = if ("ci_lower" %in% names(est_always)) est_always$ci_lower else rep(NA, nrow(est_always)),
      always_ci_hi = if ("ci_upper" %in% names(est_always)) est_always$ci_upper else rep(NA, nrow(est_always)),
      never_ci_lo = if ("ci_lower" %in% names(est_never)) est_never$ci_lower else rep(NA, nrow(est_never)),
      never_ci_hi = if ("ci_upper" %in% names(est_never)) est_never$ci_upper else rep(NA, nrow(est_never))
    )
  }
}

# --- Summarize ---
valid <- results[!vapply(results, is.null, logical(1))]
n_valid <- length(valid)
cat(sprintf("\n%d/%d simulations succeeded\n", n_valid, n_sim))

if (n_valid > 0) {
  always_mat <- do.call(rbind, lapply(valid, `[[`, "always"))
  never_mat  <- do.call(rbind, lapply(valid, `[[`, "never"))
  always_se_mat <- do.call(rbind, lapply(valid, `[[`, "always_se"))
  never_se_mat  <- do.call(rbind, lapply(valid, `[[`, "never_se"))
  always_lo_mat <- do.call(rbind, lapply(valid, `[[`, "always_ci_lo"))
  always_hi_mat <- do.call(rbind, lapply(valid, `[[`, "always_ci_hi"))
  never_lo_mat  <- do.call(rbind, lapply(valid, `[[`, "never_ci_lo"))
  never_hi_mat  <- do.call(rbind, lapply(valid, `[[`, "never_ci_hi"))

  cat("\n=== Always Treat ===\n")
  for (k in seq_along(true_always)) {
    est_k <- always_mat[, k]
    se_k <- always_se_mat[, k]
    lo_k <- always_lo_mat[, k]
    hi_k <- always_hi_mat[, k]
    bias <- mean(est_k, na.rm = TRUE) - true_always[k]
    emp_se <- sd(est_k, na.rm = TRUE)
    mean_se <- mean(se_k, na.rm = TRUE)
    coverage <- mean(lo_k <= true_always[k] & hi_k >= true_always[k], na.rm = TRUE)
    cat(sprintf("  t=%d: true=%.4f, mean_est=%.4f, bias=%.4f, emp_SE=%.4f, mean_SE=%.4f, coverage=%.2f\n",
                k - 1, true_always[k], mean(est_k, na.rm = TRUE), bias, emp_se, mean_se, coverage))
  }

  cat("\n=== Never Treat ===\n")
  for (k in seq_along(true_never)) {
    est_k <- never_mat[, k]
    se_k <- never_se_mat[, k]
    lo_k <- never_lo_mat[, k]
    hi_k <- never_hi_mat[, k]
    bias <- mean(est_k, na.rm = TRUE) - true_never[k]
    emp_se <- sd(est_k, na.rm = TRUE)
    mean_se <- mean(se_k, na.rm = TRUE)
    coverage <- mean(lo_k <= true_never[k] & hi_k >= true_never[k], na.rm = TRUE)
    cat(sprintf("  t=%d: true=%.4f, mean_est=%.4f, bias=%.4f, emp_SE=%.4f, mean_SE=%.4f, coverage=%.2f\n",
                k - 1, true_never[k], mean(est_k, na.rm = TRUE), bias, emp_se, mean_se, coverage))
  }
}
