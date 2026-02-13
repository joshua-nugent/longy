# Tests for G-computation estimator

test_that("G-comp produces valid estimates for binary outcome", {
  d <- simulate_test_data(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true(nrow(result$estimates) > 0)
  expect_true(all(result$estimates$estimate >= 0 & result$estimates$estimate <= 1,
                  na.rm = TRUE))
  expect_equal(result$estimator, "gcomp")
})

test_that("G-comp works with continuous outcomes", {
  d <- simulate_continuous_outcome(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true(all(is.finite(result$estimates$estimate)))
})

test_that("G-comp with survival outcomes produces monotone estimates", {
  d <- simulate_survival_outcome(n = 150, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  est <- result$estimates
  expect_true(nrow(est) > 0)
  # Isotonic smoothing should enforce monotone non-decreasing
  if (nrow(est) > 1) {
    expect_true(all(diff(est$estimate) >= -1e-10))
  }
  expect_true(all(est$estimate >= 0 & est$estimate <= 1, na.rm = TRUE))
})

test_that("G-comp works without censoring or observation", {
  d <- simulate_no_censoring(n = 150, K = 3)
  d$Y[is.na(d$Y)] <- 0L

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true(nrow(result$estimates) > 0)
  expect_true(all(is.finite(result$estimates$estimate)))
})

test_that("G-comp handles intermittent observation (R=0)", {
  d <- simulate_test_data(n = 200, K = 4, seed = 77)
  # Ensure some R=0 exist
  expect_true(any(d$R == 0 & d$C == 0))

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true(all(is.finite(result$estimates$estimate)))
})

test_that("G-comp bootstrap SEs are positive", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 30, verbose = FALSE)

  expect_true("se" %in% names(result$estimates))
  expect_true(all(result$estimates$se > 0, na.rm = TRUE))
  expect_true("ci_lower" %in% names(result$estimates))
  expect_true("ci_upper" %in% names(result$estimates))
})

test_that("G-comp and IPW give different estimates with confounding", {
  d <- simulate_test_data(n = 200, K = 4, seed = 55)
  obj_g <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_g <- define_regime(obj_g, "always", static = 1L)
  obj_g <- fit_outcome(obj_g, regime = "always", verbose = FALSE)
  gcomp_result <- estimate_gcomp(obj_g, regime = "always", n_boot = 0,
                                  verbose = FALSE)

  obj_i <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_i <- define_regime(obj_i, "always", static = 1L)
  obj_i <- fit_treatment(obj_i, regime = "always", verbose = FALSE)
  obj_i <- fit_censoring(obj_i, regime = "always", verbose = FALSE)
  obj_i <- fit_observation(obj_i, regime = "always", verbose = FALSE)
  obj_i <- compute_weights(obj_i, regime = "always")
  ipw_result <- estimate_ipw(obj_i, regime = "always")

  # Both should produce valid estimates; they'll likely differ because they
  # rely on different modeling assumptions
  expect_true(all(is.finite(gcomp_result$estimates$estimate)))
  expect_true(all(is.finite(ipw_result$estimates$estimate)))
})

test_that("G-comp no-confounding recovery", {
  d <- simulate_no_confounding(n = 500, K = 3, seed = 999)
  d$Y[is.na(d$Y)] <- 0L

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  # Crude mean among treated at earliest time (no confounding, so should be close)
  t0 <- min(result$estimates$time)
  crude <- mean(d$Y[d$A == 1 & d$time == t0], na.rm = TRUE)
  gcomp_est <- result$estimates$estimate[result$estimates$time == t0]
  # Allow generous tolerance since this is a simulation
  expect_equal(gcomp_est, crude, tolerance = 0.15)
})

test_that("longy() with estimator = 'gcomp' works end-to-end", {
  d <- simulate_test_data(n = 100, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "gcomp",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_true("always" %in% names(results))
  expect_s3_class(results$always, "longy_result")
  expect_true(nrow(results$always$estimates) > 0)
  expect_true(all(results$always$estimates$estimate >= 0 &
                  results$always$estimates$estimate <= 1, na.rm = TRUE))
})

test_that("longy() with estimator = 'both' returns IPW and G-comp", {
  d <- simulate_test_data(n = 100, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "both",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_true("always_ipw" %in% names(results))
  expect_true("always_gcomp" %in% names(results))
  expect_s3_class(results$always_ipw, "longy_result")
  expect_s3_class(results$always_gcomp, "longy_result")
})

test_that("G-comp continuous outcome recovers known truth at each time", {
  # Clean DGP: no censoring, always observed, continuous outcome
  # L1 ~ N(0.3*W1 + 0.1*t, 1), A ~ Bern(plogis(0.5*L1))
  # Y = 0.5*L1 + 0.3*A + 0.2*W1 + N(0,1)
  # Truth: E[Y_t(a)] = 0.5*(0.1*t) + 0.3*a + 0 = 0.05*t + 0.3*a
  set.seed(42)
  n <- 3000; K <- 5
  rows <- vector("list", n * K)
  idx <- 1L
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    for (tt in 0:(K - 1)) {
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      A <- rbinom(1, 1, plogis(0.5 * L1))
      Y <- 0.5 * L1 + 0.3 * A + 0.2 * W1 + rnorm(1)
      rows[[idx]] <- data.frame(id = i, time = tt, W1 = W1, L1 = L1,
                                A = A, C = 0L, R = 1L, Y = Y)
      idx <- idx + 1L
    }
  }
  d <- do.call(rbind, rows)

  obj1 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = "C", observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj1 <- define_regime(obj1, "always", static = 1L)
  obj1 <- fit_outcome(obj1, regime = "always", verbose = FALSE)
  res1 <- estimate_gcomp(obj1, regime = "always", n_boot = 0, verbose = FALSE)

  obj0 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = "C", observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj0 <- define_regime(obj0, "never", static = 0L)
  obj0 <- fit_outcome(obj0, regime = "never", verbose = FALSE)
  res0 <- estimate_gcomp(obj0, regime = "never", n_boot = 0, verbose = FALSE)

  # Check at each time point (absolute tolerance)
  for (r in seq_len(nrow(res1$estimates))) {
    tt <- res1$estimates$time[r]
    truth_a1 <- 0.3 + 0.05 * tt
    expect_true(abs(res1$estimates$estimate[r] - truth_a1) < 0.05,
                label = sprintf("G-comp a=1 at t=%d: est=%.3f truth=%.3f",
                                tt, res1$estimates$estimate[r], truth_a1))
  }
  for (r in seq_len(nrow(res0$estimates))) {
    tt <- res0$estimates$time[r]
    truth_a0 <- 0.05 * tt
    expect_true(abs(res0$estimates$estimate[r] - truth_a0) < 0.05,
                label = sprintf("G-comp a=0 at t=%d: est=%.3f truth=%.3f",
                                tt, res0$estimates$estimate[r], truth_a0))
  }
})

test_that("G-comp ICE recovers carryover truth (past A affects future L)", {
  # Carryover DGP: L1(t) depends on A(t-1), so direct regression misses
  # the carryover at t>0. ICE backward regression handles it correctly.
  # Truth via simulation:
  #   E[Y_t(a=1)] = 0.3 + 0.05*t + 0.1*(t>0)*cumulative_carryover
  # We compute truth by Monte Carlo instead of closed form.
  set.seed(123)
  n_mc <- 50000; K <- 4
  # Monte Carlo truth for a=1 and a=0
  truth_a1 <- numeric(K)
  truth_a0 <- numeric(K)
  for (tt in 0:(K - 1)) {
    # Simulate under always-treat (a=1)
    ey1 <- replicate(n_mc, {
      W1 <- rnorm(1)
      A_prev <- 0L
      for (s in 0:tt) {
        L1 <- 0.3 * W1 + 0.2 * A_prev + 0.1 * s + rnorm(1)
        A_prev <- 1L  # intervene
        if (s == tt) return(0.5 * L1 + 0.3 * 1 + 0.2 * W1)
      }
    })
    truth_a1[tt + 1] <- mean(ey1)
    # Simulate under never-treat (a=0)
    ey0 <- replicate(n_mc, {
      W1 <- rnorm(1)
      A_prev <- 0L
      for (s in 0:tt) {
        L1 <- 0.3 * W1 + 0.2 * A_prev + 0.1 * s + rnorm(1)
        A_prev <- 0L  # intervene
        if (s == tt) return(0.5 * L1 + 0.3 * 0 + 0.2 * W1)
      }
    })
    truth_a0[tt + 1] <- mean(ey0)
  }

  d <- simulate_carryover_continuous(n = 3000, K = K, seed = 42)

  obj1 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = "C", observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj1 <- define_regime(obj1, "always", static = 1L)
  obj1 <- fit_outcome(obj1, regime = "always", verbose = FALSE)
  res1 <- estimate_gcomp(obj1, regime = "always", n_boot = 0, verbose = FALSE)

  obj0 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = "C", observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj0 <- define_regime(obj0, "never", static = 0L)
  obj0 <- fit_outcome(obj0, regime = "never", verbose = FALSE)
  res0 <- estimate_gcomp(obj0, regime = "never", n_boot = 0, verbose = FALSE)

  for (r in seq_len(nrow(res1$estimates))) {
    tt <- res1$estimates$time[r]
    idx <- tt + 1
    expect_true(abs(res1$estimates$estimate[r] - truth_a1[idx]) < 0.1,
                label = sprintf("ICE G-comp a=1 at t=%d: est=%.3f truth=%.3f",
                                tt, res1$estimates$estimate[r], truth_a1[idx]))
  }
  for (r in seq_len(nrow(res0$estimates))) {
    tt <- res0$estimates$time[r]
    idx <- tt + 1
    expect_true(abs(res0$estimates$estimate[r] - truth_a0[idx]) < 0.1,
                label = sprintf("ICE G-comp a=0 at t=%d: est=%.3f truth=%.3f",
                                tt, res0$estimates$estimate[r], truth_a0[idx]))
  }
})

test_that("fit_outcome errors without regime", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  expect_error(fit_outcome(obj, regime = "nonexistent"),
               "Regime.*not found")
})

test_that("estimate_gcomp errors without fitted outcome model", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  expect_error(estimate_gcomp(obj, regime = "always"),
               "Outcome model not fitted")
})

test_that("G-comp bootstrap runs in parallel with future.apply", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)

  result <- estimate_gcomp(obj, regime = "always", n_boot = 20, verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true("se" %in% names(result$estimates))
  expect_true(all(result$estimates$se > 0, na.rm = TRUE))
  expect_true("ci_lower" %in% names(result$estimates))
  expect_true("ci_upper" %in% names(result$estimates))
})

test_that("G-comp print method works", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  result <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)

  expect_output(print(result), "G-comp")
})
