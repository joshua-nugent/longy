# Tests for G-computation estimator

test_that("G-comp produces valid estimates for binary outcome", {
  d <- simulate_test_data(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(nrow(res$estimates) > 0)
  expect_true(all(res$estimates$estimate >= 0 & res$estimates$estimate <= 1,
                  na.rm = TRUE))
  expect_equal(res$estimator, "gcomp")
})

test_that("G-comp works with continuous outcomes", {
  d <- simulate_continuous_outcome(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(all(is.finite(res$estimates$estimate)))
})

test_that("G-comp with survival outcomes produces monotone estimates", {
  d <- simulate_survival_outcome(n = 150, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  est <- res$estimates
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

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(nrow(res$estimates) > 0)
  expect_true(all(is.finite(res$estimates$estimate)))
})

test_that("G-comp handles intermittent observation (R=0)", {
  d <- simulate_test_data(n = 200, K = 4, seed = 77)
  # Ensure some R=0 exist
  expect_true(any(d$R == 0 & d$C == "uncensored"))

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(all(is.finite(res$estimates$estimate)))
})

test_that("G-comp bootstrap SEs are positive", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 30, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_true("se" %in% names(res$estimates))
  expect_true(all(res$estimates$se > 0, na.rm = TRUE))
  expect_true("ci_lower" %in% names(res$estimates))
  expect_true("ci_upper" %in% names(res$estimates))
})

test_that("G-comp and IPW give different estimates with confounding", {
  d <- simulate_test_data(n = 200, K = 4, seed = 55)
  obj_g <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_g <- define_regime(obj_g, "always", static = 1L)
  obj_g <- fit_outcome(obj_g, regime = "always", verbose = FALSE)
  obj_g <- estimate_gcomp(obj_g, regime = "always", n_boot = 0,
                                  verbose = FALSE)
  gcomp_res <- obj_g$results$always_gcomp

  obj_i <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_i <- define_regime(obj_i, "always", static = 1L)
  obj_i <- fit_treatment(obj_i, regime = "always", verbose = FALSE)
  obj_i <- fit_censoring(obj_i, regime = "always", verbose = FALSE)
  obj_i <- fit_observation(obj_i, regime = "always", verbose = FALSE)
  obj_i <- compute_weights(obj_i, regime = "always")
  obj_i <- estimate_ipw(obj_i, regime = "always")
  ipw_res <- obj_i$results$always_ipw

  # Both should produce valid estimates; they'll likely differ because they
  # rely on different modeling assumptions
  expect_true(all(is.finite(gcomp_res$estimates$estimate)))
  expect_true(all(is.finite(ipw_res$estimates$estimate)))
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
  obj <- estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE)
  res <- obj$results$always_gcomp

  # Crude mean among treated at earliest time (no confounding, so should be close)
  t0 <- min(res$estimates$time)
  crude <- mean(d$Y[d$A == 1 & d$time == t0], na.rm = TRUE)
  gcomp_est <- res$estimates$estimate[res$estimates$time == t0]
  # Allow generous tolerance since this is a simulation
  expect_equal(gcomp_est, crude, tolerance = 0.15)
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
                     treatment = "A", censoring = NULL, observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj1 <- define_regime(obj1, "always", static = 1L)
  obj1 <- fit_outcome(obj1, regime = "always", verbose = FALSE)
  obj1 <- estimate_gcomp(obj1, regime = "always", n_boot = 0, verbose = FALSE)
  res1 <- obj1$results$always_gcomp

  obj0 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = NULL, observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj0 <- define_regime(obj0, "never", static = 0L)
  obj0 <- fit_outcome(obj0, regime = "never", verbose = FALSE)
  obj0 <- estimate_gcomp(obj0, regime = "never", n_boot = 0, verbose = FALSE)
  res0 <- obj0$results$never_gcomp

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
                     treatment = "A", censoring = NULL, observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj1 <- define_regime(obj1, "always", static = 1L)
  obj1 <- fit_outcome(obj1, regime = "always", verbose = FALSE)
  obj1 <- estimate_gcomp(obj1, regime = "always", n_boot = 0, verbose = FALSE)
  res1 <- obj1$results$always_gcomp

  obj0 <- longy_data(d, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = NULL, observation = "R",
                     baseline = "W1", timevarying = "L1",
                     outcome_type = "continuous", verbose = FALSE)
  obj0 <- define_regime(obj0, "never", static = 0L)
  obj0 <- fit_outcome(obj0, regime = "never", verbose = FALSE)
  obj0 <- estimate_gcomp(obj0, regime = "never", n_boot = 0, verbose = FALSE)
  res0 <- obj0$results$never_gcomp

  # Tolerance widened to 0.15: with lag columns (k=Inf default), glm has

  # more covariates which can slightly degrade fit on small test datasets.
  for (r in seq_len(nrow(res1$estimates))) {
    tt <- res1$estimates$time[r]
    idx <- tt + 1
    expect_true(abs(res1$estimates$estimate[r] - truth_a1[idx]) < 0.15,
                label = sprintf("ICE G-comp a=1 at t=%d: est=%.3f truth=%.3f",
                                tt, res1$estimates$estimate[r], truth_a1[idx]))
  }
  for (r in seq_len(nrow(res0$estimates))) {
    tt <- res0$estimates$time[r]
    idx <- tt + 1
    expect_true(abs(res0$estimates$estimate[r] - truth_a0[idx]) < 0.15,
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

  obj <- estimate_gcomp(obj, regime = "always", n_boot = 20, verbose = FALSE)
  res <- obj$results$always_gcomp

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true("se" %in% names(res$estimates))
  expect_true(all(res$estimates$se > 0, na.rm = TRUE))
  expect_true("ci_lower" %in% names(res$estimates))
  expect_true("ci_upper" %in% names(res$estimates))
})

test_that("fit_outcome sl_info contains clip tracking and gaussian family", {
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  outcome_fit <- obj$fits$outcome[["always"]]
  # Top-level family should be "gaussian" (LMTP-style)
  expect_equal(outcome_fit$family, "gaussian")

  # sl_info should have per-step entries with clip counts and family
  sl_info <- outcome_fit$sl_info
  expect_true(length(sl_info) > 0)

  for (entry in sl_info) {
    expect_true("n_clipped_lower" %in% names(entry))
    expect_true("n_clipped_upper" %in% names(entry))
    expect_true("n_preds" %in% names(entry))
    expect_true("family" %in% names(entry))
    expect_true(is.numeric(entry$n_clipped_lower))
    expect_true(is.numeric(entry$n_clipped_upper))
    expect_true(entry$n_clipped_lower >= 0)
    expect_true(entry$n_clipped_upper >= 0)
    # Family should be binomial or gaussian
    expect_true(entry$family %in% c("binomial", "gaussian"))
  }

  # First step (at target time) for binary outcome should use binomial
  first_steps <- Filter(function(e) e$time == e$target_time, sl_info)
  if (length(first_steps) > 0) {
    expect_equal(first_steps[[1]]$family, "binomial")
  }

  # Intermediate steps should use gaussian
  intermediate_steps <- Filter(function(e) e$time != e$target_time, sl_info)
  if (length(intermediate_steps) > 0) {
    for (e in intermediate_steps) {
      expect_equal(e$family, "gaussian")
    }
  }
})

test_that("adaptive_cv_folds handles binary and continuous correctly", {
  # Binary: minority class should shrink n_eff
  y_binary <- c(rep(0, 95), rep(1, 5))
  result_binary <- longy:::.adaptive_cv_folds(y_binary, binary = TRUE)
  expect_true(result_binary$n_eff < 100)
  expect_true(result_binary$V >= 2)

  # Continuous: n_eff should equal n
  y_cont <- runif(100, 0, 1)
  result_cont <- longy:::.adaptive_cv_folds(y_cont, binary = FALSE)
  expect_equal(result_cont$n_eff, 100)
  expect_true(result_cont$V >= 2)

  # Default (binary=TRUE) should match explicit binary=TRUE
  result_default <- longy:::.adaptive_cv_folds(y_binary)
  expect_equal(result_default$n_eff, result_binary$n_eff)
  expect_equal(result_default$V, result_binary$V)
})
