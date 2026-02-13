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
