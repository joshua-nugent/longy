test_that("longy() runs end-to-end with simulated data", {
  d <- simulate_test_data(n = 100, K = 4)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L, never = 0L),
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_equal(length(results), 2)
  expect_true("always" %in% names(results))
  expect_true("never" %in% names(results))
  expect_s3_class(results$always, "longy_result")
  expect_s3_class(results$never, "longy_result")
})

test_that("longy() results match manual pipeline", {
  d <- simulate_test_data(n = 80, K = 3, seed = 99)

  # Manual pipeline
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  lrn <- c("SL.glm", "SL.mean")
  obj <- fit_treatment(obj, regime = "always", learners = lrn, verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", learners = lrn, verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", learners = lrn, verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")
  manual_result <- estimate_ipw(obj, regime = "always")

  # Wrapper (uses same default learners)
  wrapper_results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    verbose = FALSE
  )

  # Estimates should match
  manual_est <- manual_result$estimates$estimate
  wrapper_est <- wrapper_results$always$estimates$estimate

  # Match the times that both have
  shared_times <- intersect(manual_result$estimates$time,
                            wrapper_results$always$estimates$time)
  for (tt in shared_times) {
    m <- manual_result$estimates$estimate[manual_result$estimates$time == tt]
    w <- wrapper_results$always$estimates$estimate[wrapper_results$always$estimates$time == tt]
    expect_equal(m, w, tolerance = 0.02,
                 info = sprintf("time %d mismatch", tt))
  }
})

test_that("longy() works without censoring or observation", {
  d <- simulate_no_censoring(n = 100, K = 3)
  d$Y[is.na(d$Y)] <- 0L

  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = NULL, observation = NULL,
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_true(nrow(results$always$estimates) > 0)
})

test_that("longy() print works", {
  d <- simulate_test_data(n = 50, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    verbose = FALSE
  )

  expect_output(print(results), "longy results")
})

test_that("longy() with truncation", {
  d <- simulate_test_data(n = 80, K = 4)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    truncation = 10,
    verbose = FALSE
  )

  # Check weights were truncated
  w <- results$always$obj$weights$weights_dt$.final_weight
  expect_true(all(w <= 10))
})

test_that("longy() with continuous outcomes", {
  d <- simulate_continuous_outcome(n = 100, K = 4)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    outcome_type = "continuous",
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  est <- results$always$estimates
  expect_true(nrow(est) > 0)
  expect_true(all(is.finite(est$estimate)))
})

test_that("longy() with survival outcomes produces monotone estimates", {
  d <- simulate_survival_outcome(n = 100, K = 5)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    outcome_type = "survival",
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  est <- results$always$estimates
  expect_true(nrow(est) > 0)
  # Estimates should be monotonically non-decreasing
  if (nrow(est) > 1) {
    expect_true(all(diff(est$estimate) >= -1e-10))
  }
  # Estimates should be in [0, 1]
  expect_true(all(est$estimate >= 0 & est$estimate <= 1, na.rm = TRUE))
})

test_that("longy() end-to-end with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 200, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    learners = c("SL.glm", "SL.mean"),
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  est <- results$always$estimates
  expect_true(nrow(est) > 0)
  expect_true(all(est$estimate >= 0 & est$estimate <= 1))

  # SL info should be stored
  sl_info <- results$always$obj$fits$treatment$sl_info
  expect_true(length(sl_info) > 0)
})
