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

  expect_s3_class(results, "longy_data")
  expect_equal(length(results$results), 2)
  expect_true("always_ipw" %in% names(results$results))
  expect_true("never_ipw" %in% names(results$results))
  expect_s3_class(results$results$always_ipw, "longy_result")
  expect_s3_class(results$results$never_ipw, "longy_result")
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
  obj <- estimate_ipw(obj, regime = "always")
  manual_result <- obj$results$always_ipw

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
  wrapper_est <- wrapper_results$results$always_ipw$estimates$estimate

  # Match the times that both have
  shared_times <- intersect(manual_result$estimates$time,
                            wrapper_results$results$always_ipw$estimates$time)
  for (tt in shared_times) {
    m <- manual_result$estimates$estimate[manual_result$estimates$time == tt]
    w <- wrapper_results$results$always_ipw$estimates$estimate[wrapper_results$results$always_ipw$estimates$time == tt]
    expect_equal(m, w, tolerance = 0.02,
                 info = sprintf("time %d mismatch", tt))
  }
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

  expect_output(print(results), "longy")
})

test_that("longy() with estimator = 'gcomp'", {
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

  expect_s3_class(results, "longy_data")
  expect_true("always_gcomp" %in% names(results$results))
  expect_s3_class(results$results$always_gcomp, "longy_result")
  expect_true(nrow(results$results$always_gcomp$estimates) > 0)
  expect_true(all(results$results$always_gcomp$estimates$estimate >= 0 &
                  results$results$always_gcomp$estimates$estimate <= 1, na.rm = TRUE))
})

test_that("longy() with estimator = 'tmle'", {
  # Binary outcome
  d <- simulate_test_data(n = 150, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "tmle",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_data")
  expect_s3_class(results$results$always_tmle, "longy_result")
  expect_equal(results$results$always_tmle$estimator, "tmle")
  expect_true(nrow(results$results$always_tmle$estimates) > 0)
  expect_true(all(results$results$always_tmle$estimates$estimate >= 0 &
                  results$results$always_tmle$estimates$estimate <= 1, na.rm = TRUE))
  expect_true("se" %in% names(results$results$always_tmle$estimates))

  # Continuous outcome
  d_cont <- simulate_continuous_outcome(n = 200, K = 3)
  results_cont <- longy(
    data = d_cont,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "tmle",
    outcome_type = "continuous",
    n_boot = 0,
    verbose = FALSE
  )
  expect_equal(results_cont$results$always_tmle$estimator, "tmle")
  expect_true(all(is.finite(results_cont$results$always_tmle$estimates$estimate)))

  # Never-treat regime
  results_never <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(never = 0L),
    estimator = "tmle",
    n_boot = 0,
    verbose = FALSE
  )
  expect_true("never_tmle" %in% names(results_never$results))
  expect_equal(results_never$results$never_tmle$estimator, "tmle")
  expect_true(nrow(results_never$results$never_tmle$estimates) > 0)
})

test_that("longy() with estimator = 'all'", {
  d <- simulate_test_data(n = 150, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L, never = 0L),
    estimator = "all",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_data")
  expect_equal(length(results$results), 6)  # 2 regimes x 3 estimators
  expect_true(all(c("always_ipw", "always_gcomp", "always_tmle",
                     "never_ipw", "never_gcomp", "never_tmle") %in%
                  names(results$results)))
  expect_s3_class(results$results$always_ipw, "longy_result")
  expect_s3_class(results$results$always_gcomp, "longy_result")
  expect_s3_class(results$results$always_tmle, "longy_result")
  expect_equal(results$results$always_tmle$estimator, "tmle")
  expect_equal(results$results$always_gcomp$estimator, "gcomp")
  expect_true("se" %in% names(results$results$always_tmle$estimates))
})

test_that("print methods work for all estimator types", {
  d <- simulate_test_data(n = 100, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "all",
    n_boot = 0,
    verbose = FALSE
  )

  expect_output(print(results$results$always_ipw), "IPW")
  expect_output(print(results$results$always_gcomp), "G-comp")
  expect_output(print(results$results$always_tmle), "TMLE")
})
