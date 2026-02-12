test_that("estimate_ipw produces valid estimates", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  result <- estimate_ipw(obj, regime = "always")

  expect_s3_class(result, "longy_result")
  expect_true(nrow(result$estimates) > 0)
  # Estimates should be in [0, 1] for binary outcome
  expect_true(all(result$estimates$estimate >= 0 & result$estimates$estimate <= 1,
                  na.rm = TRUE))
})

test_that("IC inference produces non-negative SEs", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  result <- estimate_ipw(obj, regime = "always", inference = "ic")

  # SEs should be non-negative (can be 0 if all outcomes are same value)
  expect_true(all(result$estimates$se >= 0, na.rm = TRUE))
  # At least some SEs should be positive (non-degenerate time points)
  expect_true(any(result$estimates$se > 0, na.rm = TRUE))
  # CIs should be valid where SE > 0
  pos_se <- result$estimates$se > 0 & !is.na(result$estimates$se)
  if (any(pos_se)) {
    expect_true(all(
      result$estimates$ci_lower[pos_se] < result$estimates$ci_upper[pos_se]
    ))
  }
})

test_that("print and summary work for longy_result", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")
  result <- estimate_ipw(obj, regime = "always")

  expect_output(print(result), "longy IPW result")
  expect_output(summary(result), "longy IPW Result Summary")
})

test_that("estimate_ipw handles specific time points", {
  d <- simulate_test_data(n = 80, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  result <- estimate_ipw(obj, regime = "always", times = c(0, 2))

  expect_equal(nrow(result$estimates), 2)
  expect_equal(result$estimates$time, c(0, 2))
})

test_that("no-confounding recovery: IPW close to truth", {
  d <- simulate_no_confounding(n = 500, K = 3)
  # Make all outcomes observed for simplicity
  d$R <- 1L
  d$Y[is.na(d$Y)] <- rbinom(sum(is.na(d$Y)), 1, 0.3)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  result <- estimate_ipw(obj, regime = "always")

  # With no confounding, IPW estimate should be close to crude mean among treated
  crude_mean_t0 <- mean(d$Y[d$time == 0 & d$A == 1], na.rm = TRUE)
  ipw_t0 <- result$estimates$estimate[result$estimates$time == 0]

  # Should be within 0.15 (Monte Carlo noise with n=500)
  expect_true(abs(ipw_t0 - crude_mean_t0) < 0.15,
              info = sprintf("IPW=%.3f, crude=%.3f", ipw_t0, crude_mean_t0))
})
