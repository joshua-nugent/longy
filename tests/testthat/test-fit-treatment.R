test_that("fit_treatment runs without error on basic data", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  expect_false(is.null(obj$fits$treatment))
  expect_true(nrow(obj$fits$treatment$predictions) > 0)
})

test_that("fit_treatment predictions are bounded", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  bounds <- c(0.01, 0.99)
  obj <- fit_treatment(obj, regime = "always", bounds = bounds, verbose = FALSE)

  preds <- obj$fits$treatment$predictions$.p_a
  expect_true(all(preds >= bounds[1]))
  expect_true(all(preds <= bounds[2]))
})

test_that("fit_treatment uses correct risk set", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$treatment$predictions
  # At time 0, everyone should be at risk
  t0 <- preds[preds$.time == 0, ]
  expect_equal(t0$.n_risk[1], sum(d$time == 0))

  # At later times, risk set should be smaller (regime-consistent + uncensored)
  t1 <- preds[preds$.time == 1, ]
  if (nrow(t1) > 0) {
    expect_true(t1$.n_risk[1] <= t0$.n_risk[1])
  }
})

test_that("fit_treatment falls back to marginal with constant treatment", {
  d <- simulate_test_data(n = 50, K = 3)
  # Make treatment constant at time 0
  d$A[d$time == 0] <- 1L
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  t0 <- obj$fits$treatment$predictions[obj$fits$treatment$predictions$.time == 0, ]
  # When all Y=1, marginal method should be used
  expect_equal(unique(t0$.method), "marginal")
})

test_that("fit_treatment works with SuperLearner library", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always",
                       learners = c("SL.glm", "SL.mean"),
                       verbose = FALSE)

  preds <- obj$fits$treatment$predictions
  expect_true(nrow(preds) > 0)

  # Should use SuperLearner, not marginal
  methods <- unique(preds$.method)
  expect_true(any(methods %in% c("SuperLearner", "glm")))

  # sl_info should contain risk and coef
  sl_info <- obj$fits$treatment$sl_info
  expect_true(length(sl_info) > 0)
  sl_entry <- sl_info[[1]]
  if (sl_entry$method %in% c("SuperLearner", "glm")) {
    expect_false(is.null(sl_entry$sl_risk))
    expect_false(is.null(sl_entry$sl_coef))
  }
})

test_that("fit_treatment errors on missing regime", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    fit_treatment(obj, regime = "nonexistent", verbose = FALSE),
    "not found"
  )
})
