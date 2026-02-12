test_that("fit_censoring runs on data with censoring", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  expect_true(length(obj$fits$censoring) > 0)
  expect_true(nrow(obj$fits$censoring[["C"]]$predictions) > 0)
})

test_that("fit_censoring risk set conditions on A(t)", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$censoring[["C"]]$predictions
  t0 <- preds[preds$.time == 0, ]
  # At time 0, risk set should only include those with A=1 (for "always" regime)
  n_treated_t0 <- sum(d$A[d$time == 0] == 1)
  expect_equal(t0$.n_risk[1], n_treated_t0)
})

test_that("fit_censoring handles no-censoring data gracefully", {
  d <- simulate_no_censoring(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  expect_equal(length(obj$fits$censoring), 0)
})

test_that("fit_censoring handles multiple censoring sources", {
  d <- simulate_test_data(n = 80, K = 4)
  # Add a second censoring variable (administrative)
  d$C2 <- rbinom(nrow(d), 1, 0.02)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = c("C", "C2"),
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  expect_true("C" %in% names(obj$fits$censoring))
  expect_true("C2" %in% names(obj$fits$censoring))
})
