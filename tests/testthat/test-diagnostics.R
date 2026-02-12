test_that("weight_diagnostics produces correct output", {
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

  diag <- weight_diagnostics(obj, by_time = TRUE)

  expect_true(data.table::is.data.table(diag))
  expect_true("mean_weight" %in% names(diag))
  expect_true("ess" %in% names(diag))
  expect_true(all(diag$ess > 0))
  expect_true(all(diag$mean_weight > 0))
})

test_that("weight_diagnostics overall summary works", {
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

  diag <- weight_diagnostics(obj, by_time = FALSE)

  expect_equal(nrow(diag), 1)
})

test_that("ESS formula is correct", {
  # Kish's ESS = (sum(w))^2 / sum(w^2)
  w <- c(1, 1, 1, 1, 1)
  expect_equal(longy:::.ess(w), 5)

  w2 <- c(2, 2, 2, 2, 2)
  expect_equal(longy:::.ess(w2), 5)

  # Extreme weights should reduce ESS
  w3 <- c(1, 1, 1, 1, 100)
  ess3 <- longy:::.ess(w3)
  expect_true(ess3 < 5)
})

test_that("positivity_diagnostics detects extreme propensities", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  # With a very strict threshold, should flag more observations
  expect_message(
    diag <- positivity_diagnostics(obj, threshold = 0.5),
    "flagged"
  )
})

test_that("weight_diagnostics errors without weights", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(weight_diagnostics(obj), "not computed")
})
