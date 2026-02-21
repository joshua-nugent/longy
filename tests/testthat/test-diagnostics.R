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

# --- Tests for longy_result / longy_results dispatch ---

test_that("weight_diagnostics works with longy_result from estimate_ipw", {
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
  result <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                         inference = "none")

  # longy_result should work
  diag <- weight_diagnostics(result)
  expect_true(data.table::is.data.table(diag))
  expect_true("ess" %in% names(diag))
})

test_that("weight_diagnostics works with longy_results from longy()", {
  d <- simulate_test_data(n = 80, K = 4)
  results <- suppressWarnings(longy(
    d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "ipw", inference = "none", verbose = FALSE
  ))

  # longy_results should work (takes first element)
  diag <- weight_diagnostics(results)
  expect_true(data.table::is.data.table(diag))
  expect_true("ess" %in% names(diag))
})

test_that("positivity_diagnostics works with longy_result", {
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
  result <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                         inference = "none")

  expect_message(
    diag <- positivity_diagnostics(result, threshold = 0.5),
    "flagged"
  )
})

test_that("diagnostics error on invalid input", {
  expect_error(.extract_longy_data("not an object"),
               "Expected a longy_data")
  expect_error(.extract_longy_data(list(a = 1)),
               "Expected a longy_data")
})

# --- Tests for sl_diagnostics ---

test_that("sl_diagnostics returns correct structure for treatment model", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "treatment")

  expect_true(data.table::is.data.table(diag))
  expect_true(all(c("model", "time", "method") %in% names(diag)))
  expect_true(all(diag$model == "treatment"))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics returns all models", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "all")

  expect_true(data.table::is.data.table(diag))
  models_present <- unique(diag$model)
  expect_true("treatment" %in% models_present)
  expect_true("censoring" %in% models_present)
  expect_true("observation" %in% models_present)
})

test_that("sl_diagnostics shows censoring submodel names", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "censoring")
  expect_true(all(diag$submodel == ".cens_censored"))
})

test_that("sl_diagnostics works with longy_result", {
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
  result <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                         inference = "none")

  diag <- sl_diagnostics(result)
  expect_true(data.table::is.data.table(diag))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics works with longy_results from longy()", {
  d <- simulate_test_data(n = 80, K = 4)
  results <- suppressWarnings(longy(
    d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "ipw", inference = "none", verbose = FALSE
  ))

  diag <- sl_diagnostics(results)
  expect_true(data.table::is.data.table(diag))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics returns empty table when no fits", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_message(diag <- sl_diagnostics(obj), "No model fits found")
  expect_equal(nrow(diag), 0)
})
