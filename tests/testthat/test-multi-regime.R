# Tests for multi-regime support in sub-functions

test_that("fit_treatment fits multiple regimes in one call", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)

  expect_false(is.null(obj$fits$treatment[["always"]]))
  expect_false(is.null(obj$fits$treatment[["never"]]))
  expect_true(nrow(obj$fits$treatment[["always"]]$predictions) > 0)
  expect_true(nrow(obj$fits$treatment[["never"]]$predictions) > 0)
})

test_that("sequential single-regime calls preserve earlier fits", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_treatment(obj, regime = "never", verbose = FALSE)

  # Both should still be present
  expect_false(is.null(obj$fits$treatment[["always"]]))
  expect_false(is.null(obj$fits$treatment[["never"]]))
})

test_that("default regime = NULL fits all defined regimes", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, verbose = FALSE)

  expect_equal(sort(names(obj$fits$treatment)), c("always", "never"))
})

test_that("missing regime errors clearly", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  expect_error(
    fit_treatment(obj, regime = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("multi-regime estimate returns longy_results", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_censoring(obj, verbose = FALSE)
  obj <- fit_observation(obj, verbose = FALSE)
  obj <- compute_weights(obj)

  results <- estimate_ipw(obj)

  expect_s3_class(results, "longy_results")
  expect_equal(length(results), 2)
  expect_s3_class(results$always, "longy_result")
  expect_s3_class(results$never, "longy_result")
})

test_that("single-regime estimate returns longy_result", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_censoring(obj, verbose = FALSE)
  obj <- fit_observation(obj, verbose = FALSE)
  obj <- compute_weights(obj)

  result <- estimate_ipw(obj, regime = "always")

  expect_s3_class(result, "longy_result")
})

test_that("compute_weights errors if treatment fit missing for regime", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  expect_error(
    compute_weights(obj, regime = "always"),
    "Treatment model not fit"
  )
})

test_that("multi-regime G-comp pipeline works", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_outcome(obj, verbose = FALSE)
  results <- estimate_gcomp(obj, n_boot = 0, verbose = FALSE)

  expect_s3_class(results, "longy_results")
  expect_equal(length(results), 2)
  expect_true(all(is.finite(results$always$estimates$estimate)))
  expect_true(all(is.finite(results$never$estimates$estimate)))
})

test_that("multi-regime TMLE pipeline works", {
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_censoring(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)

  results <- estimate_tmle(obj, inference = "none", verbose = FALSE)

  expect_s3_class(results, "longy_results")
  expect_equal(length(results), 2)
  expect_equal(results$always$estimator, "tmle")
  expect_equal(results$never$estimator, "tmle")
})

test_that(".resolve_regimes works correctly", {
  d <- simulate_test_data(n = 20, K = 2)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  # NULL returns all
  expect_equal(sort(longy:::.resolve_regimes(obj, NULL)),
               c("always", "never"))
  # Specific regime

  expect_equal(longy:::.resolve_regimes(obj, "always"), "always")
  # Bad regime errors
  expect_error(longy:::.resolve_regimes(obj, "bogus"), "not found")
  # Multiple with one bad errors
  expect_error(longy:::.resolve_regimes(obj, c("always", "bogus")), "not found")
})

test_that("re-fitting same regime overwrites previous fit", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  preds1 <- obj$fits$treatment[["always"]]$predictions

  # Re-fit (same regime) should overwrite
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  preds2 <- obj$fits$treatment[["always"]]$predictions

  # Should have the same structure (regime fit overwritten)
  expect_equal(nrow(preds1), nrow(preds2))
})

test_that("diagnostics work with multi-regime object", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_censoring(obj, verbose = FALSE)
  obj <- fit_observation(obj, verbose = FALSE)
  obj <- compute_weights(obj)

  # Weight diagnostics for specific regime
  diag_a <- weight_diagnostics(obj, regime = "always")
  diag_n <- weight_diagnostics(obj, regime = "never")
  expect_true(data.table::is.data.table(diag_a))
  expect_true(data.table::is.data.table(diag_n))

  # Positivity diagnostics for specific regime
  expect_message(
    positivity_diagnostics(obj, regime = "always", threshold = 0.5),
    "flagged"
  )

  # SL diagnostics for specific regime
  sl_a <- sl_diagnostics(obj, regime = "always")
  sl_n <- sl_diagnostics(obj, regime = "never")
  expect_true(data.table::is.data.table(sl_a))
  expect_true(data.table::is.data.table(sl_n))
})
