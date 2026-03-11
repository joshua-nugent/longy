test_that("fit_outcome runs without error on basic data", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  expect_false(is.null(obj$fits$outcome[["always"]]))
  expect_true(nrow(obj$fits$outcome[["always"]]$predictions) > 0)
})

test_that("fit_outcome predictions are bounded for binary outcomes", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  bounds <- c(0.01, 0.99)
  obj <- fit_outcome(obj, regime = "always", bounds = bounds, verbose = FALSE)

  preds <- obj$fits$outcome[["always"]]$predictions$.Q_hat
  expect_true(all(preds >= bounds[1]))
  expect_true(all(preds <= bounds[2]))
})

test_that("fit_outcome stores correct family for binary vs continuous", {
  # Binary
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  expect_equal(obj$fits$outcome[["always"]]$family, "quasibinomial")

  # Continuous
  d2 <- simulate_continuous_outcome(n = 80, K = 3)
  obj2 <- longy_data(d2, id = "id", time = "time", outcome = "Y",
                     treatment = "A", censoring = "C", observation = "R",
                     baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                     outcome_type = "continuous", verbose = FALSE)
  obj2 <- define_regime(obj2, "always", static = 1L)
  obj2 <- fit_outcome(obj2, regime = "always", verbose = FALSE)
  expect_equal(obj2$fits$outcome[["always"]]$family, "gaussian")
})

test_that("fit_outcome errors on missing regime", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    fit_outcome(obj, regime = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("fit_outcome refit protection works", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  # Should error without refit=TRUE
  expect_error(
    fit_outcome(obj, regime = "always", verbose = FALSE),
    "refit"
  )

  # Should succeed with refit=TRUE
  expect_no_error(
    fit_outcome(obj, regime = "always", verbose = FALSE, refit = TRUE)
  )
})

test_that("fit_outcome works with survival outcomes", {
  d <- simulate_survival_outcome(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$outcome[["always"]]$predictions
  expect_true(nrow(preds) > 0)
  # Predictions should be in [0,1] for survival
  expect_true(all(preds$.Q_hat >= 0 & preds$.Q_hat <= 1))
})

test_that("fit_outcome works with competing risks", {
  d <- simulate_competing_risks(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", competing = "D",
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$outcome[["always"]]$predictions
  expect_true(nrow(preds) > 0)
  # Predictions should be in [0,1]
  expect_true(all(preds$.Q_hat >= 0 & preds$.Q_hat <= 1))
})

test_that("fit_outcome handles multiple regimes", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_outcome(obj, regime = c("always", "never"), verbose = FALSE)

  expect_false(is.null(obj$fits$outcome[["always"]]))
  expect_false(is.null(obj$fits$outcome[["never"]]))
  expect_true(nrow(obj$fits$outcome[["always"]]$predictions) > 0)
  expect_true(nrow(obj$fits$outcome[["never"]]$predictions) > 0)
})

test_that("fit_outcome risk_set='followers' produces valid fits", {
  d <- simulate_test_data(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_outcome(obj, regime = "always", risk_set = "followers",
                     verbose = FALSE)

  expect_false(is.null(obj$fits$outcome[["always"]]))
  expect_equal(obj$fits$outcome[["always"]]$risk_set, "followers")
  expect_true(nrow(obj$fits$outcome[["always"]]$predictions) > 0)
})

test_that("fit_outcome times parameter restricts target times", {
  d <- simulate_test_data(n = 100, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  # Fit only at times 2 and 4
  obj <- fit_outcome(obj, regime = "always", times = c(2, 4),
                     verbose = FALSE)

  preds <- obj$fits$outcome[["always"]]$predictions
  target_times <- unique(preds$.target_time)
  expect_true(all(target_times %in% c(2, 4)))
  expect_true(length(target_times) <= 2)
})

test_that("fit_outcome metadata_only stores metadata without predictions", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_outcome(obj, regime = "always", metadata_only = TRUE,
                     verbose = FALSE)

  fit <- obj$fits$outcome[["always"]]
  expect_true(isTRUE(fit$metadata_only))
  expect_null(fit$predictions)
  # Metadata fields should be present
  expect_false(is.null(fit$covariates))
  expect_false(is.null(fit$bounds))
  expect_false(is.null(fit$risk_set))
})

test_that("estimate_gcomp errors on metadata_only outcome fit", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", metadata_only = TRUE,
                     verbose = FALSE)

  expect_error(
    estimate_gcomp(obj, regime = "always", n_boot = 0, verbose = FALSE),
    "metadata only"
  )
})

test_that("fit_outcome does not leave tracking columns on obj$data", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  cols_before <- names(obj$data)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  cols_after <- names(obj$data)

  # No .longy_ tracking columns should be added (except .longy_fold/.longy_lag_*)
  new_cols <- setdiff(cols_after, cols_before)
  tracking_leaked <- grep("^\\.longy_(regime|consist|uncens|cum|Q|first)",
                          new_cols, value = TRUE)
  expect_length(tracking_leaked, 0)
})

test_that("fit_outcome sl_info has entries for each backward step", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", times = c(3), verbose = FALSE)

  sl_info <- obj$fits$outcome[["always"]]$sl_info
  expect_true(length(sl_info) > 0)

  # Each entry should have the expected diagnostic fields
  for (entry in sl_info) {
    expect_true("time" %in% names(entry))
    expect_true("target_time" %in% names(entry))
    expect_true("method" %in% names(entry))
    expect_true("n_risk" %in% names(entry))
    expect_true("n_train" %in% names(entry))
  }
})
