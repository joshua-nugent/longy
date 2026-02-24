test_that("fit_censoring runs on data with censoring", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  expect_true(length(obj$fits$censoring[["always"]]) > 0)
  expect_true(nrow(obj$fits$censoring[["always"]][[".cens_censored"]]$predictions) > 0)
})

test_that("fit_censoring risk set uses full uncensored sample", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$censoring[["always"]][[".cens_censored"]]$predictions
  t0 <- preds[preds$.time == 0, ]
  # At time 0, risk set is ALL subjects (everyone uncensored at baseline)
  n_subjects <- length(unique(d$id))
  expect_equal(t0$.n_risk[1], n_subjects)
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

test_that("fit_censoring handles multiple censoring causes", {
  d <- simulate_multi_censoring(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  expect_true(".cens_death" %in% names(obj$fits$censoring[["always"]]))
  expect_true(".cens_ltfu" %in% names(obj$fits$censoring[["always"]]))
})

test_that("fit_censoring works with SuperLearner library", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  obj <- fit_censoring(obj, regime = "always",
                       learners = c("SL.glm", "SL.mean"),
                       verbose = FALSE)

  preds <- obj$fits$censoring[["always"]][[".cens_censored"]]$predictions
  expect_true(nrow(preds) > 0)

  sl_info <- obj$fits$censoring[["always"]][[".cens_censored"]]$sl_info
  expect_true(length(sl_info) > 0)
})
