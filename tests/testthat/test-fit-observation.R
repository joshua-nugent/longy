test_that("fit_observation runs with intermittent data", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  expect_false(is.null(obj$fits$observation[["always"]]))
  expect_true(nrow(obj$fits$observation[["always"]]$predictions) > 0)
})

test_that("fit_observation risk set is most restrictive", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  preds_r <- obj$fits$observation[["always"]]$predictions
  preds_c <- obj$fits$censoring[["always"]][[".cens_censored"]]$predictions

  # g_R risk set should be <= g_C risk set at each time
  for (tt in unique(preds_r$.time)) {
    n_r <- preds_r[preds_r$.time == tt, ]$.n_risk[1]
    n_c_rows <- preds_c[preds_c$.time == tt, ]
    if (nrow(n_c_rows) > 0) {
      n_c <- n_c_rows$.n_risk[1]
      expect_true(n_r <= n_c,
                  info = sprintf("time %d: g_R risk (%d) > g_C risk (%d)", tt, n_r, n_c))
    }
  }
})

test_that("fit_observation skips when observation=NULL", {
  d <- simulate_always_observed(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  expect_null(obj$fits$observation[["always"]])
})

test_that("fit_observation predictions are bounded", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  bounds <- c(0.01, 0.99)
  obj <- fit_observation(obj, regime = "always", bounds = bounds, verbose = FALSE)

  preds <- obj$fits$observation[["always"]]$predictions$.p_r
  expect_true(all(preds >= bounds[1]))
  expect_true(all(preds <= bounds[2]))
})

test_that("fit_observation works with SuperLearner library", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  obj <- fit_observation(obj, regime = "always",
                         learners = c("SL.glm", "SL.mean"),
                         verbose = FALSE)

  preds <- obj$fits$observation[["always"]]$predictions
  expect_true(nrow(preds) > 0)

  sl_info <- obj$fits$observation[["always"]]$sl_info
  expect_true(length(sl_info) > 0)
})
