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

test_that("fit_observation refit protection works", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  # refit=FALSE (default) should error

  expect_error(
    fit_observation(obj, regime = "always", verbose = FALSE),
    "already fitted"
  )

  # refit=TRUE should succeed
  expect_no_error(
    fit_observation(obj, regime = "always", refit = TRUE, verbose = FALSE)
  )
})

test_that("fit_observation multi-regime shares identical predictions", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)

  preds_a <- obj$fits$observation[["always"]]$predictions
  preds_n <- obj$fits$observation[["never"]]$predictions

  # Predictions should be identical (model fit once, copied)
  expect_equal(preds_a$.p_r, preds_n$.p_r)
  expect_equal(preds_a$id, preds_n$id)

  # Regime labels should differ
  expect_equal(obj$fits$observation[["always"]]$regime, "always")
  expect_equal(obj$fits$observation[["never"]]$regime, "never")
})

test_that("fit_observation handles zero at-risk data without crashing", {
  # Create data where observation indicator is NA everywhere -> 0 at risk for g_R
  set.seed(99)
  n <- 30
  d <- data.frame(
    id = rep(seq_len(n), each = 3),
    time = rep(0:2, n),
    W1 = rep(rnorm(n), each = 3),
    W2 = rep(rbinom(n, 1, 0.5), each = 3),
    L1 = rnorm(n * 3),
    L2 = rbinom(n * 3, 1, 0.5),
    A = rbinom(n * 3, 1, 0.5),
    C = "uncensored",
    R = NA_integer_,        # no observation data -> all filtered out
    Y = NA_real_,
    stringsAsFactors = FALSE
  )

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  # Should warn but not crash (rbindlist on empty list)
  expect_warning(
    obj <- fit_observation(obj, regime = "always", verbose = FALSE),
    "No observations at risk"
  )

  # Predictions should be an empty data.table, not NULL
  expect_true(nrow(obj$fits$observation[["always"]]$predictions) == 0)
})

test_that("fit_observation constant R=1 warns appropriately", {
  # Everyone always observed — g_R has no variation
  d <- simulate_always_observed(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  # Should warn about marginal fallback (constant outcome summarized as "constant")
  expect_warning(
    obj <- fit_observation(obj, regime = "always", verbose = FALSE),
    "marginal fallback"
  )

  # All predictions should use marginal method
  preds <- obj$fits$observation[["always"]]$predictions
  expect_true(all(preds$.method == "marginal"))
})

test_that("fit_observation rare-events marginal fallback", {
  # Create data where observation is extremely rare (< min_events)
  set.seed(42)
  n <- 200
  K <- 3
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    for (tt in 0:(K - 1)) {
      L1 <- rnorm(1)
      A <- rbinom(1, 1, 0.5)
      # Very rare observation: ~0.5% rate
      R <- rbinom(1, 1, 0.005)
      Y <- if (R == 1) rbinom(1, 1, 0.3) else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, L1 = L1,
        A = A, C = "uncensored", R = as.integer(R), Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }
  d <- do.call(rbind, rows)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = "W1", timevarying = "L1",
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  # Should trigger minority-class fallback
  expect_warning(
    obj <- fit_observation(obj, regime = "always", verbose = FALSE),
    "minority-class"
  )

  preds <- obj$fits$observation[["always"]]$predictions
  marginal_rows <- preds[preds$.method == "marginal", ]
  expect_true(nrow(marginal_rows) > 0)
})
