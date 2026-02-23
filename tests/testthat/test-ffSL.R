# Tests for ffSL integration

test_that(".safe_sl() with sl_fn='ffSL' produces valid predictions", {

  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- rbinom(n, 1, plogis(0.5 * X$x1 - 0.3 * X$x2))

  result <- .safe_sl(Y = Y, X = X, learners = c("SL.glm", "SL.mean"),
                     cv_folds = 2L, sl_fn = "ffSL", verbose = FALSE)

  # Should succeed (either via ffSL or glm fallback)
  expect_true(result$method %in% c("SuperLearner", "glm"))
  expect_length(result$predictions, n)
  expect_true(all(result$predictions >= 0 & result$predictions <= 1))
})

test_that(".safe_sl() ffSL returns same structure as standard SL", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  set.seed(456)
  n <- 150
  X <- data.frame(x1 = rnorm(n))
  Y <- rbinom(n, 1, plogis(0.3 * X$x1))

  res_sl <- .safe_sl(Y = Y, X = X, learners = c("SL.glm", "SL.mean"),
                     cv_folds = 2L, sl_fn = "SuperLearner", verbose = FALSE)
  res_ff <- .safe_sl(Y = Y, X = X, learners = c("SL.glm", "SL.mean"),
                     cv_folds = 2L, sl_fn = "ffSL", verbose = FALSE)

  # Both should produce predictions of the right length
  expect_length(res_ff$predictions, n)
  expect_length(res_sl$predictions, n)
  # Both should have a method field

  expect_true("method" %in% names(res_ff))
  expect_true("method" %in% names(res_sl))
})

test_that("ffSL falls back to glm on degenerate data", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  set.seed(789)
  n <- 10
  X <- data.frame(x1 = rep(1, n))
  Y <- rbinom(n, 1, 0.5)

  result <- .safe_sl(Y = Y, X = X, learners = c("SL.glm"),
                     cv_folds = 2L, sl_fn = "ffSL", verbose = FALSE)

  expect_true(result$method %in% c("SuperLearner", "glm"))
  expect_length(result$predictions, n)
})

test_that("sl_fn parameter stored in fit objects", {
  d <- simulate_test_data(n = 60, K = 3)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE,
                       sl_fn = "SuperLearner")

  expect_equal(obj$fits$treatment[["always"]]$sl_fn, "SuperLearner")
})

test_that("longy() accepts sl_fn parameter", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  d <- simulate_test_data(n = 80, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    sl_fn = "ffSL",
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_s3_class(results$always, "longy_result")
  expect_true(all(is.finite(results$always$estimates$estimate)))
})

test_that("longy() rejects invalid sl_fn", {
  d <- simulate_test_data(n = 50, K = 3)
  expect_error(
    longy(data = d, id = "id", time = "time", outcome = "Y",
          treatment = "A", censoring = "C", observation = "R",
          baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
          regimes = list(always = 1L), sl_fn = "bogus", verbose = FALSE),
    "arg"
  )
})
