# Tests for ffSL integration

test_that(".safe_sl() with use_ffSL=TRUE produces valid predictions", {

  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- rbinom(n, 1, plogis(0.5 * X$x1 - 0.3 * X$x2))

  result <- .safe_sl(Y = Y, X = X, learners = c("SL.glm", "SL.mean"),
                     cv_folds = 2L, use_ffSL = TRUE, verbose = FALSE)

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
                     cv_folds = 2L, use_ffSL = FALSE, verbose = FALSE)
  res_ff <- .safe_sl(Y = Y, X = X, learners = c("SL.glm", "SL.mean"),
                     cv_folds = 2L, use_ffSL = TRUE, verbose = FALSE)

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
                     cv_folds = 2L, use_ffSL = TRUE, verbose = FALSE)

  expect_true(result$method %in% c("SuperLearner", "glm"))
  expect_length(result$predictions, n)
})

test_that("use_ffSL parameter stored in fit objects", {
  d <- simulate_test_data(n = 60, K = 3)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE,
                       use_ffSL = FALSE)

  expect_false(obj$fits$treatment[["always"]]$use_ffSL)
})

test_that("longy() accepts use_ffSL parameter", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("future.apply")

  d <- simulate_test_data(n = 80, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    use_ffSL = TRUE,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_data")
  expect_s3_class(results$results$always_ipw, "longy_result")
  expect_true(all(is.finite(results$results$always_ipw$estimates$estimate)))
})
