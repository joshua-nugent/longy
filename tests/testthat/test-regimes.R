test_that("define_regime creates static regimes", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  obj <- define_regime(obj, name = "always", static = 1L)
  obj <- define_regime(obj, name = "never", static = 0L)

  expect_equal(length(obj$regimes), 2)
  expect_equal(obj$regimes$always$type, "static")
  expect_equal(obj$regimes$always$value, 1L)
  expect_equal(obj$regimes$never$value, 0L)
})

test_that("define_regime rejects duplicate names", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, name = "always", static = 1L)

  expect_error(
    define_regime(obj, name = "always", static = 0L),
    "already exists"
  )
})

test_that("define_regime creates dynamic regimes", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  dyn_fn <- function(row) as.integer(row[["L1"]] > 0)
  obj <- define_regime(obj, name = "dynamic_L1", dynamic = dyn_fn)

  expect_equal(obj$regimes$dynamic_L1$type, "dynamic")
  expect_true(is.function(obj$regimes$dynamic_L1$value))
})

test_that("define_regime creates stochastic regimes", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  stoch_fn <- function(row) 0.7  # 70% chance of treatment
  obj <- define_regime(obj, name = "stoch70", stochastic = stoch_fn)

  expect_equal(obj$regimes$stoch70$type, "stochastic")
})

test_that("define_regime rejects invalid static values", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(define_regime(obj, name = "bad", static = 2L), "0 or 1")
})

test_that("define_regime requires exactly one type", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    define_regime(obj, name = "none"),
    "Exactly one"
  )
  expect_error(
    define_regime(obj, name = "both", static = 1L, dynamic = identity),
    "Exactly one"
  )
})

test_that("dynamic regime evaluation error gives helpful message", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  bad_fn <- function(row) stop("custom error in regime")
  obj <- define_regime(obj, name = "bad_dyn", dynamic = bad_fn)

  expect_error(
    longy:::.evaluate_regime(obj$regimes$bad_dyn, obj$data),
    "Regime.*bad_dyn.*failed at row.*custom error"
  )
})

test_that(".evaluate_regime works for static regimes", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, name = "always", static = 1L)

  vals <- longy:::.evaluate_regime(obj$regimes$always, obj$data)
  expect_length(vals, nrow(obj$data))
  expect_true(all(vals == 1))
})
