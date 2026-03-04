# Tests for outer-level parallelization

test_that(".parallel_or_sequential() falls back to sequential when parallel=FALSE", {
  tasks <- 1:5
  result <- .parallel_or_sequential(tasks, function(x) x^2, parallel = FALSE)
  expect_equal(result, as.list((1:5)^2))
})

test_that(".parallel_or_sequential() falls back to sequential without future", {
  # Even with parallel=TRUE, falls back when future::plan() is sequential
  tasks <- 1:3
  result <- .parallel_or_sequential(tasks, function(x) x + 1, parallel = TRUE)
  expect_equal(result, as.list(2:4))
})

test_that("fit_treatment() works with parallel=TRUE (sequential fallback)", {
  d <- simulate_test_data(n = 60, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  # parallel=TRUE but no future plan set -> sequential fallback
  obj_par <- fit_treatment(obj, regime = "always", verbose = FALSE,
                           parallel = TRUE)
  expect_true(!is.null(obj_par$fits$treatment[["always"]]))
  expect_true(nrow(obj_par$fits$treatment[["always"]]$predictions) > 0)
})

test_that("fit_censoring() works with parallel=TRUE (sequential fallback)", {
  d <- simulate_test_data(n = 60, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj_par <- fit_censoring(obj, regime = "always", verbose = FALSE,
                           parallel = TRUE)
  expect_true(!is.null(obj_par$fits$censoring[["always"]]))
})

test_that("fit_observation() works with parallel=TRUE (sequential fallback)", {
  d <- simulate_test_data(n = 60, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj_par <- fit_observation(obj, regime = "always", verbose = FALSE,
                             parallel = TRUE)
  expect_true(!is.null(obj_par$fits$observation[["always"]]))
})

test_that("fit_outcome() works with parallel=TRUE (sequential fallback)", {
  d <- simulate_test_data(n = 60, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj_par <- fit_outcome(obj, regime = "always", verbose = FALSE,
                         parallel = TRUE)
  expect_true(!is.null(obj_par$fits$outcome[["always"]]))
  expect_true(nrow(obj_par$fits$outcome[["always"]]$predictions) > 0)
})

test_that("longy() parallel=TRUE with sequential plan produces valid results", {
  d <- simulate_test_data(n = 80, K = 3, seed = 42)

  # parallel=TRUE but no future plan => sequential fallback
  obj_par <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    parallel = TRUE,
    verbose = FALSE
  )

  expect_s3_class(obj_par, "longy_data")
  expect_s3_class(obj_par$results$always_ipw, "longy_result")
  expect_true(all(is.finite(obj_par$results$always_ipw$estimates$estimate)))
})

test_that("G-comp bootstrap conflict warning fires", {
  d <- simulate_test_data(n = 80, K = 3)

  expect_warning(
    longy(
      data = d,
      id = "id", time = "time", outcome = "Y",
      treatment = "A", censoring = "C", observation = "R",
      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
      regimes = list(always = 1L),
      estimator = "gcomp",
      n_boot = 10,
      parallel = TRUE,
      verbose = FALSE
    ),
    "parallel=TRUE with G-comp bootstrap"
  )
})

test_that("parallel=TRUE with multisession runs without error", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  d <- simulate_test_data(n = 80, K = 3, seed = 42)

  # Run with actual parallel workers
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)

  obj_par <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    parallel = TRUE,
    verbose = FALSE
  )

  # Should complete successfully and produce finite estimates
  expect_s3_class(obj_par, "longy_data")
  expect_true(all(is.finite(obj_par$results$always_ipw$estimates$estimate)))
})
