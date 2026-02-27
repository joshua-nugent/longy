test_that("longy_data survives saveRDS/readRDS roundtrip", {
  skip_if_not_installed("SuperLearner")
  set.seed(42)
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L, never = 0L),
               verbose = FALSE)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  obj2 <- readRDS(tmp)

  # Class should be preserved
  expect_s3_class(obj2, "longy_data")

  # All methods should work
  expect_no_error(print(obj2))
  expect_no_error(weight_diagnostics(obj2))
  expect_no_error(prediction_diagnostics(obj2))
  expect_no_error(sl_diagnostics(obj2))
  expect_no_error(positivity_diagnostics(obj2))

  # Plot should dispatch correctly
  expect_no_error(plot(obj2))

  # Results should be intact
  expect_false(is.null(obj2$results$always_ipw))
  expect_equal(obj2$results$always_ipw$estimates$estimate,
               obj$results$always_ipw$estimates$estimate)
})

test_that("as_longy_data() restores class when lost after serialization", {
  skip_if_not_installed("SuperLearner")
  set.seed(42)
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L, never = 0L),
               verbose = FALSE)

  # Simulate class loss (strip class to plain list)
  obj_stripped <- unclass(obj)
  expect_false(inherits(obj_stripped, "longy_data"))

  # as_longy_data should detect by structure and restore
  obj_fixed <- as_longy_data(obj_stripped)
  expect_s3_class(obj_fixed, "longy_data")

  # All methods work on the restored object
  expect_no_error(weight_diagnostics(obj_fixed))
  expect_no_error(prediction_diagnostics(obj_fixed))
  expect_no_error(plot(obj_fixed))
})

test_that("as_longy_data() is a no-op on valid longy_data", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = c("W1", "W2"),
                    verbose = FALSE)

  obj2 <- as_longy_data(obj)
  expect_s3_class(obj2, "longy_data")
  expect_identical(obj$meta, obj2$meta)
})

test_that("as_longy_data() errors on non-longy objects", {
  expect_error(as_longy_data(list(a = 1)), "Expected")
  expect_error(as_longy_data(42), "Expected")
  expect_error(as_longy_data("foo"), "Expected")
})

test_that(".as_longy_data() repairs data.table selfref after deserialization", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = c("W1", "W2"),
                    verbose = FALSE)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  obj2 <- readRDS(tmp)

  # After as_longy_data repair, := should work without warning
  obj2 <- as_longy_data(obj2)
  expect_no_warning(obj2$data[, .tmp_test := 1L])
  obj2$data[, .tmp_test := NULL]
})

test_that("plot_sl_diagnostics works on object with lost class", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("ggplot2")
  set.seed(42)
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               verbose = FALSE)

  obj_stripped <- unclass(obj)
  expect_no_error(plot_sl_diagnostics(obj_stripped))
})
