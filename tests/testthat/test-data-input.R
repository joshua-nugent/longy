test_that("longy_data creates valid object from simple data", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(
    data = d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    verbose = FALSE
  )

  expect_s3_class(obj, "longy_data")
  expect_true(data.table::is.data.table(obj$data))
  expect_equal(obj$nodes$id, "id")
  expect_equal(obj$nodes$treatment, "A")
  expect_equal(obj$nodes$censoring, "C")
  expect_equal(obj$nodes$observation, "R")
  expect_equal(obj$meta$n_subjects, 50)
  expect_equal(obj$meta$min_time, 0)
})

test_that("longy_data rejects missing columns", {
  d <- simulate_test_data(n = 20, K = 3)
  expect_error(
    longy_data(d, id = "id", time = "time", outcome = "Y",
               treatment = "MISSING_COL", verbose = FALSE),
    "not found"
  )
})

test_that("longy_data rejects non-binary treatment", {
  d <- simulate_test_data(n = 20, K = 3)
  d$A <- d$A + 1  # now {1, 2}
  expect_error(
    longy_data(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", verbose = FALSE),
    "binary"
  )
})

test_that("longy_data rejects duplicate id-time pairs", {
  d <- simulate_test_data(n = 20, K = 3)
  d <- rbind(d, d[1, ])  # duplicate
  expect_error(
    longy_data(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", verbose = FALSE),
    "Duplicate"
  )
})

test_that("longy_data converts factors to dummies", {
  d <- simulate_test_data(n = 30, K = 3)
  d$W_cat <- sample(c("low", "med", "high"), nrow(d), replace = TRUE)
  # Make W_cat constant within subject
  for (i in unique(d$id)) {
    val <- d$W_cat[d$id == i][1]
    d$W_cat[d$id == i] <- val
  }
  obj <- longy_data(
    data = d, id = "id", time = "time", outcome = "Y",
    treatment = "A", baseline = c("W1", "W2", "W_cat"),
    timevarying = c("L1", "L2"), verbose = FALSE
  )
  # Original factor column should be gone; dummies should exist
  expect_false("W_cat" %in% names(obj$data))
  dummy_cols <- grep("^W_cat_", names(obj$data), value = TRUE)
  expect_true(length(dummy_cols) > 0)
})

test_that("longy_data rejects non-constant baseline", {
  d <- simulate_test_data(n = 20, K = 3)
  # Make W1 vary within a subject
  idx <- which(d$id == 1)
  if (length(idx) > 1) d$W1[idx[2]] <- d$W1[idx[1]] + 10
  expect_error(
    longy_data(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", baseline = c("W1", "W2"), verbose = FALSE),
    "varies within"
  )
})

test_that("longy_data works without censoring or observation", {
  d <- simulate_no_censoring(n = 30, K = 3)
  d$C <- NULL
  d$R <- NULL
  d$Y[is.na(d$Y)] <- 0  # fill NAs
  obj <- longy_data(
    data = d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = NULL, observation = NULL,
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    verbose = FALSE
  )
  expect_s3_class(obj, "longy_data")
  expect_null(obj$nodes$censoring)
  expect_null(obj$nodes$observation)
})

test_that("print.longy_data runs without error", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  expect_output(print(obj), "longy_data object")
})

test_that("summary.longy_data runs without error", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  expect_output(summary(obj), "longy_data summary")
})

test_that("longy_data rejects continuous Y when outcome_type = 'binary'", {
  d <- simulate_continuous_outcome(n = 20, K = 3)
  expect_error(
    longy_data(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               outcome_type = "binary", verbose = FALSE),
    "binary.*0.*1"
  )
})

test_that("set_crossfit assigns folds at subject level", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- set_crossfit(obj, n_folds = 5, seed = 1)

  expect_true(obj$crossfit$enabled)
  expect_equal(obj$crossfit$n_folds, 5)
  expect_true(".longy_fold" %in% names(obj$data))

  # Each subject should have one fold value
  fold_by_id <- obj$data[, .(nfold = data.table::uniqueN(.longy_fold)), by = id]
  expect_true(all(fold_by_id$nfold == 1))
})
