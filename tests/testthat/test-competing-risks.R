# Tests for competing risks support in longy

test_that("longy_data() accepts competing column and stores in nodes", {
  dat <- simulate_competing_risks(n = 50, K = 5, seed = 111)
  obj <- longy_data(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", competing = "D", verbose = FALSE
  )
  expect_s3_class(obj, "longy_data")
  expect_equal(obj$nodes$competing, "D")
  expect_true(is.integer(obj$data$D))
})

test_that("competing with non-survival outcome_type errors", {
  dat <- simulate_competing_risks(n = 50, K = 5, seed = 111)
  expect_error(
    longy_data(
      data = dat, id = "id", time = "time", outcome = "Y",
      treatment = "A", censoring = "C",
      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
      outcome_type = "binary", competing = "D", verbose = FALSE
    ),
    "outcome_type = 'survival'"
  )
})

test_that("non-absorbing competing column errors", {
  dat <- simulate_competing_risks(n = 50, K = 5, seed = 111)
  # Break absorbing property: set a D=1 row back to D=0
  ids_with_d <- unique(dat$id[dat$D == 1 & !is.na(dat$D)])
  if (length(ids_with_d) > 0) {
    bad_id <- ids_with_d[1]
    bad_rows <- which(dat$id == bad_id & dat$D == 1 & !is.na(dat$D))
    if (length(bad_rows) > 1) {
      # Set the last D=1 row to 0 (breaks absorbing)
      dat$D[bad_rows[length(bad_rows)]] <- 0L
    } else {
      # Only 1 D=1 row; add another row after it with D=0
      # Find a row after the event and set to 0
      all_rows_id <- which(dat$id == bad_id)
      after <- all_rows_id[all_rows_id > bad_rows[1]]
      if (length(after) > 0) {
        dat$D[after[1]] <- 0L
      }
    }
  }
  expect_error(
    longy_data(
      data = dat, id = "id", time = "time", outcome = "Y",
      treatment = "A", censoring = "C",
      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
      outcome_type = "survival", competing = "D", verbose = FALSE
    ),
    "absorbing"
  )
})

test_that("simultaneous Y=1 and D=1 errors", {
  dat <- simulate_competing_risks(n = 50, K = 5, seed = 111)
  # Force a row where both Y=1 and D=1 while keeping D absorbing
  first_y1 <- which(dat$Y == 1 & !is.na(dat$Y))[1]
  if (!is.na(first_y1)) {
    bad_id <- dat$id[first_y1]
    bad_time <- dat$time[first_y1]
    # Set D=1 for this row and all subsequent rows for this subject
    subsequent <- dat$id == bad_id & dat$time >= bad_time & !is.na(dat$D)
    dat$D[subsequent] <- 1L
    # Also make sure Y stays absorbing (Y=1 at and after the event)
    # but we need at least one row where both Y=1 and D=1
  }
  expect_error(
    longy_data(
      data = dat, id = "id", time = "time", outcome = "Y",
      treatment = "A", censoring = "C",
      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
      outcome_type = "survival", competing = "D", verbose = FALSE
    ),
    "mutually exclusive"
  )
})

test_that("competing = NULL (default) is backward compatible", {
  dat <- simulate_survival_outcome(n = 50, K = 5, seed = 654)
  obj <- longy_data(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", verbose = FALSE
  )
  expect_null(obj$nodes$competing)
})

test_that("TMLE produces estimates in [0,1] with competing risks", {
  dat <- simulate_competing_risks(n = 200, K = 5, seed = 111)
  result <- longy(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", competing = "D",
    regimes = list(always = 1L),
    estimator = "tmle", learners = NULL,
    inference = "eif", verbose = FALSE
  )
  est <- result$results$always_tmle$estimates
  expect_true(all(est$estimate >= 0 & est$estimate <= 1, na.rm = TRUE))
  # CIF should be monotone non-decreasing
  diffs <- diff(est$estimate)
  expect_true(all(diffs >= -1e-10))  # allow numerical noise
})

test_that("G-comp works with competing risks", {
  dat <- simulate_competing_risks(n = 200, K = 5, seed = 111)
  result <- longy(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", competing = "D",
    regimes = list(always = 1L),
    estimator = "gcomp", learners = NULL,
    n_boot = 0L, verbose = FALSE
  )
  est <- result$results$always_gcomp$estimates
  expect_true(all(est$estimate >= 0 & est$estimate <= 1, na.rm = TRUE))
})

test_that("all-zeros competing column matches standard survival", {
  dat <- simulate_survival_outcome(n = 100, K = 5, seed = 654)
  # Add D column that's always 0
  dat$D <- 0L

  # With competing
  res_comp <- longy(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", competing = "D",
    regimes = list(always = 1L),
    estimator = "gcomp", learners = NULL,
    n_boot = 0L, verbose = FALSE
  )

  # Without competing
  res_std <- longy(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival",
    regimes = list(always = 1L),
    estimator = "gcomp", learners = NULL,
    n_boot = 0L, verbose = FALSE
  )

  # Estimates should be identical
  expect_equal(res_comp$results$always_gcomp$estimates$estimate,
               res_std$results$always_gcomp$estimates$estimate, tolerance = 1e-10)
})

test_that("print.longy_data shows competing column", {
  dat <- simulate_competing_risks(n = 50, K = 5, seed = 111)
  obj <- longy_data(
    data = dat, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    outcome_type = "survival", competing = "D", verbose = FALSE
  )
  out <- capture.output(print(obj))
  expect_true(any(grepl("Competing", out)))
})
