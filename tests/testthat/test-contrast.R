test_that("contrast() produces valid difference contrast with IPW", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr <- contrast(obj, regime = c("always", "never"))

  expect_s3_class(ctr, "longy_contrast")
  expect_equal(ctr$regime, c("always", "never"))
  expect_equal(ctr$estimator, "ipw")
  expect_equal(ctr$scale, "difference")
  expect_equal(ctr$inference, "delta_method")
  expect_true(nrow(ctr$estimates) > 0)
  expect_true(all(c("time", "estimate", "se", "ci_lower", "ci_upper") %in%
                    names(ctr$estimates)))
  # At least some SEs should be non-NA (delta method available)
  expect_true(any(!is.na(ctr$estimates$se)))
})

test_that("contrast() ref= argument works", {
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
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr1 <- contrast(obj, regime = c("always", "never"))
  ctr2 <- contrast(obj, regime = "always", ref = "never")

  expect_equal(ctr1$estimates$estimate, ctr2$estimates$estimate)
  expect_equal(ctr1$estimates$se, ctr2$estimates$se)
})

test_that("contrast() difference equals regime1 - regime2 estimates", {
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
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr <- contrast(obj, regime = c("always", "never"), scale = "difference")

  est_a <- obj$results$always_ipw$estimates
  est_n <- obj$results$never_ipw$estimates
  common_t <- intersect(est_a$time, est_n$time)

  for (tt in common_t) {
    expected <- est_a$estimate[est_a$time == tt] - est_n$estimate[est_n$time == tt]
    actual <- ctr$estimates$estimate[ctr$estimates$time == tt]
    expect_equal(actual, expected, tolerance = 1e-10)
  }
})

test_that("contrast() ratio scale works", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr <- contrast(obj, regime = c("always", "never"), scale = "ratio")

  expect_equal(ctr$scale, "ratio")

  est_a <- obj$results$always_ipw$estimates
  est_n <- obj$results$never_ipw$estimates
  common_t <- intersect(est_a$time, est_n$time)

  for (tt in common_t) {
    expected <- est_a$estimate[est_a$time == tt] / est_n$estimate[est_n$time == tt]
    actual <- ctr$estimates$estimate[ctr$estimates$time == tt]
    expect_equal(actual, expected, tolerance = 1e-10)
  }
})

test_that("contrast() odds_ratio scale works", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr <- contrast(obj, regime = c("always", "never"), scale = "odds_ratio")

  expect_equal(ctr$scale, "odds_ratio")

  est_a <- obj$results$always_ipw$estimates
  est_n <- obj$results$never_ipw$estimates
  common_t <- intersect(est_a$time, est_n$time)

  for (tt in common_t) {
    pa <- est_a$estimate[est_a$time == tt]
    pn <- est_n$estimate[est_n$time == tt]
    expected <- (pa / (1 - pa)) / (pn / (1 - pn))
    actual <- ctr$estimates$estimate[ctr$estimates$time == tt]
    expect_equal(actual, expected, tolerance = 1e-10)
  }
})

test_that("contrast() works with TMLE EIF", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_treatment(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_outcome(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- estimate_tmle(obj, regime = c("always", "never"), inference = "eif",
                       verbose = FALSE)

  ctr <- contrast(obj, regime = c("always", "never"), estimator = "tmle")

  expect_s3_class(ctr, "longy_contrast")
  expect_equal(ctr$estimator, "tmle")
  expect_equal(ctr$inference, "delta_method")
  expect_true(all(!is.na(ctr$estimates$se)))
})

test_that("contrast() auto-detects estimator", {
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
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  # Should auto-detect IPW
  ctr <- contrast(obj, regime = c("always", "never"))
  expect_equal(ctr$estimator, "ipw")
})

test_that("contrast() errors on invalid inputs", {
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
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  # Same regime
  expect_error(contrast(obj, regime = c("always", "always")), "itself")

  # Missing regime
  expect_error(contrast(obj, regime = c("always", "nonexistent")),
               "nonexistent")

  # ref with two regimes
  expect_error(contrast(obj, regime = c("always", "never"), ref = "never"),
               "single regime name")
})

test_that("contrast() G-comp returns point estimates without SEs", {
  d <- simulate_no_confounding(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)
  obj <- fit_outcome(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- estimate_gcomp(obj, regime = c("always", "never"), n_boot = 0,
                        verbose = FALSE)

  expect_message(
    ctr <- contrast(obj, regime = c("always", "never"), estimator = "gcomp"),
    "influence curves"
  )

  expect_equal(ctr$inference, "none")
  expect_true(all(is.na(ctr$estimates$se)))
  # But point estimates should be present
  expect_true(all(!is.na(ctr$estimates$estimate)))
})

test_that("contrast() print/summary/plot methods work", {
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
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"), inference = "ic")

  ctr <- contrast(obj, regime = c("always", "never"))

  expect_output(print(ctr), "Risk Difference")
  expect_output(summary(ctr), "Contrast Summary")
  expect_no_error(plot(ctr))
})

test_that("longy() with contrast=TRUE auto-computes contrasts", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L, never = 0L),
               estimator = "ipw", contrast = TRUE, verbose = FALSE)

  expect_true(length(obj$contrasts) > 0)
  expect_s3_class(obj$contrasts[[1]], "longy_contrast")
})

test_that("IPW result stores per-subject ICs", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")
  obj <- estimate_ipw(obj, regime = "always", inference = "ic")

  res <- obj$results$always_ipw
  expect_false(is.null(res$ic))
  expect_true(nrow(res$ic) > 0)
  expect_true(".ic" %in% names(res$ic))
  expect_true(".time" %in% names(res$ic))

  # IC should have one row per subject per time point
  n_subjects <- obj$meta$n_subjects
  for (tt in unique(res$ic$.time)) {
    ic_t <- res$ic[res$ic$.time == tt, ]
    expect_equal(nrow(ic_t), n_subjects)
  }

  # ICs should sum to approximately 0
  for (tt in unique(res$ic$.time)) {
    ic_t <- res$ic[res$ic$.time == tt, ]
    expect_equal(mean(ic_t$.ic), 0, tolerance = 1e-10)
  }
})

test_that("TMLE result stores per-subject EIF", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  obj <- estimate_tmle(obj, regime = "always", inference = "eif",
                       verbose = FALSE)

  res <- obj$results$always_tmle
  expect_false(is.null(res$ic))
  expect_true(nrow(res$ic) > 0)
  expect_true(".ic" %in% names(res$ic))
  expect_true(".time" %in% names(res$ic))
})
