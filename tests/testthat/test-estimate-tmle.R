# Tests for TMLE estimator

test_that("TMLE produces valid estimates for binary outcome", {
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_equal(result$estimator, "tmle")
  expect_true(nrow(result$estimates) > 0)
  expect_true(all(result$estimates$estimate >= 0 & result$estimates$estimate <= 1,
                  na.rm = TRUE))
})

test_that("TMLE works with continuous outcomes", {
  d <- simulate_continuous_outcome(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_equal(result$estimator, "tmle")
  expect_true(all(is.finite(result$estimates$estimate)))
})

test_that("TMLE works without censoring", {
  d <- simulate_no_censoring(n = 200, K = 3)
  d$Y[is.na(d$Y)] <- 0L

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  expect_s3_class(result, "longy_result")
  expect_true(nrow(result$estimates) > 0)
  expect_true(all(is.finite(result$estimates$estimate)))
})

test_that("TMLE with survival outcomes produces monotone estimates", {
  d <- simulate_survival_outcome(n = 200, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  est <- result$estimates
  expect_true(nrow(est) > 0)
  if (nrow(est) > 1) {
    expect_true(all(diff(est$estimate) >= -1e-10))
  }
  expect_true(all(est$estimate >= 0 & est$estimate <= 1, na.rm = TRUE))
})

test_that("TMLE EIF inference produces valid SEs and CIs", {
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "eif",
                          verbose = FALSE)

  expect_true("se" %in% names(result$estimates))
  expect_true(all(result$estimates$se > 0, na.rm = TRUE))
  expect_true("ci_lower" %in% names(result$estimates))
  expect_true("ci_upper" %in% names(result$estimates))
  # CIs contain the point estimate
  expect_true(all(result$estimates$ci_lower <= result$estimates$estimate, na.rm = TRUE))
  expect_true(all(result$estimates$ci_upper >= result$estimates$estimate, na.rm = TRUE))
})

test_that("TMLE is close to G-comp under well-specified models", {
  d <- simulate_test_data(n = 300, K = 3, seed = 88)
  # G-comp
  obj_g <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_g <- define_regime(obj_g, "always", static = 1L)
  obj_g <- fit_outcome(obj_g, regime = "always", verbose = FALSE)
  gcomp_result <- estimate_gcomp(obj_g, regime = "always", n_boot = 0,
                                  verbose = FALSE)

  # TMLE
  obj_t <- longy_data(d, id = "id", time = "time", outcome = "Y",
                      treatment = "A", censoring = "C", observation = "R",
                      baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                      verbose = FALSE)
  obj_t <- define_regime(obj_t, "always", static = 1L)
  obj_t <- fit_treatment(obj_t, regime = "always", verbose = FALSE)
  obj_t <- fit_censoring(obj_t, regime = "always", verbose = FALSE)
  obj_t <- fit_outcome(obj_t, regime = "always", verbose = FALSE)
  tmle_result <- estimate_tmle(obj_t, regime = "always", inference = "none",
                                verbose = FALSE)

  # TMLE should be reasonably close to G-comp (both valid)
  shared_times <- intersect(gcomp_result$estimates$time,
                            tmle_result$estimates$time)
  for (tt in shared_times) {
    gc_est <- gcomp_result$estimates$estimate[gcomp_result$estimates$time == tt]
    tm_est <- tmle_result$estimates$estimate[tmle_result$estimates$time == tt]
    expect_true(abs(gc_est - tm_est) < 0.2,
                label = sprintf("TMLE vs G-comp at t=%d: gc=%.3f tm=%.3f",
                                tt, gc_est, tm_est))
  }
})

test_that("TMLE no-confounding recovery", {
  d <- simulate_no_confounding(n = 500, K = 3, seed = 777)
  d$Y[is.na(d$Y)] <- 0L

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  # Crude mean among treated at earliest time (no confounding)
  t0 <- min(result$estimates$time)
  crude <- mean(d$Y[d$A == 1 & d$time == t0], na.rm = TRUE)
  tmle_est <- result$estimates$estimate[result$estimates$time == t0]
  # Allow generous tolerance since this is a simulation
  expect_equal(tmle_est, crude, tolerance = 0.15)
})

test_that("estimate_tmle errors without treatment model", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  expect_error(estimate_tmle(obj, regime = "always"),
               "Treatment model not fit")
})

test_that("estimate_tmle errors without outcome model", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  expect_error(estimate_tmle(obj, regime = "always"),
               "Outcome model not fitted")
})

test_that("longy() with estimator = 'tmle' works end-to-end", {
  d <- simulate_test_data(n = 150, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "tmle",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_true("always" %in% names(results))
  expect_s3_class(results$always, "longy_result")
  expect_equal(results$always$estimator, "tmle")
  expect_true(nrow(results$always$estimates) > 0)
  expect_true(all(results$always$estimates$estimate >= 0 &
                  results$always$estimates$estimate <= 1, na.rm = TRUE))
  # EIF inference should be present by default for TMLE
  expect_true("se" %in% names(results$always$estimates))
})

test_that("TMLE print method outputs 'TMLE'", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)
  result <- estimate_tmle(obj, regime = "always", inference = "none",
                          verbose = FALSE)

  expect_output(print(result), "TMLE")
})

test_that(".compute_cumulative_g returns valid results", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  g_dt <- longy:::.compute_cumulative_g(obj, regime = "always")

  expect_true(nrow(g_dt) > 0)
  expect_true(".g_a" %in% names(g_dt))
  expect_true(".g_c" %in% names(g_dt))
  expect_true(".g_cum" %in% names(g_dt))

  # g_cum should be bounded below at 0.01 (default)
  expect_true(all(g_dt$.g_cum >= 0.01))
  # g_cum should be non-increasing per subject
  g_dt_ord <- g_dt[order(g_dt$id, g_dt$.time), ]
  for (uid in unique(g_dt_ord$id)) {
    sub <- g_dt_ord[g_dt_ord$id == uid, ]
    if (nrow(sub) > 1) {
      diffs <- diff(sub$.g_cum)
      # Allow small numerical tolerance
      expect_true(all(diffs <= 1e-10),
                  label = sprintf("g_cum non-increasing for id=%s", uid))
    }
  }
})

test_that("longy() with estimator = 'tmle' and continuous outcome", {
  d <- simulate_continuous_outcome(n = 200, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "tmle",
    outcome_type = "continuous",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_equal(results$always$estimator, "tmle")
  expect_true(all(is.finite(results$always$estimates$estimate)))
})

test_that("longy() with estimator = 'tmle' and never-treat regime", {
  d <- simulate_test_data(n = 150, K = 3)
  results <- longy(
    data = d,
    id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(never = 0L),
    estimator = "tmle",
    n_boot = 0,
    verbose = FALSE
  )

  expect_s3_class(results, "longy_results")
  expect_true("never" %in% names(results))
  expect_equal(results$never$estimator, "tmle")
  expect_true(nrow(results$never$estimates) > 0)
})

test_that("TMLE EIF SEs for continuous outcomes", {
  d <- simulate_continuous_outcome(n = 200, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_outcome(obj, regime = "always", verbose = FALSE)

  result <- estimate_tmle(obj, regime = "always", inference = "eif",
                          verbose = FALSE)

  expect_true("se" %in% names(result$estimates))
  expect_true(all(result$estimates$se > 0, na.rm = TRUE))
  expect_true(all(result$estimates$ci_lower < result$estimates$ci_upper, na.rm = TRUE))
})
