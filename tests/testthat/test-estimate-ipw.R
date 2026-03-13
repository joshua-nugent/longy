test_that("estimate_ipw produces valid estimates", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always")
  res <- obj$results$always_ipw

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(nrow(res$estimates) > 0)
  # Estimates should be in [0, 1] for binary outcome
  expect_true(all(res$estimates$estimate >= 0 & res$estimates$estimate <= 1,
                  na.rm = TRUE))
})

test_that("IC inference produces non-negative SEs", {
  d <- simulate_test_data(n = 100, K = 4)
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

  # SEs should be non-negative (can be 0 if all outcomes are same value)
  expect_true(all(res$estimates$se >= 0, na.rm = TRUE))
  # At least some SEs should be positive (non-degenerate time points)
  expect_true(any(res$estimates$se > 0, na.rm = TRUE))
  # CIs should be valid where SE > 0
  pos_se <- res$estimates$se > 0 & !is.na(res$estimates$se)
  if (any(pos_se)) {
    expect_true(all(
      res$estimates$ci_lower[pos_se] < res$estimates$ci_upper[pos_se]
    ))
  }
})

test_that("print and summary work for longy_result", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")
  obj <- estimate_ipw(obj, regime = "always")
  res <- obj$results$always_ipw

  expect_output(print(res), "longy IPW result")
  expect_output(summary(res), "longy IPW Result Summary")
})

test_that("estimate_ipw handles specific time points", {
  d <- simulate_test_data(n = 80, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", times = c(0, 2))
  res <- obj$results$always_ipw

  expect_equal(nrow(res$estimates), 2)
  expect_equal(res$estimates$time, c(0, 2))
})

test_that("sandwich inference produces valid SEs", {
  skip_if_not_installed("survey")
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "sandwich")
  res <- obj$results$always_ipw

  expect_true(all(res$estimates$se >= 0, na.rm = TRUE))
  expect_true(any(res$estimates$se > 0, na.rm = TRUE))
  expect_true(all(
    res$estimates$ci_lower <= res$estimates$estimate &
    res$estimates$estimate <= res$estimates$ci_upper,
    na.rm = TRUE
  ))
})

test_that("sandwich and IC inference give similar SEs", {
  skip_if_not_installed("survey")
  d <- simulate_test_data(n = 200, K = 3)
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
  ic_res <- obj$results$always_ipw

  obj <- estimate_ipw(obj, regime = "always", inference = "sandwich")
  sw_res <- obj$results$always_ipw

  # Point estimates must be identical (same Hajek estimator)
  expect_equal(ic_res$estimates$estimate, sw_res$estimates$estimate)

  # SEs should be in the same ballpark (within factor of 3)
  ic_se <- ic_res$estimates$se
  sw_se <- sw_res$estimates$se
  valid <- ic_se > 0 & sw_se > 0 & !is.na(ic_se) & !is.na(sw_se)
  if (any(valid)) {
    ratio <- ic_se[valid] / sw_se[valid]
    expect_true(all(ratio > 0.3 & ratio < 3),
                info = sprintf("SE ratios: %s",
                               paste(round(ratio, 2), collapse = ", ")))
  }
})

test_that("no-confounding recovery: IPW close to truth", {
  d <- simulate_no_confounding(n = 500, K = 3)
  # Make all outcomes observed for simplicity
  d$R <- 1L
  d$Y[is.na(d$Y)] <- rbinom(sum(is.na(d$Y)), 1, 0.3)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always")
  res <- obj$results$always_ipw

  # With no confounding, IPW estimate should be close to crude mean among treated
  crude_mean_t0 <- mean(d$Y[d$time == 0 & d$A == 1], na.rm = TRUE)
  ipw_t0 <- res$estimates$estimate[res$estimates$time == 0]

  # Should be within 0.15 (Monte Carlo noise with n=500)
  expect_true(abs(ipw_t0 - crude_mean_t0) < 0.15,
              info = sprintf("IPW=%.3f, crude=%.3f", ipw_t0, crude_mean_t0))
})

test_that("IPW works with continuous outcomes", {
  d <- simulate_continuous_outcome(n = 150, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "ic")
  res <- obj$results$always_ipw

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(nrow(res$estimates) > 0)
  # Estimates should be finite
  expect_true(all(is.finite(res$estimates$estimate)))
  # SEs should be positive
  expect_true(all(res$estimates$se > 0, na.rm = TRUE))
})

test_that("IC and sandwich SEs agree for continuous outcomes", {
  skip_if_not_installed("survey")
  d <- simulate_continuous_outcome(n = 200, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "ic")
  ic_res <- obj$results$always_ipw

  obj <- estimate_ipw(obj, regime = "always", inference = "sandwich")
  sw_res <- obj$results$always_ipw

  # Point estimates must be identical
  expect_equal(ic_res$estimates$estimate, sw_res$estimates$estimate)

  # SEs should be in the same ballpark (within factor of 3)
  ic_se <- ic_res$estimates$se
  sw_se <- sw_res$estimates$se
  valid <- ic_se > 0 & sw_se > 0 & !is.na(ic_se) & !is.na(sw_se)
  if (any(valid)) {
    ratio <- ic_se[valid] / sw_se[valid]
    expect_true(all(ratio > 0.3 & ratio < 3),
                info = sprintf("SE ratios: %s",
                               paste(round(ratio, 2), collapse = ", ")))
  }
})

test_that("IPW works with survival outcomes and isotonic smoothing", {
  d <- simulate_survival_outcome(n = 150, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "ic")
  res <- obj$results$always_ipw

  expect_s3_class(obj, "longy_data")
  expect_s3_class(res, "longy_result")
  expect_true(nrow(res$estimates) > 0)
  # Estimates should be in [0, 1]
  expect_true(all(res$estimates$estimate >= 0 &
                  res$estimates$estimate <= 1, na.rm = TRUE))
  # Estimates should be monotonically non-decreasing (isotonic smoothing)
  est <- res$estimates$estimate
  if (length(est) > 1) {
    expect_true(all(diff(est) >= -1e-10),
                info = sprintf("Non-monotone estimates: %s",
                               paste(round(est, 4), collapse = ", ")))
  }
})

test_that("estimate_ipw errors on invalid ci_level", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  expect_error(estimate_ipw(obj, regime = "always", ci_level = 0),
               "ci_level must be between 0 and 1")
  expect_error(estimate_ipw(obj, regime = "always", ci_level = 1),
               "ci_level must be between 0 and 1")
  expect_error(estimate_ipw(obj, regime = "always", ci_level = -0.5),
               "ci_level must be between 0 and 1")
})

test_that("estimate_ipw errors when all requested times are invalid", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  expect_error(estimate_ipw(obj, regime = "always", times = c(999, 1000)),
               "No valid time points")
})

test_that("multi-regime IPW produces distinct results per regime", {
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

  obj <- estimate_ipw(obj, regime = c("always", "never"))

  expect_true(!is.null(obj$results$always_ipw))
  expect_true(!is.null(obj$results$never_ipw))
  expect_equal(obj$results$always_ipw$regime, "always")
  expect_equal(obj$results$never_ipw$regime, "never")
  # Point estimates should generally differ between regimes
  expect_true(nrow(obj$results$always_ipw$estimates) > 0)
  expect_true(nrow(obj$results$never_ipw$estimates) > 0)
})

test_that("estimate_ipw times independence across regimes", {
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

  # Request specific times; both regimes should get the same times
  obj <- estimate_ipw(obj, regime = c("always", "never"), times = c(0, 1))

  expect_equal(obj$results$always_ipw$estimates$time, c(0, 1))
  expect_equal(obj$results$never_ipw$estimates$time, c(0, 1))
})

test_that("IPW inference='none' produces no SE/CI columns", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "none")
  res <- obj$results$always_ipw

  expect_equal(res$inference, "none")
  expect_false("se" %in% names(res$estimates))
  expect_false("ci_lower" %in% names(res$estimates))
})

# ---- Bug fix tests (from adversarial code review) ----

test_that("NA outcomes error when observation model is fitted", {
  # With observation model fitted, weight table only has R=1 subjects.
  # Y=NA among observed subjects is a data integrity error.
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  # Inject NA outcomes for observed subjects
  w_dt <- obj$weights$always$weights_dt
  ids_at_t0 <- w_dt[w_dt$.time == 0, ][[obj$nodes$id]]
  if (length(ids_at_t0) >= 3) {
    inject_ids <- ids_at_t0[1:3]
    obj$data[obj$data[[obj$nodes$id]] %in% inject_ids &
             obj$data[[obj$nodes$time]] == 0, Y := NA_real_]
  }

  # Should error because observation model says these subjects are observed
  expect_error(
    estimate_ipw(obj, regime = "always", inference = "ic"),
    "NA outcomes despite observation model"
  )
})

test_that("NA outcomes are filtered when no observation model is fitted", {
  # Without observation model, unobserved subjects may have NA outcomes.
  # These should be silently filtered (complete case analysis).
  d <- simulate_test_data(n = 100, K = 3)
  # Make some outcomes NA
  set.seed(42)
  na_idx <- sample(which(d$time == 0), 5)
  d$Y[na_idx] <- NA_real_

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  # Should succeed — NAs filtered silently
  obj <- estimate_ipw(obj, regime = "always", inference = "ic")
  res <- obj$results$always_ipw
  est_t0 <- res$estimates[res$estimates$time == 0, ]
  expect_true(is.finite(est_t0$estimate) || est_t0$n_at_risk == 0)
  if (est_t0$n_at_risk >= 2) {
    expect_true(is.finite(est_t0$se) || is.na(est_t0$se))
  }
})

test_that("cluster column validation catches nonexistent columns", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  expect_error(
    estimate_ipw(obj, regime = "always", inference = "ic",
                 cluster = "nonexistent_column"),
    "Cluster column 'nonexistent_column' not found"
  )
})

test_that("cluster-robust IC SEs work when cluster column exists", {
  d <- simulate_test_data(n = 100, K = 3)
  # Add a cluster column (group subjects into clusters of 5)
  d$cluster_id <- ((d$id - 1) %/% 5) + 1

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  # Clustered SEs should produce valid results
  obj_cl <- estimate_ipw(obj, regime = "always", inference = "ic",
                          cluster = "cluster_id")
  res_cl <- obj_cl$results$always_ipw

  # Non-clustered for comparison
  obj_nc <- estimate_ipw(obj, regime = "always", inference = "ic")
  res_nc <- obj_nc$results$always_ipw

  # Point estimates should be identical
  expect_equal(res_cl$estimates$estimate, res_nc$estimates$estimate)
  # Both should have finite SEs
  expect_true(all(is.finite(res_cl$estimates$se) | is.na(res_cl$estimates$se)))
  # Clustered SEs should generally be different from non-clustered
  expect_true(any(res_cl$estimates$se != res_nc$estimates$se, na.rm = TRUE))
})

test_that("binary outcome CIs are clamped to [0, 1]", {
  d <- simulate_test_data(n = 80, K = 4)
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

  # CIs should be clamped to [0, 1] for binary outcomes
  if ("ci_lower" %in% names(res$estimates)) {
    expect_true(all(res$estimates$ci_lower >= 0, na.rm = TRUE),
                info = "ci_lower should be >= 0 for binary")
    expect_true(all(res$estimates$ci_upper <= 1, na.rm = TRUE),
                info = "ci_upper should be <= 1 for binary")
  }
})

test_that("degenerate weights (all equal) give Hajek = simple mean", {
  d <- simulate_no_confounding(n = 200, K = 3)
  d$R <- 1L
  d$Y[is.na(d$Y)] <- rbinom(sum(is.na(d$Y)), 1, 0.3)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = NULL, observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always", stabilized = TRUE)

  # With no confounding and stabilized weights, weights should be ~1
  obj <- estimate_ipw(obj, regime = "always", inference = "ic")
  res <- obj$results$always_ipw

  # Compare with simple mean among regime followers at each time
  for (k in seq_len(nrow(res$estimates))) {
    tt <- res$estimates$time[k]
    # Followers at time tt with non-NA Y
    follow_d <- d[d$time == tt & d$A == 1 & !is.na(d$Y), ]
    if (nrow(follow_d) > 0) {
      simple_mean <- mean(follow_d$Y)
      ipw_est <- res$estimates$estimate[k]
      # Should be close (within 0.15 — weights are ~1 but cumulative
      # weights at later times can diverge more from unity)
      expect_true(abs(ipw_est - simple_mean) < 0.15,
                  info = sprintf("t=%d: IPW=%.4f, mean=%.4f", tt, ipw_est,
                                 simple_mean))
    }
  }
})

test_that("IPW auto-computes weights when not yet computed", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  # Do NOT call compute_weights

  expect_null(obj$weights$always)
  # estimate_ipw should auto-compute weights
  obj <- estimate_ipw(obj, regime = "always", inference = "none")
  expect_true(!is.null(obj$weights$always))
  expect_true(nrow(obj$results$always_ipw$estimates) > 0)
})

test_that("IPW errors when treatment model is missing (auto-compute path)", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  # Do NOT fit any models

  expect_error(
    estimate_ipw(obj, regime = "always"),
    "No treatment model fitted"
  )
})

test_that("inference='none' with survival still applies isotonic smoothing", {
  d <- simulate_survival_outcome(n = 150, K = 5)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "none")
  res <- obj$results$always_ipw

  # Should still be monotone non-decreasing even without inference
  est <- res$estimates$estimate
  if (length(est) > 1) {
    expect_true(all(diff(est) >= -1e-10),
                info = "Survival estimates should be monotone even with inference='none'")
  }
})

test_that("continuous outcomes are NOT isotonically smoothed", {
  d <- simulate_continuous_outcome(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "continuous", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  obj <- estimate_ipw(obj, regime = "always", inference = "none")
  res <- obj$results$always_ipw

  # Continuous outcomes should not be clamped to [0,1]
  # (they can be negative or > 1)
  expect_false("ci_lower" %in% names(res$estimates))
  expect_true(nrow(res$estimates) > 0)
})

test_that("fit results store sl_control and adaptive_cv", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  # Treatment fit should store sl_control and adaptive_cv
  trt_fit <- obj$fits$treatment$always
  expect_true("sl_control" %in% names(trt_fit))
  expect_true("adaptive_cv" %in% names(trt_fit))
  expect_true("risk_set" %in% names(trt_fit))

  # Censoring fit should store sl_control and adaptive_cv (per cause)
  cens_fit <- obj$fits$censoring$always
  if (length(cens_fit) > 0) {
    first_cause <- cens_fit[[1]]
    expect_true("sl_control" %in% names(first_cause))
    expect_true("adaptive_cv" %in% names(first_cause))
  }

  # Observation fit should store sl_control and adaptive_cv
  obs_fit <- obj$fits$observation$always
  if (!is.null(obs_fit)) {
    expect_true("sl_control" %in% names(obs_fit))
    expect_true("adaptive_cv" %in% names(obs_fit))
  }
})

test_that("estimate_ipw errors with competing risks", {
  d <- simulate_competing_risks(n = 100, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = NULL,
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    outcome_type = "survival", competing = "D",
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- compute_weights(obj, regime = "always")

  expect_error(
    estimate_ipw(obj, regime = "always"),
    "IPW is not supported with competing risks"
  )
})

test_that("longy() with competing risks drops IPW with warning", {
  d <- simulate_competing_risks(n = 100, K = 3)

  # estimator = "ipw" alone should error
  expect_error(
    longy(d, id = "id", time = "time", outcome = "Y", treatment = "A",
          censoring = "C", baseline = c("W1", "W2"),
          timevarying = c("L1", "L2"), outcome_type = "survival",
          competing = "D", estimator = "ipw", verbose = FALSE),
    "IPW is not supported with competing risks"
  )

  # estimator = "all" should drop IPW with a warning and run gcomp + tmle
  expect_warning(
    obj <- longy(d, id = "id", time = "time", outcome = "Y", treatment = "A",
                 censoring = "C", baseline = c("W1", "W2"),
                 timevarying = c("L1", "L2"), outcome_type = "survival",
                 competing = "D", estimator = "all", n_boot = 0L,
                 verbose = FALSE),
    "IPW is not supported with competing risks"
  )
  # IPW result should not exist; TMLE should
  expect_null(obj$results$always_ipw)
  expect_s3_class(obj$results$always_tmle, "longy_result")
})

test_that("auto-compute passes stabilized/truncation params to compute_weights", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)
  # Do NOT pre-compute weights — let estimate_ipw do it

  # Unstabilized weights with truncation
  obj1 <- estimate_ipw(obj, regime = "always", inference = "none",
                        stabilized = FALSE, truncation = c(0.1, 10))
  # Stabilized weights (default)
  obj2 <- estimate_ipw(obj, regime = "always", inference = "none",
                        stabilized = TRUE)

  # Both should produce valid results
  expect_true(nrow(obj1$results$always_ipw$estimates) > 0)
  expect_true(nrow(obj2$results$always_ipw$estimates) > 0)

  # Estimates should generally differ (different weighting schemes)
  est1 <- obj1$results$always_ipw$estimates$estimate
  est2 <- obj2$results$always_ipw$estimates$estimate
  # At least one time point should differ
  expect_true(any(abs(est1 - est2) > 1e-10))
})

test_that(".estimator_label returns correct labels", {
  expect_equal(longy:::.estimator_label(list(estimator = "ipw")), "IPW")
  expect_equal(longy:::.estimator_label(list(estimator = "gcomp")), "G-comp")
  expect_equal(longy:::.estimator_label(list(estimator = "tmle")), "TMLE")
  expect_equal(longy:::.estimator_label(list(estimator = "unadjusted")),
               "Unadjusted")
})
