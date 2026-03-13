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

# --- Time-varying static regimes ---

test_that("define_regime accepts named vector static", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  obj <- define_regime(obj, name = "early",
                       static = c("1" = 1, "2" = 1, "3" = 0))

  expect_equal(obj$regimes$early$type, "static")
  expect_true(obj$regimes$early$time_varying)
  expect_equal(obj$regimes$early$value, c("1" = 1L, "2" = 1L, "3" = 0L))
})

test_that("define_regime accepts function static", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  obj <- define_regime(obj, name = "early_fn",
                       static = function(t) as.integer(t <= 2))

  expect_equal(obj$regimes$early_fn$type, "static")
  expect_true(obj$regimes$early_fn$time_varying)
  expect_true(is.function(obj$regimes$early_fn$value))
})

test_that("define_regime rejects bad named vector values", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  # Values not 0/1
  expect_error(
    define_regime(obj, name = "bad", static = c("1" = 2, "2" = 1)),
    "0 or 1"
  )

  # Missing names
  expect_error(
    define_regime(obj, name = "bad2", static = c(1, 0, 1)),
    "named"
  )
})

test_that(".resolve_static_at_time handles scalar", {
  result <- longy:::.resolve_static_at_time(1L, c(1, 2, 3))
  expect_equal(result, c(1L, 1L, 1L))

  result <- longy:::.resolve_static_at_time(0L, c(1, 2, 3))
  expect_equal(result, c(0L, 0L, 0L))
})

test_that(".resolve_static_at_time handles named vector", {
  val <- c("1" = 1L, "2" = 1L, "3" = 0L)
  result <- longy:::.resolve_static_at_time(val, c(1, 2, 3))
  expect_equal(result, c(1L, 1L, 0L))

  # Repeated times (as in real data with multiple subjects)
  result <- longy:::.resolve_static_at_time(val, c(1, 1, 2, 2, 3, 3))
  expect_equal(result, c(1L, 1L, 1L, 1L, 0L, 0L))
})

test_that(".resolve_static_at_time handles function", {
  fn <- function(t) as.integer(t <= 2)
  result <- longy:::.resolve_static_at_time(fn, c(1, 2, 3))
  expect_equal(result, c(1L, 1L, 0L))
})

test_that(".resolve_static_at_time errors on missing time values", {
  val <- c("1" = 1L, "2" = 0L)
  expect_error(
    longy:::.resolve_static_at_time(val, c(1, 2, 3)),
    "no value defined for time.*3"
  )
})

test_that(".resolve_static_at_time errors on non-0/1 from function", {
  fn <- function(t) t  # returns 1, 2, 3 — not 0/1
  expect_error(
    longy:::.resolve_static_at_time(fn, c(1, 2, 3)),
    "must all be 0 or 1"
  )
})

test_that(".evaluate_regime works for time-varying static (vector)", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  # Times are 0, 1, 2
  obj <- define_regime(obj, name = "early",
                       static = c("0" = 1, "1" = 1, "2" = 0))

  vals <- longy:::.evaluate_regime(obj$regimes$early, obj$data,
                                   time_col = "time")
  expect_length(vals, nrow(obj$data))
  expect_true(all(vals[obj$data$time <= 1] == 1L))
  expect_true(all(vals[obj$data$time == 2] == 0L))
})

test_that(".evaluate_regime works for time-varying static (function)", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, name = "early_fn",
                       static = function(t) as.integer(t <= 1))

  vals <- longy:::.evaluate_regime(obj$regimes$early_fn, obj$data,
                                   time_col = "time")
  expect_length(vals, nrow(obj$data))
  expect_true(all(vals[obj$data$time <= 1] == 1L))
  expect_true(all(vals[obj$data$time == 2] == 0L))
})

test_that("time-varying static requires time_col in .evaluate_regime", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, name = "early",
                       static = c("1" = 1, "2" = 1, "3" = 0))

  expect_error(
    longy:::.evaluate_regime(obj$regimes$early, obj$data),
    "time_col required"
  )
})

test_that("time-varying static regime end-to-end with IPW", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  # Define time-varying: treat at times 0-1, don't treat at time 2
  obj <- define_regime(obj, name = "early_treat",
                       static = c("0" = 1, "1" = 1, "2" = 0))
  # Also define constant for comparison
  obj <- define_regime(obj, name = "always", static = 1L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- compute_weights(obj, stabilized = TRUE)

  # Both regimes should have weights
  expect_false(is.null(obj$weights$early_treat))
  expect_false(is.null(obj$weights$always))

  obj <- estimate_ipw(obj, inference = "none")

  # Both should produce estimates
  expect_true("early_treat_ipw" %in% names(obj$results))
  expect_true("always_ipw" %in% names(obj$results))

  # Estimates should be finite
  expect_true(all(is.finite(obj$results$early_treat_ipw$estimates$estimate)))
  expect_true(all(is.finite(obj$results$always_ipw$estimates$estimate)))

  # Estimates should differ (different regimes)
  early_est <- obj$results$early_treat_ipw$estimates$estimate
  always_est <- obj$results$always_ipw$estimates$estimate
  expect_false(all(early_est == always_est))
})

test_that("time-varying static function regime end-to-end with IPW", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  # Function form: treat at times 0-1 only
  obj <- define_regime(obj, name = "early_fn",
                       static = function(t) as.integer(t <= 1))

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- compute_weights(obj, stabilized = TRUE)
  obj <- estimate_ipw(obj, inference = "none")

  expect_true("early_fn_ipw" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$early_fn_ipw$estimates$estimate)))
})

test_that("time-varying static regime via longy() named vector", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)

  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", outcome_type = "continuous",
               baseline = c("W1", "W2"),
               timevarying = c("L1", "L2"),
               regimes = list(
                 always = 1L,
                 early = c("0" = 1, "1" = 1, "2" = 0)
               ),
               estimator = "ipw",
               verbose = FALSE)

  expect_true("early_ipw" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$early_ipw$estimates$estimate)))
})

test_that("time-varying static regime with G-comp", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "early",
                       static = c("0" = 1, "1" = 1, "2" = 0))
  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)
  obj <- estimate_gcomp(obj, n_boot = 0, verbose = FALSE)

  expect_true("early_gcomp" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$early_gcomp$estimates$estimate)))
})

test_that("time-varying static regime with TMLE", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "early",
                       static = c("0" = 1, "1" = 1, "2" = 0))
  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)

  expect_no_error(
    obj <- estimate_tmle(obj, inference = "none", verbose = FALSE)
  )
  expect_true("early_tmle" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$early_tmle$estimates$estimate)))
})

# --- Shifted (pre-computed) regimes ---

test_that("define_regime accepts shifted column", {
  d <- simulate_test_data(n = 20, K = 3)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  obj <- define_regime(obj, name = "custom", shifted = "A_cf")

  expect_equal(obj$regimes$custom$type, "shifted")
  expect_equal(obj$regimes$custom$value, "A_cf")
})

test_that("define_regime rejects missing shifted column", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    define_regime(obj, name = "bad", shifted = "nonexistent"),
    "not found in data"
  )
})

test_that("define_regime rejects non-0/1 shifted column", {
  d <- simulate_test_data(n = 20, K = 3)
  d$A_cf <- d$time  # values 0, 1, 2 — 2 is invalid
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    define_regime(obj, name = "bad", shifted = "A_cf"),
    "non-0/1"
  )
})

test_that("define_regime rejects non-character shifted", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    define_regime(obj, name = "bad", shifted = 42),
    "single column name"
  )
})

test_that(".evaluate_regime works for shifted regimes", {
  d <- simulate_test_data(n = 20, K = 3)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, name = "custom", shifted = "A_cf")

  vals <- longy:::.evaluate_regime(obj$regimes$custom, obj$data)
  expect_length(vals, nrow(obj$data))
  expect_true(all(vals[obj$data$time <= 1] == 1L))
  expect_true(all(vals[obj$data$time == 2] == 0L))
})

test_that("shifted regime end-to-end with IPW", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)
  # Pre-compute: treat at times 0-1, don't treat at time 2
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "shifted_early", shifted = "A_cf")
  obj <- define_regime(obj, name = "always", static = 1L)

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- compute_weights(obj, stabilized = TRUE)
  obj <- estimate_ipw(obj, inference = "none")

  expect_true("shifted_early_ipw" %in% names(obj$results))
  expect_true("always_ipw" %in% names(obj$results))

  # Estimates should be finite
  expect_true(all(is.finite(obj$results$shifted_early_ipw$estimates$estimate)))

  # Shifted and always should produce different estimates
  shifted_est <- obj$results$shifted_early_ipw$estimates$estimate
  always_est <- obj$results$always_ipw$estimates$estimate
  expect_false(all(shifted_est == always_est))
})

test_that("shifted regime via longy() wrapper", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)

  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", outcome_type = "continuous",
               baseline = c("W1", "W2"),
               timevarying = c("L1", "L2"),
               regimes = list(always = 1L, custom = "A_cf"),
               estimator = "ipw",
               verbose = FALSE)

  expect_true("custom_ipw" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$custom_ipw$estimates$estimate)))
})

test_that("shifted regime with G-comp", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "custom", shifted = "A_cf")
  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)
  obj <- estimate_gcomp(obj, n_boot = 0, verbose = FALSE)

  expect_true("custom_gcomp" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$custom_gcomp$estimates$estimate)))
})

test_that("shifted regime with TMLE", {
  set.seed(42)
  n <- 100
  K <- 3
  d <- simulate_test_data(n = n, K = K)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "custom", shifted = "A_cf")
  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)

  expect_no_error(
    obj <- estimate_tmle(obj, inference = "none", verbose = FALSE)
  )
  expect_true("custom_tmle" %in% names(obj$results))
  expect_true(all(is.finite(obj$results$custom_tmle$estimates$estimate)))
})

test_that("shifted regime matches equivalent time-varying static", {
  # A shifted column that is purely a function of time should give the

  # same estimates as the equivalent time-varying static regime.
  set.seed(42)
  n <- 150
  K <- 3
  d <- simulate_test_data(n = n, K = K)
  # Pre-compute column: treat at 0-1, not at 2 (same as static below)
  d$A_cf <- ifelse(d$time <= 1, 1L, 0L)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", outcome_type = "continuous",
                    baseline = c("W1", "W2"),
                    timevarying = c("L1", "L2"), verbose = FALSE)

  obj <- define_regime(obj, name = "shifted_ver", shifted = "A_cf")
  obj <- define_regime(obj, name = "static_ver",
                       static = c("0" = 1, "1" = 1, "2" = 0))

  obj <- fit_treatment(obj, verbose = FALSE)
  obj <- fit_outcome(obj, verbose = FALSE)
  obj <- estimate_gcomp(obj, n_boot = 0, verbose = FALSE)

  est_shifted <- obj$results$shifted_ver_gcomp$estimates$estimate
  est_static  <- obj$results$static_ver_gcomp$estimates$estimate

  expect_equal(est_shifted, est_static, tolerance = 1e-10)
})

# =============================================================================
# add_regime() tests
# =============================================================================

test_that("add_regime adds a regime and produces IPW results", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  # Add never-treat
  obj <- add_regime(obj, name = "never", static = 0L, verbose = FALSE)

  expect_true("never" %in% names(obj$regimes))
  expect_true("never_ipw" %in% names(obj$results))
  expect_true("never_unadjusted" %in% names(obj$results))
  est <- obj$results$never_ipw$estimates
  expect_true(nrow(est) > 0)
  expect_true(all(c("time", "estimate") %in% names(est)))
})

test_that("add_regime reuses nuisance fits from donor regime", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  obj <- add_regime(obj, name = "never", static = 0L, verbose = FALSE)

  # Treatment, censoring, observation fits should exist for "never"
  expect_false(is.null(obj$fits$treatment[["never"]]))
  expect_false(is.null(obj$fits$censoring[["never"]]))
  expect_false(is.null(obj$fits$observation[["never"]]))

  # Predictions should be identical (reused, not re-fitted)
  # Treatment: same predictions (fit once on full population, risk_set="all")
  expect_identical(
    obj$fits$treatment[["always"]]$predictions,
    obj$fits$treatment[["never"]]$predictions
  )
  # Censoring: always identical
  expect_identical(
    obj$fits$censoring[["always"]],
    obj$fits$censoring[["never"]]
  )
})

test_that("add_regime auto-detects estimator from existing results", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "gcomp", n_boot = 0,
               learners = NULL, verbose = FALSE)

  obj <- add_regime(obj, name = "never", static = 0L,
                    n_boot = 0, verbose = FALSE)

  # Should auto-detect gcomp
  expect_true("never_gcomp" %in% names(obj$results))
  expect_false("never_ipw" %in% names(obj$results))
})

test_that("add_regime works with time-varying static regimes", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  obj <- add_regime(obj, name = "early_stop",
                    static = c("0" = 1, "1" = 1, "2" = 0),
                    verbose = FALSE)

  expect_equal(obj$regimes$early_stop$type, "static")
  expect_true(obj$regimes$early_stop$time_varying)
  expect_true("early_stop_ipw" %in% names(obj$results))
})

test_that("add_regime works with shifted regimes", {
  d <- simulate_test_data(n = 100, K = 3)
  d$A_cf <- ifelse(d$time == 0, 1L, 0L)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  obj <- add_regime(obj, name = "shifted", shifted = "A_cf",
                    verbose = FALSE)

  expect_equal(obj$regimes$shifted$type, "shifted")
  expect_true("shifted_ipw" %in% names(obj$results))
})

test_that("add_regime with multiple estimators", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = c("ipw", "tmle"), learners = NULL,
               verbose = FALSE)

  obj <- add_regime(obj, name = "never", static = 0L, verbose = FALSE)

  expect_true("never_ipw" %in% names(obj$results))
  expect_true("never_tmle" %in% names(obj$results))
})

test_that("add_regime with explicit estimator override", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  # Existing results are IPW, but we want gcomp for the new regime
  obj <- add_regime(obj, name = "never", static = 0L,
                    estimator = "gcomp", n_boot = 0, verbose = FALSE)

  expect_true("never_gcomp" %in% names(obj$results))
  expect_false("never_ipw" %in% names(obj$results))
})

test_that("add_regime errors on duplicate regime name", {
  d <- simulate_test_data(n = 50, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  expect_error(
    add_regime(obj, name = "always", static = 0L, verbose = FALSE),
    "already exists"
  )
})

test_that("add_regime works without prior fits (no donor)", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)

  # No prior longy() call — no fits to reuse. Pass learners explicitly
  # since there's no donor to inherit from.
  obj <- add_regime(obj, name = "always", static = 1L,
                    estimator = "ipw", learners = NULL,
                    verbose = FALSE)

  expect_true("always_ipw" %in% names(obj$results))
})

test_that("add_regime errors on learner mismatch with donor", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  # Existing fit used learners = NULL (glm). Passing SL libraries should error.
  expect_error(
    add_regime(obj, name = "never", static = 0L,
               learners = c("SL.glm", "SL.mean"),
               verbose = FALSE),
    "Learner mismatch"
  )
})

test_that("add_regime accepts matching learners without error", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  # Same learners as the original call — should succeed
  obj <- add_regime(obj, name = "never", static = 0L,
                    learners = NULL,
                    verbose = FALSE)
  expect_true("never_ipw" %in% names(obj$results))
})

test_that("add_regime accepts named list learners matching donor", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw",
               learners = c("SL.glm", "SL.mean"),
               verbose = FALSE)

  # Same library — should pass validation
  obj <- add_regime(obj, name = "never", static = 0L,
                    learners = c("SL.glm", "SL.mean"),
                    verbose = FALSE)
  expect_true("never_ipw" %in% names(obj$results))

  # Different library — should error
  expect_error(
    add_regime(obj, name = "early", static = c("0" = 1, "1" = 0, "2" = 0),
               learners = c("SL.glm"),
               verbose = FALSE),
    "Learner mismatch"
  )
})

test_that("add_regime gives same results as running regime through longy()", {
  d <- simulate_test_data(n = 200, K = 3, seed = 42)

  # Run longy() with both regimes together
  obj_both <- longy(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    regimes = list(always = 1L, never = 0L),
                    estimator = "ipw", learners = NULL,
                    verbose = FALSE)

  # Run longy() with one regime, then add_regime
  obj_add <- longy(d, id = "id", time = "time", outcome = "Y",
                   treatment = "A", censoring = "C",
                   baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                   regimes = list(always = 1L),
                   estimator = "ipw", learners = NULL,
                   verbose = FALSE)
  obj_add <- add_regime(obj_add, name = "never", static = 0L,
                        verbose = FALSE)

  # The "never" IPW estimates should be identical because the treatment
  # model was fit once on full population (risk_set="all") and reused
  est_both <- obj_both$results$never_ipw$estimates$estimate
  est_add  <- obj_add$results$never_ipw$estimates$estimate

  expect_equal(est_add, est_both, tolerance = 1e-10)
})

test_that("add_regime can add multiple regimes sequentially", {
  d <- simulate_test_data(n = 100, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L),
               estimator = "ipw", learners = NULL,
               verbose = FALSE)

  obj <- add_regime(obj, name = "never", static = 0L, verbose = FALSE)
  obj <- add_regime(obj, name = "early",
                    static = c("0" = 1, "1" = 0, "2" = 0),
                    verbose = FALSE)

  expect_equal(length(obj$regimes), 3)
  expect_true(all(c("always_ipw", "never_ipw", "early_ipw") %in%
                    names(obj$results)))
})
