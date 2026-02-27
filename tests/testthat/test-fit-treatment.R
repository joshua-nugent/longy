test_that("fit_treatment runs without error on basic data", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  expect_false(is.null(obj$fits$treatment[["always"]]))
  expect_true(nrow(obj$fits$treatment[["always"]]$predictions) > 0)
})

test_that("fit_treatment predictions are bounded", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  bounds <- c(0.01, 0.99)
  obj <- fit_treatment(obj, regime = "always", bounds = bounds, verbose = FALSE)

  preds <- obj$fits$treatment[["always"]]$predictions$.p_a
  expect_true(all(preds >= bounds[1]))
  expect_true(all(preds <= bounds[2]))
})

test_that("fit_treatment uses correct risk set", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  preds <- obj$fits$treatment[["always"]]$predictions
  # At time 0, everyone should be at risk
  t0 <- preds[preds$.time == 0, ]
  expect_equal(t0$.n_risk[1], sum(d$time == 0))

  # At later times, risk set should be smaller (regime-consistent + uncensored)
  t1 <- preds[preds$.time == 1, ]
  if (nrow(t1) > 0) {
    expect_true(t1$.n_risk[1] <= t0$.n_risk[1])
  }
})

test_that("fit_treatment falls back to marginal with constant treatment", {
  d <- simulate_test_data(n = 50, K = 3)
  # Make treatment constant at time 0
  d$A[d$time == 0] <- 1L
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  t0 <- obj$fits$treatment[["always"]]$predictions[obj$fits$treatment[["always"]]$predictions$.time == 0, ]
  # When all Y=1, marginal method should be used
  expect_equal(unique(t0$.method), "marginal")
})

test_that("fit_treatment works with SuperLearner library", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always",
                       learners = c("SL.glm", "SL.mean"),
                       verbose = FALSE)

  preds <- obj$fits$treatment[["always"]]$predictions
  expect_true(nrow(preds) > 0)

  # Should use SuperLearner, not marginal
  methods <- unique(preds$.method)
  expect_true(any(methods %in% c("SuperLearner", "glm")))

  # sl_info should contain risk and coef
  sl_info <- obj$fits$treatment[["always"]]$sl_info
  expect_true(length(sl_info) > 0)
  sl_entry <- sl_info[[1]]
  if (sl_entry$method %in% c("SuperLearner", "glm")) {
    expect_false(is.null(sl_entry$sl_risk))
    expect_false(is.null(sl_entry$sl_coef))
  }
})

test_that("fit_treatment errors on missing regime", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(
    fit_treatment(obj, regime = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("fit_* functions accept named learner lists", {
  skip_if_not_installed("SuperLearner")
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  # Named list (longy() per-model format)
  lib <- list(treatment = c("SL.glm", "SL.mean"),
              censoring = c("SL.glm", "SL.mean"),
              observation = c("SL.glm", "SL.mean"),
              outcome = c("SL.glm", "SL.mean"))

  obj <- fit_treatment(obj, regime = "always", learners = lib, verbose = FALSE)
  expect_false(is.null(obj$fits$treatment[["always"]]))

  obj <- fit_censoring(obj, regime = "always", learners = lib, verbose = FALSE)
  expect_true(length(obj$fits$censoring[["always"]]) > 0)

  obj <- fit_observation(obj, regime = "always", learners = lib, verbose = FALSE)

  obj <- fit_outcome(obj, regime = "always", learners = lib, verbose = FALSE)
  expect_false(is.null(obj$fits$outcome[["always"]]))
})

test_that(".resolve_learners extracts correct model", {
  lib <- list(treatment = c("SL.glm"), outcome = c("SL.ranger", "SL.mean"))
  expect_equal(.resolve_learners(lib, "treatment"), "SL.glm")
  expect_equal(.resolve_learners(lib, "outcome"), c("SL.ranger", "SL.mean"))
  # Missing model falls back to default element
  lib2 <- list(default = c("SL.mean"), treatment = c("SL.glm"))
  expect_equal(.resolve_learners(lib2, "outcome"), "SL.mean")
  # No default falls back to SL.glm + SL.mean
  expect_equal(.resolve_learners(lib, "censoring"), c("SL.glm", "SL.mean"))
  # NULL passes through
  expect_null(.resolve_learners(NULL, "treatment"))
  # Character vector passes through
  expect_equal(.resolve_learners(c("SL.glm"), "treatment"), "SL.glm")
})

test_that("fit_treatment uses marginal for rare events with large dataset", {
  # 1000 subjects, treatment rate ~0.5% at time 0 (5 events out of 1000)
  set.seed(99)
  n <- 1000
  d <- data.frame(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, n),
    W1 = rep(rnorm(n), each = 2),
    A = 0L,
    R = 1L,
    Y = rbinom(2 * n, 1, 0.3)
  )
  # Make treatment very rare at time 0: only 5 treated
  d$A[d$time == 0] <- c(rep(1L, 5), rep(0L, n - 5))
  # time 1: normal treatment rate
  d$A[d$time == 1] <- rbinom(n, 1, 0.5)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = "W1", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  expect_warning(
    obj <- fit_treatment(obj, regime = "always", min_events = 20L,
                         verbose = FALSE),
    "minority-class events"
  )
  t0 <- obj$fits$treatment[["always"]]$predictions[obj$fits$treatment[["always"]]$predictions$.time == 0, ]
  expect_equal(unique(t0$.method), "marginal")
})

test_that("fit_treatment still fits SL when rate > 1% even with few events", {

  # 50 subjects, 15 events = 30% rate -> should NOT trigger rare-events fallback
  set.seed(100)
  n <- 80
  d <- data.frame(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, n),
    W1 = rep(rnorm(n), each = 2),
    A = rbinom(2 * n, 1, 0.3),
    R = 1L,
    Y = rbinom(2 * n, 1, 0.3)
  )

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = "W1", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  obj <- fit_treatment(obj, regime = "always", min_obs = 10L,
                       min_events = 20L, verbose = FALSE)
  t0 <- obj$fits$treatment[["always"]]$predictions[obj$fits$treatment[["always"]]$predictions$.time == 0, ]
  # With 30% event rate and enough obs, should not be marginal
  expect_true(unique(t0$.method) != "marginal")
})

test_that("fit_treatment respects custom min_events value", {
  # With min_events = 3, even very rare events should pass
  set.seed(101)
  n <- 500
  d <- data.frame(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, n),
    W1 = rep(rnorm(n), each = 2),
    A = 0L,
    R = 1L,
    Y = rbinom(2 * n, 1, 0.3)
  )
  # 5 events out of 500 = 1% rate
  d$A[d$time == 0] <- c(rep(1L, 5), rep(0L, n - 5))
  d$A[d$time == 1] <- rbinom(n, 1, 0.5)

  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = "W1", verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)

  # With min_events = 3, 5 events should be enough -> SL used
  obj <- fit_treatment(obj, regime = "always", min_events = 3L,
                       verbose = FALSE)
  t0 <- obj$fits$treatment[["always"]]$predictions[obj$fits$treatment[["always"]]$predictions$.time == 0, ]
  expect_true(unique(t0$.method) != "marginal")
})

# ---------------------------------------------------------------------------
# risk_set = "followers" tests
# ---------------------------------------------------------------------------

test_that("fit_treatment with risk_set='followers' produces valid fits", {
  set.seed(42)
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, regime = c("always", "never"),
                       risk_set = "followers", verbose = FALSE)

  expect_false(is.null(obj$fits$treatment[["always"]]))
  expect_false(is.null(obj$fits$treatment[["never"]]))
  expect_true(nrow(obj$fits$treatment[["always"]]$predictions) > 0)
  expect_true(nrow(obj$fits$treatment[["never"]]$predictions) > 0)
  # risk_set stored in metadata
  expect_equal(obj$fits$treatment[["always"]]$risk_set, "followers")
  expect_equal(obj$fits$treatment[["never"]]$risk_set, "followers")
})

test_that("risk_set='followers' has smaller or equal risk set vs 'all' at t > 0", {
  set.seed(42)
  d <- simulate_test_data(n = 200, K = 4)

  obj_all <- longy_data(d, id = "id", time = "time", outcome = "Y",
                        treatment = "A", censoring = "C", observation = "R",
                        baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                        verbose = FALSE)
  obj_all <- define_regime(obj_all, "always", static = 1L)
  obj_all <- fit_treatment(obj_all, regime = "always", risk_set = "all",
                           verbose = FALSE)

  obj_fol <- longy_data(d, id = "id", time = "time", outcome = "Y",
                        treatment = "A", censoring = "C", observation = "R",
                        baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                        verbose = FALSE)
  obj_fol <- define_regime(obj_fol, "always", static = 1L)
  obj_fol <- fit_treatment(obj_fol, regime = "always", risk_set = "followers",
                           verbose = FALSE)

  preds_all <- obj_all$fits$treatment[["always"]]$predictions
  preds_fol <- obj_fol$fits$treatment[["always"]]$predictions

  # At t=0, everyone is a follower (no prior treatment), so risk sets equal
  n_all_t0 <- preds_all[preds_all$.time == 0, ]$.n_risk[1]
  n_fol_t0 <- preds_fol[preds_fol$.time == 0, ]$.n_risk[1]
  expect_equal(n_fol_t0, n_all_t0)

  # At t > 0, followers risk set should be <= all risk set
  later_times <- unique(preds_all$.time[preds_all$.time > 0])
  for (tt in later_times) {
    n_all_t <- preds_all[preds_all$.time == tt, ]$.n_risk[1]
    n_fol_t <- preds_fol[preds_fol$.time == tt, ]$.n_risk[1]
    if (!is.null(n_all_t) && !is.null(n_fol_t)) {
      expect_true(n_fol_t <= n_all_t,
                  info = sprintf("Followers risk set at t=%d should be <= all", tt))
    }
  }
})

test_that("risk_set='followers' produces different fits per regime", {
  set.seed(123)
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, regime = c("always", "never"),
                       risk_set = "followers", verbose = FALSE)

  # With followers, "always" and "never" should have different risk set sizes
  # at later times (different subsets follow each regime)
  preds_always <- obj$fits$treatment[["always"]]$predictions
  preds_never <- obj$fits$treatment[["never"]]$predictions

  later_times <- unique(preds_always$.time[preds_always$.time > 0])
  any_differ <- FALSE
  for (tt in later_times) {
    n_always <- preds_always[preds_always$.time == tt, ]$.n_risk[1]
    n_never <- preds_never[preds_never$.time == tt, ]$.n_risk[1]
    if (!is.null(n_always) && !is.null(n_never) && n_always != n_never) {
      any_differ <- TRUE
      break
    }
  }
  # With random binary treatment, followers of "always" vs "never" will differ
  expect_true(any_differ)
})

test_that("default risk_set is 'all' and stored in metadata", {
  d <- simulate_test_data(n = 80, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", baseline = c("W1", "W2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  expect_equal(obj$fits$treatment[["always"]]$risk_set, "all")
})

test_that("risk_set='followers' works in full IPW pipeline", {
  skip_if_not_installed("SuperLearner")
  set.seed(456)
  d <- simulate_test_data(n = 200, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- define_regime(obj, "never", static = 0L)

  obj <- fit_treatment(obj, regime = c("always", "never"),
                       risk_set = "followers", verbose = FALSE)
  obj <- fit_censoring(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- fit_observation(obj, regime = c("always", "never"), verbose = FALSE)
  obj <- compute_weights(obj, regime = c("always", "never"))
  obj <- estimate_ipw(obj, regime = c("always", "never"))

  expect_false(is.null(obj$results$always_ipw))
  expect_false(is.null(obj$results$never_ipw))
  expect_true(nrow(obj$results$always_ipw$estimates) > 0)
  expect_true(nrow(obj$results$never_ipw$estimates) > 0)
})

test_that("longy() forwards risk_set_treatment parameter", {
  skip_if_not_installed("SuperLearner")
  set.seed(789)
  d <- simulate_test_data(n = 150, K = 3)
  obj <- longy(d, id = "id", time = "time", outcome = "Y",
               treatment = "A", censoring = "C", observation = "R",
               baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
               regimes = list(always = 1L, never = 0L),
               risk_set_treatment = "followers", verbose = FALSE)

  expect_equal(obj$fits$treatment[["always"]]$risk_set, "followers")
  expect_equal(obj$fits$treatment[["never"]]$risk_set, "followers")
  expect_false(is.null(obj$results$always_ipw))
})
