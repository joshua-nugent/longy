test_that("weight_diagnostics produces correct output", {
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

  diag <- weight_diagnostics(obj, by_time = TRUE)

  expect_true(data.table::is.data.table(diag))
  expect_true("mean_weight" %in% names(diag))
  expect_true("ess" %in% names(diag))
  expect_true(all(diag$ess > 0))
  expect_true(all(diag$mean_weight > 0))
})

test_that("weight_diagnostics overall summary works", {
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

  diag <- weight_diagnostics(obj, by_time = FALSE)

  expect_equal(nrow(diag), 1)
})

test_that("ESS formula is correct", {
  # Kish's ESS = (sum(w))^2 / sum(w^2)
  w <- c(1, 1, 1, 1, 1)
  expect_equal(longy:::.ess(w), 5)

  w2 <- c(2, 2, 2, 2, 2)
  expect_equal(longy:::.ess(w2), 5)

  # Extreme weights should reduce ESS
  w3 <- c(1, 1, 1, 1, 100)
  ess3 <- longy:::.ess(w3)
  expect_true(ess3 < 5)
})

test_that("positivity_diagnostics detects extreme propensities", {
  d <- simulate_test_data(n = 100, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  # With a very strict threshold, should flag more observations
  expect_message(
    diag <- positivity_diagnostics(obj, threshold = 0.5),
    "flagged"
  )
})

test_that("weight_diagnostics errors without weights", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_error(weight_diagnostics(obj), "not computed")
})

# --- Tests for longy_data dispatch (estimate_* returns longy_data) ---

test_that("weight_diagnostics works with longy_data from estimate_ipw", {
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
  obj <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                      inference = "none")

  # longy_data with weights should work
  diag <- weight_diagnostics(obj)
  expect_true(data.table::is.data.table(diag))
  expect_true("ess" %in% names(diag))
})

test_that("weight_diagnostics works with longy_data from longy()", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- suppressWarnings(longy(
    d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "ipw", inference = "none", verbose = FALSE
  ))

  # longy_data from longy() should work
  diag <- weight_diagnostics(obj)
  expect_true(data.table::is.data.table(diag))
  expect_true("ess" %in% names(diag))
})

test_that("positivity_diagnostics works with longy_data from estimate_ipw", {
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
  obj <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                      inference = "none")

  expect_message(
    diag <- positivity_diagnostics(obj, threshold = 0.5),
    "flagged"
  )
})

test_that("diagnostics error on invalid input", {
  expect_error(longy:::.as_longy_data("not an object"),
               "Expected a longy_data")
  expect_error(longy:::.as_longy_data(list(a = 1)),
               "Expected a longy_data")
})

# --- Tests for sl_diagnostics ---

test_that("sl_diagnostics returns correct structure for treatment model", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "treatment")

  expect_true(data.table::is.data.table(diag))
  expect_true(all(c("model", "time", "method") %in% names(diag)))
  expect_true(all(diag$model == "treatment"))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics returns all models", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "all")

  expect_true(data.table::is.data.table(diag))
  models_present <- unique(diag$model)
  expect_true("treatment" %in% models_present)
  expect_true("censoring" %in% models_present)
  expect_true("observation" %in% models_present)
})

test_that("sl_diagnostics shows censoring submodel names", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "censoring")
  expect_true(all(diag$submodel == ".cens_censored"))
})

test_that("sl_diagnostics works with longy_data from estimate_ipw", {
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
  obj <- estimate_ipw(obj, regime = "always", times = c(2, 3),
                      inference = "none")

  diag <- sl_diagnostics(obj)
  expect_true(data.table::is.data.table(diag))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics works with longy_data from longy()", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- suppressWarnings(longy(
    d, id = "id", time = "time", outcome = "Y",
    treatment = "A", censoring = "C", observation = "R",
    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
    regimes = list(always = 1L),
    estimator = "ipw", inference = "none", verbose = FALSE
  ))

  diag <- sl_diagnostics(obj)
  expect_true(data.table::is.data.table(diag))
  expect_true(nrow(diag) > 0)
})

test_that("sl_diagnostics returns empty table when no fits", {
  d <- simulate_test_data(n = 20, K = 3)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", verbose = FALSE)

  expect_message(diag <- sl_diagnostics(obj), "No model fits found")
  expect_equal(nrow(diag), 0)
})

test_that("sl_diagnostics includes target_time column", {
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  diag <- sl_diagnostics(obj, model = "treatment")
  expect_true("target_time" %in% names(diag))
  expect_true(all(is.na(diag$target_time)))
})

# --- Tests for plot_sl_diagnostics ---

test_that("plot_sl_diagnostics messages when no SL fits", {
  skip_if_not_installed("ggplot2")
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)

  # All glm/marginal fits -> no SL coefs to plot
  expect_message(
    p <- plot_sl_diagnostics(obj, type = "weights"),
    "No SuperLearner fits|all glm/marginal"
  )
  expect_null(p)
})

test_that("plot_sl_diagnostics method plot works with glm/marginal fits", {
  skip_if_not_installed("ggplot2")
  d <- simulate_test_data(n = 80, K = 4)
  obj <- longy_data(d, id = "id", time = "time", outcome = "Y",
                    treatment = "A", censoring = "C", observation = "R",
                    baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
                    verbose = FALSE)
  obj <- define_regime(obj, "always", static = 1L)
  obj <- fit_treatment(obj, regime = "always", verbose = FALSE)
  obj <- fit_censoring(obj, regime = "always", verbose = FALSE)
  obj <- fit_observation(obj, regime = "always", verbose = FALSE)

  # Method plot should work even without SL fits
  p <- plot_sl_diagnostics(obj, type = "method")
  expect_true(inherits(p, "gg"))
})

test_that("plot_sl_diagnostics works with synthetic SL data", {
  skip_if_not_installed("ggplot2")

  # Create a synthetic sl_diagnostics data.table with fake SL coefficients
  diag <- data.table::data.table(
    model = c("treatment", "treatment", "treatment",
              "censoring", "censoring",
              "observation", "observation"),
    submodel = c(NA, NA, NA,
                 ".cens_death", ".cens_death",
                 NA, NA),
    time = c(1L, 2L, 3L, 1L, 2L, 1L, 2L),
    target_time = rep(NA_integer_, 7),
    method = rep("SuperLearner", 7),
    sl_risk = list(
      c(SL.glm_All = 0.20, SL.mean_All = 0.25),
      c(SL.glm_All = 0.18, SL.mean_All = 0.24),
      c(SL.glm_All = 0.22, SL.mean_All = 0.26),
      c(SL.glm_All = 0.19, SL.mean_All = 0.23),
      c(SL.glm_All = 0.21, SL.mean_All = 0.25),
      c(SL.glm_All = 0.15, SL.mean_All = 0.20),
      c(SL.glm_All = 0.16, SL.mean_All = 0.21)
    ),
    sl_coef = list(
      c(SL.glm_All = 0.7, SL.mean_All = 0.3),
      c(SL.glm_All = 0.6, SL.mean_All = 0.4),
      c(SL.glm_All = 0.8, SL.mean_All = 0.2),
      c(SL.glm_All = 0.5, SL.mean_All = 0.5),
      c(SL.glm_All = 0.4, SL.mean_All = 0.6),
      c(SL.glm_All = 0.9, SL.mean_All = 0.1),
      c(SL.glm_All = 0.85, SL.mean_All = 0.15)
    )
  )

  # Weights plot
  p <- plot_sl_diagnostics(diag, type = "weights")
  expect_true(inherits(p, "gg"))

  # Risk plot
  p <- plot_sl_diagnostics(diag, type = "risk")
  expect_true(inherits(p, "gg"))

  # Heatmap
  p <- plot_sl_diagnostics(diag, type = "heatmap")
  expect_true(inherits(p, "gg"))

  # Method plot
  p <- plot_sl_diagnostics(diag, type = "method")
  expect_true(inherits(p, "gg"))
})

test_that("plot_sl_diagnostics handles drop_zero correctly", {
  skip_if_not_installed("ggplot2")

  diag <- data.table::data.table(
    model = c("treatment", "treatment"),
    submodel = c(NA, NA),
    time = c(1L, 2L),
    target_time = c(NA_integer_, NA_integer_),
    method = c("SuperLearner", "SuperLearner"),
    sl_risk = list(
      c(SL.glm_All = 0.20, SL.mean_All = 0.25),
      c(SL.glm_All = 0.18, SL.mean_All = 0.24)
    ),
    sl_coef = list(
      c(SL.glm_All = 1.0, SL.mean_All = 0.0),
      c(SL.glm_All = 0.8, SL.mean_All = 0.2)
    )
  )

  # With drop_zero = TRUE, SL.mean_All at t=1 should be excluded from bars
  p <- plot_sl_diagnostics(diag, type = "weights", drop_zero = TRUE)
  expect_true(inherits(p, "gg"))

  # Heatmap with drop_zero keeps learner if non-zero at ANY time
  p <- plot_sl_diagnostics(diag, type = "heatmap", drop_zero = TRUE)
  expect_true(inherits(p, "gg"))
})

test_that("plot_sl_diagnostics errors on invalid input", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_sl_diagnostics("not_an_object"), "must be a longy_data")
})

test_that("plot_sl_diagnostics handles outcome target_time faceting", {
  skip_if_not_installed("ggplot2")

  diag <- data.table::data.table(
    model = c("outcome", "outcome", "outcome", "outcome"),
    submodel = c(NA, NA, NA, NA),
    time = c(1L, 2L, 1L, 2L),
    target_time = c(3L, 3L, 5L, 5L),
    method = rep("SuperLearner", 4),
    sl_risk = list(
      c(SL.glm_All = 0.20), c(SL.glm_All = 0.18),
      c(SL.glm_All = 0.22), c(SL.glm_All = 0.19)
    ),
    sl_coef = list(
      c(SL.glm_All = 1.0), c(SL.glm_All = 1.0),
      c(SL.glm_All = 1.0), c(SL.glm_All = 1.0)
    )
  )

  p <- plot_sl_diagnostics(diag, type = "weights")
  expect_true(inherits(p, "gg"))
  # Should have facets "outcome (t=3)" and "outcome (t=5)"
})
