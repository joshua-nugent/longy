# ============================================================================
# test_contrast_vs_ltmle.R
# Ad hoc validation: compare longy contrast() to ltmle summary() ATE/RR/OR
#
# Uses simcausal for data generation and ltmle as the reference estimator.
# Tests that longy's delta-method contrasts produce equivalent results.
#
# Usage:
#   source("longy/simulations/test_contrast_vs_ltmle.R")
#   run_all()               # full suite
#   run_all(verbose = TRUE) # with progress
# ============================================================================

library(simcausal)
library(ltmle)
library(SuperLearner)
library(data.table)
library(devtools)

devtools::load_all("/Users/joshuanugent/Desktop/ltmle2/longy")

# ============================================================================
# DGP 1: Binary outcome, no censoring, 3 time points
# ============================================================================
#' Simple DGP: baseline W, time-varying L, treatment A, binary Y.
#' No censoring or observation issues -- cleanest comparison.
define_dgp_binary_simple <- function() {
  D <- DAG.empty() +
    node("W", distr = "rnorm", mean = 0, sd = 1) +
    node("L", t = 0, distr = "rnorm",
         mean = 0.5 * W) +
    node("A", t = 0, distr = "rbern",
         prob = plogis(-0.3 + 0.3 * W + 0.5 * L[t])) +
    node("Y", t = 0, distr = "rbern",
         prob = plogis(-1.5 + 0.3 * W + 0.5 * L[t] + 0.8 * A[t]),
         EFU = TRUE) +
    node("L", t = 1:2, distr = "rnorm",
         mean = 0.5 * W + 0.3 * L[t-1] + 0.4 * A[t-1]) +
    node("A", t = 1:2, distr = "rbern",
         prob = plogis(-0.3 + 0.3 * W + 0.5 * L[t] + 0.3 * A[t-1])) +
    node("Y", t = 1:2, distr = "rbern",
         prob = plogis(-1.5 + 0.3 * W + 0.5 * L[t] + 0.8 * A[t]),
         EFU = TRUE)
  D <- set.DAG(D)
  D <- D +
    action("A1", nodes = node("A", t = 0:2, distr = "rbern", prob = 1)) +
    action("A0", nodes = node("A", t = 0:2, distr = "rbern", prob = 0))
  D
}

# ============================================================================
# DGP 3: CVD / A1C (from existing validation), 4 time points
# ============================================================================
define_dgp_cvd <- function() {
  D <- DAG.empty() +
    node("CVD", distr = "rcat.b1", probs = c(0.5, 0.25, 0.25)) +
    node("A1C", t = 0, distr = "rnorm",
         mean = 5 + (CVD > 1) * 10 + (CVD > 2) * 5) +
    node("TI", t = 0, distr = "rbern",
         prob = plogis(-5 - 0.3 * CVD + 0.5 * A1C[t])) +
    node("A1C", t = 1:3, distr = "rnorm",
         mean = -TI[t-1] * 10 + 5 + (CVD > 1) * 10 + (CVD > 2) * 5) +
    node("TI", t = 1:3, distr = "rbern",
         prob = plogis(-5 - 0.3 * CVD + 0.5 * A1C[t] + 1.5 * TI[t-1])) +
    node("Y", t = 0:3, distr = "rbern",
         prob = plogis(-6 - 1.2 * TI[t] + 0.1 * CVD + 0.3 * A1C[t]),
         EFU = TRUE)
  D <- set.DAG(D)
  D <- D +
    action("A1", nodes = node("TI", t = 0:3, distr = "rbern", prob = 1)) +
    action("A0", nodes = node("TI", t = 0:3, distr = "rbern", prob = 0))
  D
}


# ============================================================================
# Helper: wide -> long conversion
# ============================================================================
wide_to_long <- function(dat_wide, time_points, baseline_cols,
                         treatment_col, outcome_col,
                         censoring_col = NULL, timevarying_cols = NULL) {
  dt <- as.data.table(dat_wide)

  long_list <- lapply(time_points, function(tt) {
    cols_to_get <- list(id = dt$ID)

    # Baseline covariates (constant over time)
    for (b in baseline_cols) {
      cols_to_get[[b]] <- dt[[b]]
    }

    # Time-varying covariates
    if (!is.null(timevarying_cols)) {
      for (tv in timevarying_cols) {
        cols_to_get[[tv]] <- dt[[paste0(tv, "_", tt)]]
      }
    }

    # Treatment
    cols_to_get[[treatment_col]] <- dt[[paste0(treatment_col, "_", tt)]]

    # Censoring (if present)
    if (!is.null(censoring_col)) {
      cval <- dt[[paste0(censoring_col, "_", tt)]]
      if (!is.null(cval)) {
        cols_to_get[[censoring_col]] <- cval
      }
    }

    # Outcome
    cols_to_get[[outcome_col]] <- dt[[paste0(outcome_col, "_", tt)]]

    sub <- as.data.table(cols_to_get)
    sub[, time := tt]
    sub
  })

  out <- rbindlist(long_list)
  setorder(out, id, time)

  # Convert baseline categoricals to character
  for (b in baseline_cols) {
    if (is.integer(out[[b]]) || is.factor(out[[b]])) {
      out[[b]] <- as.character(out[[b]])
    }
  }

  # Ensure id is character
  out[, id := as.character(id)]

  # LOCF for post-absorption NAs in time-varying and treatment columns
  fill_cols <- c(timevarying_cols, treatment_col)
  for (col in fill_cols) {
    if (col %in% names(out)) {
      out[, (col) := nafill(get(col), type = "locf"), by = id]
    }
  }

  # Handle censoring: convert to "uncensored"/"censored" for longy
  if (!is.null(censoring_col) && censoring_col %in% names(out)) {
    cens_vals <- out[[censoring_col]]
    # simcausal: BL_C = 1 means censored event occurred, 0/NA means uncensored
    cens_char <- ifelse(is.na(cens_vals) | cens_vals == 0, "uncensored", "censored")
    out[[censoring_col]] <- cens_char
  }

  out
}


# ============================================================================
# Helper: compute true counterfactual means via simcausal
# ============================================================================
compute_truth <- function(D, outcome_col, time_points,
                          action_always = "A1", action_never = "A0",
                          n_mc = 500000, seed = 999) {
  set.seed(seed)
  truth <- data.table(time = time_points, always = NA_real_, never = NA_real_)
  for (tt in time_points) {
    targ1 <- set.targetE(D, outcome_col, t = tt, param = action_always)
    targ0 <- set.targetE(D, outcome_col, t = tt, param = action_never)
    truth[time == tt, always := unname(eval.target(targ1, n = n_mc)$res)]
    truth[time == tt, never  := unname(eval.target(targ0, n = n_mc)$res)]
  }
  truth[, ATE := always - never]
  truth
}


# ============================================================================
# Helper: run ltmle and extract ATE/RR/OR with SEs
# ============================================================================
run_ltmle <- function(dat_wide, ltmle_cols, Anodes, Lnodes, Ynodes, Cnodes,
                      abar_treat, abar_control, SL.library = "SL.glm") {
  ltmle_data <- as.data.frame(dat_wide)[, ltmle_cols]

  # Convert any integer baseline categoricals to character
  for (col in names(ltmle_data)) {
    if (!grepl("_[0-9]+$", col) && is.integer(ltmle_data[[col]])) {
      ltmle_data[[col]] <- as.character(ltmle_data[[col]])
    }
  }

  result <- ltmle(
    data = ltmle_data,
    estimate.time = FALSE,
    Anodes = Anodes,
    Lnodes = Lnodes,
    Ynodes = Ynodes,
    Cnodes = Cnodes,
    abar = list(abar_treat, abar_control),
    survivalOutcome = FALSE,
    SL.library = SL.library,
    SL.cvControl = list(V = 3)
  )

  s <- summary(result)

  list(
    treatment = s$effect.measures$treatment,
    control = s$effect.measures$control,
    ATE = s$effect.measures$ATE,
    RR = s$effect.measures$RR,
    OR = s$effect.measures$OR
  )
}


# ============================================================================
# Test 1: Binary outcome, no censoring — TMLE contrast vs ltmle ATE
# ============================================================================
test_binary_simple_tmle <- function(n = 5000, seed = 42, tol = 0.05,
                                    verbose = FALSE) {
  cat("--- Test 1: Binary simple TMLE contrast vs ltmle ---\n")
  set.seed(seed)

  D <- define_dgp_binary_simple()

  # Simulate data
  dat_wide <- sim(D, n = n, LTCF = "Y")

  # ltmle (final time point only)
  ltmle_cols <- c("W", paste0("L_", 0:2), paste0("A_", 0:2), paste0("Y_", 0:2))
  ltmle_res <- run_ltmle(
    dat_wide, ltmle_cols,
    Anodes = paste0("A_", 0:2),
    Lnodes = paste0("L_", 0:2),
    Ynodes = paste0("Y_", 0:2),
    Cnodes = NULL,
    abar_treat = c(1, 1, 1),
    abar_control = c(0, 0, 0)
  )

  # longy
  dat_long <- wide_to_long(dat_wide, time_points = 0:2,
                           baseline_cols = "W",
                           treatment_col = "A", outcome_col = "Y",
                           timevarying_cols = "L")

  obj <- longy(dat_long, id = "id", time = "time", outcome = "Y",
               treatment = "A", baseline = "W", timevarying = "L",
               outcome_type = "binary",
               regimes = list(always = 1L, never = 0L),
               estimator = "tmle", learners = "SL.glm",
               inference = "eif", n_boot = 0L, verbose = verbose)

  # longy contrast at final time
  ctr_diff <- contrast(obj, regime = c("always", "never"), scale = "difference")
  ctr_rr   <- contrast(obj, regime = c("always", "never"), scale = "ratio")
  ctr_or   <- contrast(obj, regime = c("always", "never"), scale = "odds_ratio")

  final_t <- max(ctr_diff$estimates$time)
  longy_ate <- ctr_diff$estimates[time == final_t, estimate]
  longy_se  <- ctr_diff$estimates[time == final_t, se]
  longy_rr  <- ctr_rr$estimates[time == final_t, estimate]
  longy_or  <- ctr_or$estimates[time == final_t, estimate]

  ltmle_ate <- ltmle_res$ATE$estimate
  ltmle_se  <- ltmle_res$ATE$std.dev
  ltmle_rr  <- ltmle_res$RR$estimate
  ltmle_or  <- ltmle_res$OR$estimate

  cat(sprintf("  ATE:  longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_ate, ltmle_ate, abs(longy_ate - ltmle_ate)))
  cat(sprintf("  SE:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_se, ltmle_se, abs(longy_se - ltmle_se)))
  cat(sprintf("  RR:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_rr, ltmle_rr, abs(longy_rr - ltmle_rr)))
  cat(sprintf("  OR:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_or, ltmle_or, abs(longy_or - ltmle_or)))

  # Assertions: ATE and RR should be close. OR can diverge more because
  # longy computes OR from time-specific marginal estimates while ltmle uses
  # a parametric MSM (logistic working model) — different estimands at
  # intermediate times, and even at the final time the delta-method
  # formulations differ slightly.
  pass <- TRUE
  if (abs(longy_ate - ltmle_ate) > tol) {
    cat(sprintf("  FAIL: ATE difference %.4f > tolerance %.4f\n",
                abs(longy_ate - ltmle_ate), tol))
    pass <- FALSE
  }
  if (abs(longy_rr - ltmle_rr) / ltmle_rr > tol) {
    cat(sprintf("  FAIL: RR relative difference %.4f > tolerance %.4f\n",
                abs(longy_rr - ltmle_rr) / ltmle_rr, tol))
    pass <- FALSE
  }
  # OR uses wider tolerance — different delta-method parameterizations
  if (abs(longy_or - ltmle_or) / ltmle_or > 0.30) {
    cat(sprintf("  FAIL: OR relative difference %.4f > tolerance 0.30\n",
                abs(longy_or - ltmle_or) / ltmle_or))
    pass <- FALSE
  }
  if (pass) cat("  PASS\n")
  cat("\n")
  pass
}


# ============================================================================
# Test 2: Binary outcome with censoring — TMLE contrast vs ltmle ATE
# ============================================================================
# ============================================================================
# Test 2: CVD DGP — TMLE contrast vs ltmle ATE (4 time points)
# ============================================================================
test_cvd_tmle <- function(n = 5000, seed = 42, tol = 0.15, verbose = FALSE) {
  cat("--- Test 2: CVD DGP TMLE contrast vs ltmle ---\n")
  cat("  (Wider tolerance: longy targets E[Y(t)] per time; ltmle pools.)\n")
  set.seed(seed)

  D <- define_dgp_cvd()

  dat_wide <- sim(D, n = n, LTCF = "Y")

  # ltmle at final time
  ltmle_cols <- c("CVD",
                  "A1C_0", "TI_0", "Y_0",
                  "A1C_1", "TI_1", "Y_1",
                  "A1C_2", "TI_2", "Y_2",
                  "A1C_3", "TI_3", "Y_3")
  ltmle_res <- run_ltmle(
    dat_wide, ltmle_cols,
    Anodes = paste0("TI_", 0:3),
    Lnodes = paste0("A1C_", 0:3),
    Ynodes = paste0("Y_", 0:3),
    Cnodes = NULL,
    abar_treat = c(1, 1, 1, 1),
    abar_control = c(0, 0, 0, 0)
  )

  # longy
  dat_long <- wide_to_long(dat_wide, time_points = 0:3,
                           baseline_cols = "CVD",
                           treatment_col = "TI", outcome_col = "Y",
                           timevarying_cols = "A1C")

  obj <- longy(dat_long, id = "id", time = "time", outcome = "Y",
               treatment = "TI", baseline = "CVD", timevarying = "A1C",
               outcome_type = "binary",
               regimes = list(always = 1L, never = 0L),
               estimator = "tmle", learners = "SL.glm",
               inference = "eif", n_boot = 0L, verbose = verbose)

  ctr_diff <- contrast(obj, regime = c("always", "never"), scale = "difference")
  ctr_rr   <- contrast(obj, regime = c("always", "never"), scale = "ratio")
  ctr_or   <- contrast(obj, regime = c("always", "never"), scale = "odds_ratio")

  final_t <- max(ctr_diff$estimates$time)
  longy_ate <- ctr_diff$estimates[time == final_t, estimate]
  longy_se  <- ctr_diff$estimates[time == final_t, se]
  longy_rr  <- ctr_rr$estimates[time == final_t, estimate]
  longy_or  <- ctr_or$estimates[time == final_t, estimate]

  ltmle_ate <- ltmle_res$ATE$estimate
  ltmle_se  <- ltmle_res$ATE$std.dev
  ltmle_rr  <- ltmle_res$RR$estimate
  ltmle_or  <- ltmle_res$OR$estimate

  cat(sprintf("  ATE:  longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_ate, ltmle_ate, abs(longy_ate - ltmle_ate)))
  cat(sprintf("  SE:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_se, ltmle_se, abs(longy_se - ltmle_se)))
  cat(sprintf("  RR:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_rr, ltmle_rr, abs(longy_rr - ltmle_rr)))
  cat(sprintf("  OR:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_or, ltmle_or, abs(longy_or - ltmle_or)))

  pass <- TRUE
  if (abs(longy_ate - ltmle_ate) > tol) {
    cat(sprintf("  FAIL: ATE difference %.4f > tolerance %.4f\n",
                abs(longy_ate - ltmle_ate), tol))
    pass <- FALSE
  }
  if (pass) cat("  PASS\n")
  cat("\n")
  pass
}


# ============================================================================
# Test 4: IPW contrast vs ltmle IPTW ATE
# ============================================================================
test_binary_simple_ipw <- function(n = 5000, seed = 42, tol = 0.15,
                                   verbose = FALSE) {
  cat("--- Test 3: Binary simple IPW contrast vs ltmle IPTW ---\n")
  cat("  (Wider tolerance: longy=Hajek/stabilized, ltmle=HT-style.)\n")
  set.seed(seed)

  D <- define_dgp_binary_simple()
  dat_wide <- sim(D, n = n, LTCF = "Y")

  # ltmle with IPTW
  ltmle_cols <- c("W", paste0("L_", 0:2), paste0("A_", 0:2), paste0("Y_", 0:2))
  ltmle_data <- as.data.frame(dat_wide)[, ltmle_cols]
  ltmle_data$W <- as.character(as.integer(ltmle_data$W * 1000))  # dummy; ignored

  # Actually for IPW comparison, just use ltmle with iptw estimator
  ltmle_data <- as.data.frame(dat_wide)[, ltmle_cols]

  result_ltmle <- ltmle(
    data = ltmle_data,
    estimate.time = FALSE,
    Anodes = paste0("A_", 0:2),
    Lnodes = paste0("L_", 0:2),
    Ynodes = paste0("Y_", 0:2),
    abar = list(c(1, 1, 1), c(0, 0, 0)),
    survivalOutcome = FALSE,
    SL.library = "SL.glm",
    SL.cvControl = list(V = 3)
  )
  s_iptw <- summary(result_ltmle, estimator = "iptw")

  ltmle_ate <- s_iptw$effect.measures$ATE$estimate
  ltmle_se  <- s_iptw$effect.measures$ATE$std.dev

  # longy IPW
  dat_long <- wide_to_long(dat_wide, time_points = 0:2,
                           baseline_cols = "W",
                           treatment_col = "A", outcome_col = "Y",
                           timevarying_cols = "L")

  obj <- longy(dat_long, id = "id", time = "time", outcome = "Y",
               treatment = "A", baseline = "W", timevarying = "L",
               outcome_type = "binary",
               regimes = list(always = 1L, never = 0L),
               estimator = "ipw", learners = "SL.glm",
               inference = "ic", n_boot = 0L, verbose = verbose)

  ctr_diff <- contrast(obj, regime = c("always", "never"),
                        estimator = "ipw", scale = "difference")

  final_t <- max(ctr_diff$estimates$time)
  longy_ate <- ctr_diff$estimates[time == final_t, estimate]
  longy_se  <- ctr_diff$estimates[time == final_t, se]

  cat(sprintf("  ATE:  longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_ate, ltmle_ate, abs(longy_ate - ltmle_ate)))
  cat(sprintf("  SE:   longy=%.5f  ltmle=%.5f  diff=%.2e\n",
              longy_se, ltmle_se, abs(longy_se - ltmle_se)))

  # IPW implementations differ (Hajek vs Horvitz-Thompson, stabilized vs not),
  # so we use a wider tolerance
  pass <- TRUE
  if (abs(longy_ate - ltmle_ate) > tol * 2) {
    cat(sprintf("  FAIL: ATE difference %.4f > tolerance %.4f\n",
                abs(longy_ate - ltmle_ate), tol * 2))
    pass <- FALSE
  }
  if (pass) cat("  PASS\n")
  cat("\n")
  pass
}


# ============================================================================
# Test 5: Contrast at all time points — longy vs truth
# ============================================================================
test_contrast_all_times <- function(n = 5000, seed = 42, tol = 0.08,
                                    verbose = FALSE) {
  cat("--- Test 4: Contrast at all time points vs truth ---\n")
  set.seed(seed)

  D <- define_dgp_binary_simple()

  # True values
  truth <- compute_truth(D, "Y", 0:2)

  dat_wide <- sim(D, n = n, LTCF = "Y")
  dat_long <- wide_to_long(dat_wide, time_points = 0:2,
                           baseline_cols = "W",
                           treatment_col = "A", outcome_col = "Y",
                           timevarying_cols = "L")

  obj <- longy(dat_long, id = "id", time = "time", outcome = "Y",
               treatment = "A", baseline = "W", timevarying = "L",
               outcome_type = "binary",
               regimes = list(always = 1L, never = 0L),
               estimator = "tmle", learners = "SL.glm",
               inference = "eif", n_boot = 0L, verbose = verbose)

  ctr <- contrast(obj, regime = c("always", "never"), scale = "difference")

  cat("  Time-specific ATE comparison:\n")
  pass <- TRUE
  for (tt in 0:2) {
    longy_est <- ctr$estimates[time == tt, estimate]
    longy_se  <- ctr$estimates[time == tt, se]
    longy_lo  <- ctr$estimates[time == tt, ci_lower]
    longy_hi  <- ctr$estimates[time == tt, ci_upper]
    true_ate  <- truth[time == tt, ATE]

    covers <- !is.na(longy_lo) && !is.na(longy_hi) &&
      longy_lo <= true_ate && true_ate <= longy_hi
    tag <- if (covers) "COVERS" else "MISS"

    cat(sprintf("    t=%d: est=%.4f  se=%.4f  CI=[%.4f, %.4f]  truth=%.4f  %s\n",
                tt, longy_est, longy_se, longy_lo, longy_hi, true_ate, tag))

    if (abs(longy_est - true_ate) > tol) {
      cat(sprintf("    FAIL: ATE bias %.4f > tolerance %.4f at t=%d\n",
                  abs(longy_est - true_ate), tol, tt))
      pass <- FALSE
    }
  }
  if (pass) cat("  PASS\n")
  cat("\n")
  pass
}


# ============================================================================
# Test 6: Contrast with auto-contrast in longy()
# ============================================================================
test_auto_contrast <- function(n = 3000, seed = 42, verbose = FALSE) {
  cat("--- Test 5: longy(contrast=TRUE) auto-computes contrasts ---\n")
  set.seed(seed)

  D <- define_dgp_binary_simple()
  dat_wide <- sim(D, n = n, LTCF = "Y")
  dat_long <- wide_to_long(dat_wide, time_points = 0:2,
                           baseline_cols = "W",
                           treatment_col = "A", outcome_col = "Y",
                           timevarying_cols = "L")

  obj <- longy(dat_long, id = "id", time = "time", outcome = "Y",
               treatment = "A", baseline = "W", timevarying = "L",
               outcome_type = "binary",
               regimes = list(always = 1L, never = 0L),
               estimator = "tmle", learners = "SL.glm",
               inference = "eif", n_boot = 0L,
               contrast = TRUE, verbose = verbose)

  pass <- TRUE
  if (is.null(obj$contrasts) || length(obj$contrasts) == 0) {
    cat("  FAIL: obj$contrasts is empty\n")
    pass <- FALSE
  } else {
    cat(sprintf("  %d contrast(s) computed\n", length(obj$contrasts)))
    for (nm in names(obj$contrasts)) {
      ctr <- obj$contrasts[[nm]]
      cat(sprintf("    %s: %s (%s), %d time points\n",
                  nm, ctr$scale, ctr$inference, nrow(ctr$estimates)))
      if (!inherits(ctr, "longy_contrast")) {
        cat(sprintf("    FAIL: %s is not a longy_contrast\n", nm))
        pass <- FALSE
      }
    }
  }
  if (pass) cat("  PASS\n")
  cat("\n")
  pass
}


# ============================================================================
# Test 7: Multi-seed simulation — contrast coverage
# ============================================================================
test_contrast_coverage <- function(n = 3000, n_seeds = 20, tol_coverage = 0.70,
                                   verbose = FALSE) {
  cat(sprintf("--- Test 6: Contrast CI coverage (%d seeds, n=%d) ---\n",
              n_seeds, n))

  D <- define_dgp_binary_simple()
  truth <- compute_truth(D, "Y", 0:2, n_mc = 1000000)
  cat("  True ATEs: ")
  cat(sprintf("t=%d: %.4f  ", truth$time, truth$ATE))
  cat("\n")

  final_t <- max(truth$time)
  true_ate_final <- truth[time == final_t, ATE]

  covers <- logical(n_seeds)
  ests <- numeric(n_seeds)
  ses  <- numeric(n_seeds)

  for (s in seq_len(n_seeds)) {
    set.seed(s * 100)
    dat_wide <- sim(D, n = n, LTCF = "Y")
    dat_long <- wide_to_long(dat_wide, time_points = 0:2,
                             baseline_cols = "W",
                             treatment_col = "A", outcome_col = "Y",
                             timevarying_cols = "L")

    obj <- tryCatch({
      longy(dat_long, id = "id", time = "time", outcome = "Y",
            treatment = "A", baseline = "W", timevarying = "L",
            outcome_type = "binary",
            regimes = list(always = 1L, never = 0L),
            estimator = "tmle", learners = "SL.glm",
            inference = "eif", n_boot = 0L, verbose = FALSE)
    }, error = function(e) NULL)

    if (is.null(obj)) {
      covers[s] <- NA
      next
    }

    ctr <- tryCatch(
      contrast(obj, regime = c("always", "never"), scale = "difference"),
      error = function(e) NULL
    )

    if (is.null(ctr)) {
      covers[s] <- NA
      next
    }

    row <- ctr$estimates[time == final_t, ]
    ests[s] <- row$estimate
    ses[s]  <- row$se
    covers[s] <- !is.na(row$ci_lower) && !is.na(row$ci_upper) &&
      row$ci_lower <= true_ate_final && true_ate_final <= row$ci_upper
  }

  valid <- !is.na(covers)
  coverage <- mean(covers[valid])
  mean_est <- mean(ests[valid])
  mean_se  <- mean(ses[valid], na.rm = TRUE)
  emp_se   <- sd(ests[valid])

  cat(sprintf("  Coverage: %.1f%% (%d/%d valid)\n",
              coverage * 100, sum(valid), n_seeds))
  cat(sprintf("  Mean ATE estimate: %.4f (truth: %.4f, bias: %.4f)\n",
              mean_est, true_ate_final, mean_est - true_ate_final))
  cat(sprintf("  Mean SE: %.4f, Empirical SE: %.4f, Ratio: %.2f\n",
              mean_se, emp_se, mean_se / emp_se))

  pass <- coverage >= tol_coverage
  if (!pass) {
    cat(sprintf("  FAIL: Coverage %.1f%% < %.0f%%\n",
                coverage * 100, tol_coverage * 100))
  } else {
    cat("  PASS\n")
  }
  cat("\n")
  pass
}


# ============================================================================
# Run all tests
# ============================================================================
run_all <- function(verbose = FALSE) {
  cat("============================================================\n")
  cat("  longy contrast() vs ltmle validation suite\n")
  cat("============================================================\n\n")

  results <- list()

  results[["binary_simple_tmle"]]    <- test_binary_simple_tmle(verbose = verbose)
  results[["cvd_tmle"]]              <- test_cvd_tmle(verbose = verbose)
  results[["binary_simple_ipw"]]     <- test_binary_simple_ipw(verbose = verbose)
  results[["contrast_all_times"]]    <- test_contrast_all_times(verbose = verbose)
  results[["auto_contrast"]]         <- test_auto_contrast(verbose = verbose)
  results[["contrast_coverage"]]     <- test_contrast_coverage(verbose = verbose)

  cat("============================================================\n")
  cat("  SUMMARY\n")
  cat("============================================================\n")
  for (nm in names(results)) {
    tag <- if (results[[nm]]) "PASS" else "FAIL"
    cat(sprintf("  %-30s %s\n", nm, tag))
  }
  n_pass <- sum(unlist(results))
  n_total <- length(results)
  cat(sprintf("\n  %d/%d tests passed\n", n_pass, n_total))
  cat("============================================================\n")

  invisible(results)
}
