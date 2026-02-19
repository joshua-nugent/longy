# ==========================================================================
# Comparison: longy vs lmtp TMLE with competing risks
#
# - n = 10,000, K = 5 time points
# - Time-varying treatment, time-varying covariates
# - Censoring + competing risks + primary event
# - Cross-fitting (5 folds) in both
# - future parallelization in both
# - GLM learners only (fair comparison)
# ==========================================================================

library(longy)
library(lmtp)
library(data.table)
library(future)

# Use multicore parallelization
plan(multisession, workers = 8)
cat("future plan:", class(plan())[1], "\n\n")

# --- DGP: generate in long format, convert to wide for lmtp ---
generate_cr_data <- function(n = 10000, K = 5, seed = 2026) {
  set.seed(seed)
  rows <- vector("list", n * K)
  idx <- 0L
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.4)
    alive <- TRUE
    primary <- FALSE
    competing <- FALSE
    for (tt in 0:(K - 1)) {
      if (!alive) {
        # After censoring: still emit rows (needed for wide reshape)
        idx <- idx + 1L
        rows[[idx]] <- list(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = NA_real_, A = NA_integer_, C = 1L,
          Y = NA_integer_, D = NA_integer_
        )
        next
      }

      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      p_a <- plogis(-0.2 + 0.3 * L1 + 0.2 * W1 - 0.1 * W2)
      A <- rbinom(1, 1, p_a)

      # Censoring (~5% per period)
      p_c <- plogis(-3.5 + 0.15 * L1 - 0.1 * A + 0.08 * tt)
      C <- rbinom(1, 1, p_c)
      if (C == 1) {
        idx <- idx + 1L
        rows[[idx]] <- list(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, A = A, C = 1L,
          Y = NA_integer_, D = NA_integer_
        )
        alive <- FALSE
        next
      }

      if (primary) {
        Y <- 1L; D_val <- 0L
      } else if (competing) {
        Y <- 0L; D_val <- 1L
      } else {
        p_primary <- plogis(-3.5 + 0.3 * L1 + 0.25 * A + 0.12 * tt + 0.15 * W1)
        p_competing <- plogis(-4 + 0.2 * L1 - 0.15 * A + 0.1 * tt + 0.1 * W2)
        u <- runif(1)
        if (u < p_primary) {
          Y <- 1L; D_val <- 0L; primary <- TRUE
        } else if (u < p_primary + p_competing) {
          Y <- 0L; D_val <- 1L; competing <- TRUE
        } else {
          Y <- 0L; D_val <- 0L
        }
      }

      idx <- idx + 1L
      rows[[idx]] <- list(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, A = A, C = 0L,
        Y = Y, D = D_val
      )
    }
  }

  rbindlist(rows[seq_len(idx)])
}

# --- Generate data ---
cat("Generating data (n=10,000, K=5)...\n")
dat_long <- generate_cr_data(n = 10000, K = 5, seed = 2026)

cat(sprintf("  Long format: %d rows, %d subjects\n", nrow(dat_long), uniqueN(dat_long$id)))
cat(sprintf("  Primary events: %d\n", sum(dat_long$Y == 1, na.rm = TRUE)))
cat(sprintf("  Competing events: %d\n", sum(dat_long$D == 1, na.rm = TRUE)))
cat(sprintf("  Censored rows: %d\n", sum(dat_long$C == 1, na.rm = TRUE)))

# --- Convert to wide format for lmtp ---
# lmtp conventions:
#   - Wide: one row per subject
#   - C = 1 means "not censored" (opposite of longy)
#   - Treatment, covariates, outcome, competing: separate columns per time
#   - Outcome/competing LOCF: once Y=1, all subsequent Y=1; once D=1, all subsequent D=1
#   - After censoring (C=0 in lmtp): subsequent Y/D should remain as their LOCF'd values

cat("\nConverting to wide format for lmtp...\n")

K <- 5
time_vals <- 0:(K - 1)

# Pivot each variable to wide
make_wide <- function(dt, var, prefix) {
  w <- dcast(dt, id ~ time, value.var = var)
  old_names <- as.character(time_vals)
  new_names <- paste0(prefix, seq_len(K))
  setnames(w, old_names, new_names)
  w
}

wide_A <- make_wide(dat_long, "A", "A")
wide_L1 <- make_wide(dat_long, "L1", "L1_")
wide_C <- make_wide(dat_long, "C", "C")
wide_Y <- make_wide(dat_long, "Y", "Y")
wide_D <- make_wide(dat_long, "D", "D")

# Baseline covariates (constant within subject)
wide_W <- unique(dat_long[, .(id, W1, W2)])

# Merge all
dat_wide <- Reduce(function(a, b) merge(a, b, by = "id"),
                   list(wide_W, wide_A, wide_L1, wide_C, wide_Y, wide_D))

# Flip censoring: longy C=1 (censored) -> lmtp C=1 (not censored)
for (k in seq_len(K)) {
  cname <- paste0("C", k)
  dat_wide[, (cname) := 1L - get(cname)]
}

# LOCF for outcome after censoring (lmtp expects this)
# After C=0 (censored in lmtp), carry forward last Y and D
for (i in 2:K) {
  y_cur <- paste0("Y", i)
  y_prev <- paste0("Y", i - 1)
  d_cur <- paste0("D", i)
  d_prev <- paste0("D", i - 1)
  c_prev <- paste0("C", i - 1)

  # If censored at prev time (C_prev = 0), carry forward
  dat_wide[get(c_prev) == 0 & is.na(get(y_cur)), (y_cur) := get(y_prev)]
  dat_wide[get(c_prev) == 0 & is.na(get(d_cur)), (d_cur) := get(d_prev)]
}

# LOCF for outcome after event (lmtp expects absorbing states carried forward)
for (i in 2:K) {
  y_cur <- paste0("Y", i)
  y_prev <- paste0("Y", i - 1)
  d_cur <- paste0("D", i)
  d_prev <- paste0("D", i - 1)

  dat_wide[get(y_prev) == 1 & (is.na(get(y_cur)) | get(y_cur) == 0), (y_cur) := 1L]
  dat_wide[get(d_prev) == 1 & (is.na(get(d_cur)) | get(d_cur) == 0), (d_cur) := 1L]
}

# Fill remaining NAs in Y/D with 0 (subjects who were censored before any event)
for (k in seq_len(K)) {
  yname <- paste0("Y", k)
  dname <- paste0("D", k)
  dat_wide[is.na(get(yname)), (yname) := 0L]
  dat_wide[is.na(get(dname)), (dname) := 0L]
}

# Fill NA in A and L1 after censoring (lmtp needs complete data for treatment/covariates)
# Use 0 for treatment, carry forward for L1
for (k in seq_len(K)) {
  aname <- paste0("A", k)
  lname <- paste0("L1_", k)
  dat_wide[is.na(get(aname)), (aname) := 0L]
  if (k > 1) {
    lprev <- paste0("L1_", k - 1)
    dat_wide[is.na(get(lname)), (lname) := get(lprev)]
  } else {
    dat_wide[is.na(get(lname)), (lname) := 0]
  }
}

# Convert to data.frame (lmtp doesn't accept data.table)
dat_wide_df <- as.data.frame(dat_wide)
dat_wide_df$id <- NULL  # lmtp doesn't need an id column in the data (unless clusters)

cat(sprintf("  Wide format: %d rows, %d cols\n", nrow(dat_wide_df), ncol(dat_wide_df)))

# --- For longy: clean long format (drop post-censoring rows) ---
# longy expects: after censoring (C=1), no more rows OR rows with C=1 and Y=NA
# Keep all rows but longy handles NAs after censoring naturally
dat_long_df <- as.data.frame(dat_long)

# =========================================================================
# Run longy
# =========================================================================
cat("\n=== Running longy TMLE (cross-fit=5) ===\n")
t_longy_start <- Sys.time()

res_longy <- longy(
  data = dat_long_df,
  id = "id", time = "time", outcome = "Y",
  treatment = "A", censoring = "C",
  baseline = c("W1", "W2"), timevarying = "L1",
  outcome_type = "survival", competing = "D",
  regimes = list(always = 1L, never = 0L),
  estimator = "tmle",
  learners = NULL,  # glm only
  cross_fit = 5L, cross_fit_seed = 42L,
  inference = "eif",
  sl_fn = "SuperLearner",  # doesn't matter with NULL learners
  verbose = TRUE
)

t_longy_end <- Sys.time()
longy_time <- difftime(t_longy_end, t_longy_start, units = "secs")

cat(sprintf("\nlongy completed in %.1f seconds\n", longy_time))

# =========================================================================
# Run lmtp
# =========================================================================
cat("\n=== Running lmtp TMLE (folds=5) ===\n")

trt_cols <- paste0("A", 1:K)
outcome_cols <- paste0("Y", 1:K)
cens_cols <- paste0("C", 1:K)
compete_cols <- paste0("D", 1:K)
baseline_cols <- c("W1", "W2")
timevarying_cols <- lapply(1:K, function(k) paste0("L1_", k))

# lmtp shift function: always treat (A = 1)
shift_always <- function(data, trt) rep(1L, nrow(data))
shift_never  <- function(data, trt) rep(0L, nrow(data))

t_lmtp_start <- Sys.time()

# Use progressr for progress bars
progressr::handlers(global = TRUE)

# Always treat
cat("  lmtp: always treat...\n")
res_lmtp_always <- lmtp_survival(
  data = dat_wide_df,
  trt = trt_cols,
  outcomes = outcome_cols,
  baseline = baseline_cols,
  time_vary = timevarying_cols,
  cens = cens_cols,
  compete = compete_cols,
  shift = shift_always,
  estimator = "lmtp_tmle",
  mtp = FALSE,
  folds = 5,
  learners_outcome = "SL.glm",
  learners_trt = "SL.glm"
)

# Never treat
cat("  lmtp: never treat...\n")
res_lmtp_never <- lmtp_survival(
  data = dat_wide_df,
  trt = trt_cols,
  outcomes = outcome_cols,
  baseline = baseline_cols,
  time_vary = timevarying_cols,
  cens = cens_cols,
  compete = compete_cols,
  shift = shift_never,
  estimator = "lmtp_tmle",
  mtp = FALSE,
  folds = 5,
  learners_outcome = "SL.glm",
  learners_trt = "SL.glm"
)

t_lmtp_end <- Sys.time()
lmtp_time <- difftime(t_lmtp_end, t_lmtp_start, units = "secs")

cat(sprintf("\nlmtp completed in %.1f seconds\n", lmtp_time))

# =========================================================================
# Compare results
# =========================================================================
cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON: longy vs lmtp TMLE with competing risks\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

# Extract longy estimates
longy_always <- res_longy$always$estimates
longy_never  <- res_longy$never$estimates

# Extract lmtp estimates
# NOTE: lmtp_survival() returns event-free probabilities: P(Y=0 by t) = 1 - CIF
# longy returns CIF directly: P(Y=1 by t)
# So we convert lmtp to CIF scale: 1 - lmtp_estimate
lmtp_always_est <- sapply(res_lmtp_always, function(x) {
  tidy_x <- lmtp::tidy(x)
  1 - tidy_x$estimate  # convert to CIF
})
lmtp_always_se <- sapply(res_lmtp_always, function(x) {
  tidy_x <- lmtp::tidy(x)
  tidy_x$std.error  # SE is the same on either scale
})

lmtp_never_est <- sapply(res_lmtp_never, function(x) {
  tidy_x <- lmtp::tidy(x)
  1 - tidy_x$estimate  # convert to CIF
})
lmtp_never_se <- sapply(res_lmtp_never, function(x) {
  tidy_x <- lmtp::tidy(x)
  tidy_x$std.error
})

cat("\n--- Always Treat (CIF) ---\n")
cat("Note: lmtp estimates converted from event-free P to CIF = 1 - P\n")
cat(sprintf("%-6s  %-12s  %-12s  %-10s  %-12s  %-12s\n",
            "Time", "longy", "lmtp", "Diff", "longy_SE", "lmtp_SE"))
cat(rep("-", 76), "\n", sep = "")
for (k in seq_len(K)) {
  le <- longy_always$estimate[k]
  lse <- if ("se" %in% names(longy_always)) longy_always$se[k] else NA
  me <- lmtp_always_est[k]
  mse <- lmtp_always_se[k]
  cat(sprintf("%-6d  %-12.6f  %-12.6f  %-10.6f  %-12.6f  %-12.6f\n",
              k - 1, le, me, le - me, lse, mse))
}

cat("\n--- Never Treat (CIF) ---\n")
cat(sprintf("%-6s  %-12s  %-12s  %-10s  %-12s  %-12s\n",
            "Time", "longy", "lmtp", "Diff", "longy_SE", "lmtp_SE"))
cat(rep("-", 76), "\n", sep = "")
for (k in seq_len(K)) {
  le <- longy_never$estimate[k]
  lse <- if ("se" %in% names(longy_never)) longy_never$se[k] else NA
  me <- lmtp_never_est[k]
  mse <- lmtp_never_se[k]
  cat(sprintf("%-6d  %-12.6f  %-12.6f  %-10.6f  %-12.6f  %-12.6f\n",
              k - 1, le, me, le - me, lse, mse))
}

cat(sprintf("\nTiming: longy=%.1fs, lmtp=%.1fs\n", longy_time, lmtp_time))

# Reset future plan
plan(sequential)
