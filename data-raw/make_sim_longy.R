# Generate sim_longy dataset with known true causal effects
# DGP: 2000 subjects, 10 time points
# Baseline: W1 ~ N(0,1), W2 ~ Bernoulli(0.3)
# Time-varying: L1(t), L2(t) confounded by past treatment
# Treatment: A(t) confounded
# Censoring: C(t) informative dropout
# Observation: R(t) intermittent ~80% observed
# Outcome: Y(t) binary

set.seed(2026)

# --- Parameters ---
n <- 2000
K <- 10  # time points 0:9

# --- Helper ---
expit <- function(x) 1 / (1 + exp(-x))

# --- Generate observed data ---
generate_one_subject <- function(id) {
  W1 <- rnorm(1)
  W2 <- rbinom(1, 1, 0.3)

  rows <- list()
  alive <- TRUE
  A_prev <- 0  # treatment at t-1

  for (tt in 0:(K - 1)) {
    if (!alive) break

    # Time-varying confounders (affected by past treatment)
    L1 <- rnorm(1, mean = 0.3 * W1 + 0.2 * A_prev + 0.05 * tt, sd = 1)
    L2 <- rbinom(1, 1, expit(-0.5 + 0.3 * L1 + 0.1 * A_prev))

    # Treatment (confounded by L1, L2, W1, W2)
    p_a <- expit(-0.2 + 0.4 * L1 + 0.3 * W1 - 0.2 * L2 + 0.3 * W2)
    A <- rbinom(1, 1, p_a)

    # Censoring: informative dropout
    p_c <- expit(-3.5 + 0.2 * L1 - 0.1 * A + 0.15 * tt + 0.1 * W1)
    C <- rbinom(1, 1, p_c)

    if (C == 1) {
      rows[[length(rows) + 1]] <- data.frame(
        id = id, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = "censored", R = 0L, Y = NA_real_,
        stringsAsFactors = FALSE
      )
      alive <- FALSE
      next
    }

    # Observation: intermittent (~80% observed)
    p_r <- expit(1.5 - 0.2 * L1 + 0.1 * A - 0.05 * tt)
    R <- rbinom(1, 1, p_r)

    # Outcome (only observed when R=1)
    if (R == 1) {
      p_y <- expit(-2 + 0.3 * L1 + 0.4 * A + 0.2 * W1 + 0.15 * tt - 0.1 * W2)
      Y <- rbinom(1, 1, p_y)
    } else {
      Y <- NA_real_
    }

    rows[[length(rows) + 1]] <- data.frame(
      id = id, time = tt, W1 = W1, W2 = W2,
      L1 = L1, L2 = L2, A = A, C = "uncensored", R = R, Y = Y,
      stringsAsFactors = FALSE
    )

    A_prev <- A
  }

  do.call(rbind, rows)
}

sim_longy <- do.call(rbind, lapply(1:n, generate_one_subject))
sim_longy$id <- as.integer(sim_longy$id)
sim_longy$time <- as.integer(sim_longy$time)
sim_longy$A <- as.integer(sim_longy$A)
sim_longy$R <- as.integer(sim_longy$R)
sim_longy$W2 <- as.integer(sim_longy$W2)
sim_longy$L2 <- as.integer(sim_longy$L2)

# --- Compute true causal effects via large-sample Monte Carlo ---
# Simulate potential outcomes under always-treat and never-treat
set.seed(99999)
n_mc <- 100000

compute_true_effects <- function(n_mc, K) {
  W1 <- rnorm(n_mc)
  W2 <- rbinom(n_mc, 1, 0.3)

  EY1 <- numeric(K)  # E[Y_t(a=1)]
  EY0 <- numeric(K)  # E[Y_t(a=0)]

  # Always treat: A(t) = 1 for all t
  A_prev_1 <- rep(0, n_mc)
  # Never treat: A(t) = 0 for all t
  A_prev_0 <- rep(0, n_mc)

  for (tt in 0:(K - 1)) {
    # Under always-treat
    L1_1 <- rnorm(n_mc, mean = 0.3 * W1 + 0.2 * A_prev_1 + 0.05 * tt, sd = 1)
    p_y_1 <- expit(-2 + 0.3 * L1_1 + 0.4 * 1 + 0.2 * W1 + 0.15 * tt - 0.1 * W2)
    EY1[tt + 1] <- mean(p_y_1)
    A_prev_1 <- rep(1, n_mc)

    # Under never-treat
    L1_0 <- rnorm(n_mc, mean = 0.3 * W1 + 0.2 * A_prev_0 + 0.05 * tt, sd = 1)
    p_y_0 <- expit(-2 + 0.3 * L1_0 + 0.4 * 0 + 0.2 * W1 + 0.15 * tt - 0.1 * W2)
    EY0[tt + 1] <- mean(p_y_0)
    A_prev_0 <- rep(0, n_mc)
  }

  data.frame(time = 0:(K - 1), EY1 = EY1, EY0 = EY0)
}

true_effects <- compute_true_effects(n_mc, K)

attr(sim_longy, "true_effects") <- true_effects

# Save
save(sim_longy, file = file.path("..", "data", "sim_longy.rda"), compress = "xz")

cat("sim_longy created:\n")
cat(sprintf("  %d rows, %d subjects, %d time points\n",
            nrow(sim_longy), length(unique(sim_longy$id)), K))
cat("\nTrue effects:\n")
print(true_effects)
