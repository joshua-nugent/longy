# Test data generators for longy tests
# These create small datasets suitable for fast test execution

#' Simulate confounded + censored + intermittent observation data
simulate_test_data <- function(n = 100, K = 5, seed = 42) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    alive <- TRUE
    for (tt in 0:(K - 1)) {
      if (!alive) break
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.3 * L1))
      p_a <- plogis(-0.2 + 0.4 * L1 + 0.3 * W1 - 0.2 * L2)
      A <- rbinom(1, 1, p_a)
      p_c <- plogis(-3 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
      C <- rbinom(1, 1, p_c)
      if (C == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = 1L, R = 0L, Y = NA_real_
        )
        alive <- FALSE
        next
      }
      p_r <- plogis(1.5 - 0.2 * L1 + 0.1 * A)
      R <- rbinom(1, 1, p_r)
      if (R == 1) {
        p_y <- plogis(-2 + 0.3 * L1 + 0.2 * A + 0.2 * W1 + 0.1 * tt)
        Y <- rbinom(1, 1, p_y)
      } else {
        Y <- NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = 0L, R = R, Y = Y
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with NO confounding (A independent of L)
simulate_no_confounding <- function(n = 200, K = 5, seed = 123) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    for (tt in 0:(K - 1)) {
      L1 <- rnorm(1, mean = 0.1 * tt)
      L2 <- rbinom(1, 1, 0.4)
      # Treatment independent of covariates
      A <- rbinom(1, 1, 0.5)
      p_r <- plogis(2)  # ~88% observed
      R <- rbinom(1, 1, p_r)
      if (R == 1) {
        p_y <- plogis(-1 + 0.3 * A + 0.1 * W1 + 0.05 * tt)
        Y <- rbinom(1, 1, p_y)
      } else {
        Y <- NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = 0L, R = R, Y = Y
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with NO censoring (everyone stays)
simulate_no_censoring <- function(n = 200, K = 5, seed = 456) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    for (tt in 0:(K - 1)) {
      L1 <- rnorm(1, mean = 0.2 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.2 * L1))
      p_a <- plogis(-0.1 + 0.3 * L1 + 0.2 * W1)
      A <- rbinom(1, 1, p_a)
      R <- 1L  # always observed
      p_y <- plogis(-1.5 + 0.3 * L1 + 0.3 * A + 0.2 * W1 + 0.05 * tt)
      Y <- rbinom(1, 1, p_y)
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = 0L, R = 1L, Y = Y
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with censoring but outcome always observed (R=1 always)
simulate_always_observed <- function(n = 200, K = 5, seed = 789) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    alive <- TRUE
    for (tt in 0:(K - 1)) {
      if (!alive) break
      L1 <- rnorm(1, mean = 0.2 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.3 + 0.2 * L1))
      p_a <- plogis(-0.1 + 0.3 * L1 + 0.2 * W1)
      A <- rbinom(1, 1, p_a)
      p_c <- plogis(-3 + 0.15 * L1 + 0.1 * tt)
      C <- rbinom(1, 1, p_c)
      if (C == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = 1L, R = 0L, Y = NA_real_
        )
        alive <- FALSE
        next
      }
      p_y <- plogis(-1.5 + 0.3 * L1 + 0.3 * A + 0.2 * W1 + 0.05 * tt)
      Y <- rbinom(1, 1, p_y)
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = 0L, R = 1L, Y = Y
      )
    }
  }
  do.call(rbind, rows)
}
