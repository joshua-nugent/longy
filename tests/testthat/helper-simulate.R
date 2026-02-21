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
      C_event <- rbinom(1, 1, p_c)
      if (C_event == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = "censored", R = 0L, Y = NA_real_,
          stringsAsFactors = FALSE
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
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = R, Y = Y,
        stringsAsFactors = FALSE
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
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = R, Y = Y,
        stringsAsFactors = FALSE
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
        L1 = L1, L2 = L2, A = A, R = 1L, Y = Y
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with continuous outcome
simulate_continuous_outcome <- function(n = 200, K = 5, seed = 321) {
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
      C_event <- rbinom(1, 1, p_c)
      if (C_event == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = "censored", R = 0L, Y = NA_real_,
          stringsAsFactors = FALSE
        )
        alive <- FALSE
        next
      }
      p_r <- plogis(1.5 - 0.2 * L1 + 0.1 * A)
      R <- rbinom(1, 1, p_r)
      if (R == 1) {
        Y <- rnorm(1, mean = 0.5 * L1 + 0.3 * A + 0.2 * W1, sd = 1)
      } else {
        Y <- NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = R, Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with survival outcome (absorbing event)
simulate_survival_outcome <- function(n = 200, K = 5, seed = 654) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    alive <- TRUE
    event_occurred <- FALSE
    for (tt in 0:(K - 1)) {
      if (!alive) break
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.3 * L1))
      p_a <- plogis(-0.2 + 0.4 * L1 + 0.3 * W1 - 0.2 * L2)
      A <- rbinom(1, 1, p_a)
      p_c <- plogis(-3 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
      C_event <- rbinom(1, 1, p_c)
      if (C_event == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = "censored", R = 0L, Y = NA_real_,
          stringsAsFactors = FALSE
        )
        alive <- FALSE
        next
      }
      p_r <- plogis(1.5 - 0.2 * L1 + 0.1 * A)
      R <- rbinom(1, 1, p_r)
      if (R == 1) {
        if (event_occurred) {
          Y <- 1L  # absorbing: once event, stays event
        } else {
          p_event <- plogis(-3 + 0.3 * L1 + 0.2 * A + 0.15 * tt)
          Y <- rbinom(1, 1, p_event)
          if (Y == 1L) event_occurred <- TRUE
        }
      } else {
        Y <- NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = as.integer(R), Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with carryover effects (past A affects future L1)
#' No censoring, always observed, continuous outcome.
#' L1(t) = 0.3*W1 + 0.2*A(t-1) + 0.1*t + N(0,1)  [A(-1)=0]
#' A(t) ~ Bern(plogis(0.5*L1(t)))
#' Y(t) = 0.5*L1(t) + 0.3*A(t) + 0.2*W1 + N(0,1)
simulate_carryover_continuous <- function(n = 200, K = 5, seed = 42) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    A_prev <- 0L
    for (tt in 0:(K - 1)) {
      L1 <- 0.3 * W1 + 0.2 * A_prev + 0.1 * tt + rnorm(1)
      A <- rbinom(1, 1, plogis(0.5 * L1))
      Y <- 0.5 * L1 + 0.3 * A + 0.2 * W1 + rnorm(1)
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, L1 = L1,
        A = A, R = 1L, Y = Y
      )
      A_prev <- A
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
      C_event <- rbinom(1, 1, p_c)
      if (C_event == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = "censored", R = 0L, Y = NA_real_,
          stringsAsFactors = FALSE
        )
        alive <- FALSE
        next
      }
      p_y <- plogis(-1.5 + 0.3 * L1 + 0.3 * A + 0.2 * W1 + 0.05 * tt)
      Y <- rbinom(1, 1, p_y)
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = 1L, Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with competing risks (primary event + competing event)
#' Both Y and D are absorbing binary indicators. Mutually exclusive.
simulate_competing_risks <- function(n = 200, K = 5, seed = 111) {
  set.seed(seed)
  rows <- list()
  for (i in seq_len(n)) {
    W1 <- rnorm(1)
    W2 <- rbinom(1, 1, 0.3)
    alive <- TRUE
    primary_occurred <- FALSE
    competing_occurred <- FALSE
    for (tt in 0:(K - 1)) {
      if (!alive) break
      L1 <- rnorm(1, mean = 0.3 * W1 + 0.1 * tt)
      L2 <- rbinom(1, 1, plogis(-0.5 + 0.3 * L1))
      p_a <- plogis(-0.2 + 0.4 * L1 + 0.3 * W1 - 0.2 * L2)
      A <- rbinom(1, 1, p_a)
      # Censoring
      p_c <- plogis(-3 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
      C_event <- rbinom(1, 1, p_c)
      if (C_event == 1) {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = "censored", R = 0L,
          Y = NA_integer_, D = NA_integer_,
          stringsAsFactors = FALSE
        )
        alive <- FALSE
        next
      }
      R <- 1L  # always observed when uncensored
      if (primary_occurred) {
        Y <- 1L; D <- 0L  # absorbing primary
      } else if (competing_occurred) {
        Y <- 0L; D <- 1L  # absorbing competing
      } else {
        # Competing risk: primary event
        p_primary <- plogis(-3 + 0.3 * L1 + 0.2 * A + 0.15 * tt)
        # Competing event (e.g., death from other cause)
        p_competing <- plogis(-3.5 + 0.2 * L1 - 0.1 * A + 0.1 * tt)
        # Draw: at most one event per time
        u <- runif(1)
        if (u < p_primary) {
          Y <- 1L; D <- 0L
          primary_occurred <- TRUE
        } else if (u < p_primary + p_competing) {
          Y <- 0L; D <- 1L
          competing_occurred <- TRUE
        } else {
          Y <- 0L; D <- 0L
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        id = i, time = tt, W1 = W1, W2 = W2,
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = R,
        Y = Y, D = D,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Simulate data with multiple censoring causes
#' C_status has values "uncensored", "death", "ltfu"
simulate_multi_censoring <- function(n = 200, K = 5, seed = 222) {
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
      # Two censoring causes
      p_death <- plogis(-4 + 0.2 * L1 + 0.1 * tt)
      p_ltfu <- plogis(-3.5 + 0.15 * L1 - 0.1 * A + 0.1 * tt)
      u <- runif(1)
      if (u < p_death) {
        C_status <- "death"
      } else if (u < p_death + p_ltfu) {
        C_status <- "ltfu"
      } else {
        C_status <- "uncensored"
      }
      if (C_status != "uncensored") {
        rows[[length(rows) + 1]] <- data.frame(
          id = i, time = tt, W1 = W1, W2 = W2,
          L1 = L1, L2 = L2, A = A, C = C_status, R = 0L, Y = NA_real_,
          stringsAsFactors = FALSE
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
        L1 = L1, L2 = L2, A = A, C = "uncensored", R = R, Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
