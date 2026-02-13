# Internal utility functions for longy
# None of these are exported

#' Bound values to a range
#' @noRd
.bound <- function(x, lower = 0.005, upper = 0.995) {
  pmin(pmax(x, lower), upper)
}

#' Expit (inverse logit)
#' @noRd
.expit <- function(x) {
  1 / (1 + exp(-x))
}

#' Logit
#' @noRd
.logit <- function(p) {
  log(p / (1 - p))
}

#' Effective sample size (Kish's formula)
#' @noRd
.ess <- function(w) {
  sum(w)^2 / sum(w^2)
}

#' Adaptive CV fold selection based on effective sample size
#'
#' Ported from user's get_neff_and_V. Determines appropriate number of
#' CV folds based on the effective sample size of a binary outcome vector.
#' @param y Binary outcome vector
#' @return List with n_eff, V, n, p_hat, n_rare
#' @noRd
.adaptive_cv_folds <- function(y) {
  n <- length(y)
  p_hat <- mean(y)
  n_rare <- n * min(p_hat, 1 - p_hat)
  n_eff <- min(n, 5 * n_rare)

  V <- if (n_eff < 30) {
    as.integer(n_eff)
  } else if (n_eff < 500) {
    20L
  } else if (n_eff < 5000) {
    10L
  } else if (n_eff < 10000) {
    5L
  } else {
    2L
  }

  # Safety: V must be >= 2 and <= n
  V <- as.integer(max(2L, min(V, n)))

  list(n_eff = n_eff, V = V, n = n, p_hat = p_hat, n_rare = n_rare)
}

#' Safe SuperLearner wrapper with glm fallback
#'
#' Tries to fit SuperLearner. On failure, falls back to glm.
#' @param Y Binary outcome vector
#' @param X Covariate data.frame
#' @param family Model family (default binomial)
#' @param learners Character vector of SL learner names
#' @param cv_folds Number of CV folds
#' @param obs_weights Numeric vector of observation/sampling weights (same
#'   length as Y). Passed as \code{obsWeights} to SuperLearner and
#'   \code{weights} to glm. NULL means equal weights.
#' @param verbose Logical
#' @return List with predictions (numeric vector) and fit object
#' @noRd
.safe_sl <- function(Y, X, family = stats::binomial(),
                     learners = c("SL.glm", "SL.mean"),
                     cv_folds = 10L, obs_weights = NULL,
                     verbose = FALSE) {
  X <- as.data.frame(X)

  # If SuperLearner is available and learners specified, try it
  if (!is.null(learners) && length(learners) > 0 &&
      requireNamespace("SuperLearner", quietly = TRUE)) {
    # SuperLearner looks up learner and screening functions in env=;
    # point it at the SuperLearner namespace so SL.glm, All, etc. are found
    fit <- tryCatch(
      {
        sl_args <- list(
          Y = Y, X = X, family = family,
          SL.library = learners,
          cvControl = list(V = cv_folds),
          env = asNamespace("SuperLearner")
        )
        if (!is.null(obs_weights)) {
          sl_args$obsWeights <- obs_weights
        }
        sl_fit <- do.call(SuperLearner::SuperLearner, sl_args)
        list(
          predictions = as.numeric(sl_fit$SL.predict),
          fit = sl_fit,
          method = "SuperLearner",
          sl_risk = sl_fit$cvRisk,
          sl_coef = sl_fit$coef
        )
      },
      error = function(e) {
        warning(sprintf("SuperLearner failed: %s. Falling back to glm.",
                        e$message), call. = FALSE)
        NULL
      }
    )
    if (!is.null(fit)) return(fit)
  }

  # Fallback: glm
  df <- X
  df$.Y <- Y
  glm_fit <- tryCatch(
    {
      if (!is.null(obs_weights)) {
        stats::glm(.Y ~ ., data = df, family = family,
                    weights = obs_weights)
      } else {
        stats::glm(.Y ~ ., data = df, family = family)
      }
    },
    error = function(e) {
      stop(sprintf("Both SuperLearner and glm failed. glm error: %s",
                   e$message), call. = FALSE)
    }
  )
  preds <- as.numeric(stats::predict(glm_fit, newdata = X, type = "response"))
  list(predictions = preds, fit = glm_fit, method = "glm")
}

#' Verbose messaging helper
#' @noRd
.vmsg <- function(fmt, ..., verbose = TRUE) {
  if (verbose) {
    message(sprintf(fmt, ...))
  }
}

#' Convert factors/characters to dummy variables in a data.table
#' @param dt data.table
#' @param cols Character vector of column names to convert
#' @return Modified data.table with factor columns replaced by dummies
#' @noRd
.factors_to_dummies <- function(dt, cols) {
  if (length(cols) == 0) return(dt)

  for (col in cols) {
    vals <- dt[[col]]
    if (is.character(vals)) vals <- as.factor(vals)
    lvls <- levels(vals)
    if (length(lvls) <= 1) next

    # Create dummies for all levels except first (reference)
    for (i in seq_along(lvls)[-1]) {
      new_name <- paste0(col, "_", lvls[i])
      data.table::set(dt, j = new_name, value = as.integer(vals == lvls[i]))
    }
    # Remove original column
    data.table::set(dt, j = col, value = NULL)
  }
  dt
}
