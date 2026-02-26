# Internal utility functions for longy
# None of these are exported

#' Get lag covariate column names for a given time index
#'
#' Returns the \code{.longy_lag_*} column names to append to the covariate
#' set at each time point. At \code{time_index = 1} (first time), returns
#' nothing. At \code{time_index = 4} with \code{k = Inf}, returns all
#' \code{_lag1}, \code{_lag2}, \code{_lag3} columns.
#'
#' @param nodes List of node names (must contain \code{lag_vars} and
#'   \code{lag_k}).
#' @param time_index Integer. 1-based index of the current time point in the
#'   time loop.
#' @return Character vector of lag column names.
#' @noRd
.get_lag_covariates <- function(nodes, time_index) {
  k <- nodes$lag_k
  if (is.null(k) || k == 0) return(character(0))
  lag_vars <- nodes$lag_vars
  if (is.null(lag_vars) || length(lag_vars) == 0) return(character(0))
  max_lag <- min(k, time_index - 1)
  if (max_lag <= 0) return(character(0))
  unlist(lapply(lag_vars, function(col) {
    paste0(".longy_lag_", col, "_", seq_len(max_lag))
  }))
}

#' Extract or validate a longy_data object from any longy type
#'
#' Accepts \code{longy_data} (pass-through), \code{longy_result} (extracts
#' \code{$obj}), or legacy \code{longy_results} (extracts first element's
#' \code{$obj}). This replaces the older \code{.extract_longy_data()} helper
#' and all \code{stopifnot(inherits(obj, "longy_data"))} call sites.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object.
#' @return A \code{longy_data} object.
#' @noRd
.as_longy_data <- function(obj) {
  if (inherits(obj, "longy_data")) return(obj)
  if (inherits(obj, "longy_result") && !is.null(obj$obj)) return(obj$obj)
  if (inherits(obj, "longy_results") && length(obj) > 0 &&
      !is.null(obj[[1]]$obj)) return(obj[[1]]$obj)
  stop("Expected a longy_data, longy_result, or longy_results object.",
       call. = FALSE)
}

#' Resolve regime argument for fit_*/compute_weights/estimate_* functions
#'
#' When \code{regime} is NULL, defaults to all defined regimes. Validates
#' that all requested regimes exist.
#' @param obj A \code{longy_data} object.
#' @param regime Character vector of regime names, or NULL for all.
#' @return Character vector of validated regime names.
#' @noRd
.resolve_regimes <- function(obj, regime) {
  if (is.null(regime)) regime <- names(obj$regimes)
  if (length(regime) == 0) {
    stop("No regimes specified and none defined.", call. = FALSE)
  }
  bad <- setdiff(regime, names(obj$regimes))
  if (length(bad) > 0) {
    stop(sprintf("Regime(s) not found: %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  regime
}

#' Resolve learners argument for fit_* functions
#'
#' When \code{learners} is a named list (the \code{longy()} per-model format),
#' extracts the element matching \code{model_name} (falling back to
#' \code{$default}, then \code{c("SL.glm", "SL.mean")}). When it is a
#' character vector or NULL, returns it unchanged.
#' @param learners Character vector, named list, or NULL.
#' @param model_name One of "treatment", "censoring", "observation", "outcome".
#' @return Character vector of learner names, or NULL.
#' @noRd
.resolve_learners <- function(learners, model_name) {
  if (is.list(learners) && !is.null(names(learners))) {
    if (!is.null(learners[[model_name]])) return(learners[[model_name]])
    if (!is.null(learners$default))       return(learners$default)
    return(c("SL.glm", "SL.mean"))
  }
  learners
}

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
#' CV folds based on the effective sample size of an outcome vector.
#' For binary outcomes, effective sample size is based on minority class
#' prevalence. For continuous outcomes, effective sample size equals n.
#' @param y Outcome vector
#' @param binary Logical. If TRUE (default), use minority-class logic for
#'   effective sample size. If FALSE, use n directly (for continuous
#'   pseudo-outcomes in backward ICE steps).
#' @return List with n_eff, V, n, p_hat, n_rare
#' @noRd
.adaptive_cv_folds <- function(y, binary = TRUE) {
  n <- length(y)

  if (binary) {
    p_hat <- mean(y)
    n_rare <- n * min(p_hat, 1 - p_hat)
    n_eff <- min(n, 5 * n_rare)
  } else {
    p_hat <- NA_real_
    n_rare <- NA_real_
    n_eff <- n
  }

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
#' @param sl_fn Character. Which SuperLearner implementation to use:
#'   \code{"SuperLearner"} (default, sequential) or \code{"ffSL"}
#'   (future-factorial, parallelizes fold x algorithm combinations).
#' @param verbose Logical
#' @return List with predictions (numeric vector) and fit object
#' @noRd
.safe_sl <- function(Y, X, family = stats::binomial(),
                     learners = c("SL.glm", "SL.mean"),
                     cv_folds = 10L, obs_weights = NULL,
                     sl_fn = "SuperLearner",
                     context = "", verbose = FALSE) {
  X <- as.data.frame(X)

  # If SuperLearner is available and learners specified, try it
  if (!is.null(learners) && length(learners) > 0 &&
      requireNamespace("SuperLearner", quietly = TRUE)) {

    # Determine which SL function to call
    use_ffSL <- identical(sl_fn, "ffSL")
    if (use_ffSL && !requireNamespace("future.apply", quietly = TRUE)) {
      warning("future.apply not available; falling back to standard SuperLearner.",
              call. = FALSE)
      use_ffSL <- FALSE
    }

    sl_call_fn <- if (use_ffSL) .ffSL else SuperLearner::SuperLearner

    # Many SL wrappers (SL.earth, SL.nnet, SL.ranger) only handle
    # family$family == "gaussian" or "binomial". quasibinomial uses the
    # same link/variance so we can safely substitute for SL compatibility.
    sl_family <- family
    is_quasi <- identical(family$family, "quasibinomial")
    if (is_quasi) {
      sl_family <- stats::binomial()
    }

    # SuperLearner looks up learner and screening functions in env=;
    # use a combined environment so custom wrappers in globalenv() are
    # found alongside standard SL learners in the SuperLearner namespace.
    sl_env <- new.env(parent = asNamespace("SuperLearner"))
    for (nm in ls(globalenv(), pattern = "^SL\\.")) {
      assign(nm, get(nm, envir = globalenv()), envir = sl_env)
    }

    # When Y is continuous [0,1] (quasibinomial), some learners crash or
    # misbehave. Swap them for gaussian-family wrappers that truncate
    # predictions to the bounds after fitting.
    if (is_quasi) {
      q_lo <- 0.005
      q_hi <- 0.995

      if ("SL.xgboost" %in% learners) {
        xgb_reg_fn <- function(Y, X, newX, family, obsWeights, id, ...) {
          out <- SuperLearner::SL.xgboost(Y = Y, X = X, newX = newX,
                                           family = stats::gaussian(),
                                           obsWeights = obsWeights, id = id, ...)
          out$pred <- pmin(pmax(out$pred, q_lo), q_hi)
          out
        }
        assign("SL.xgboost.reg", xgb_reg_fn, envir = sl_env)
        learners[learners == "SL.xgboost"] <- "SL.xgboost.reg"
      }

      if ("SL.glmnet" %in% learners) {
        y_is_binary <- all(Y %in% c(0, 1))
        if (!y_is_binary) {
          message("SL.glmnet: switching to gaussian family for continuous ",
                  "[0,1] pseudo-outcomes (binomial requires integer Y). ",
                  "Predictions clipped to [", q_lo, ", ", q_hi, "].")
          glmnet_reg_fn <- function(Y, X, newX, family, obsWeights, id, ...) {
            out <- SuperLearner::SL.glmnet(Y = Y, X = X, newX = newX,
                                            family = stats::gaussian(),
                                            obsWeights = obsWeights, id = id, ...)
            out$pred <- pmin(pmax(out$pred, q_lo), q_hi)
            out
          }
          assign("SL.glmnet.reg", glmnet_reg_fn, envir = sl_env)
          learners[learners == "SL.glmnet"] <- "SL.glmnet.reg"
        }
      }

      if ("SL.ranger" %in% learners) {
        y_is_binary <- all(Y %in% c(0, 1))
        if (!y_is_binary) {
          message("SL.ranger: switching to gaussian family for continuous ",
                  "[0,1] pseudo-outcomes (tree-based predictions naturally bounded).")
          ranger_reg_fn <- function(Y, X, newX, family, obsWeights, id, ...) {
            out <- SuperLearner::SL.ranger(Y = Y, X = X, newX = newX,
                                            family = stats::gaussian(),
                                            obsWeights = obsWeights, id = id, ...)
            out$pred <- pmin(pmax(out$pred, q_lo), q_hi)
            out
          }
          assign("SL.ranger.reg", ranger_reg_fn, envir = sl_env)
          learners[learners == "SL.ranger"] <- "SL.ranger.reg"
        }
      }
    }

    fit <- tryCatch(
      {
        sl_args <- list(
          Y = Y, X = X, family = sl_family,
          SL.library = learners,
          cvControl = list(V = cv_folds),
          env = sl_env
        )
        if (!is.null(obs_weights)) {
          sl_args$obsWeights <- obs_weights
        }
        sl_fit <- do.call(sl_call_fn, sl_args)
        # Post-fit learner failure diagnostics
        cv_errs <- sl_fit$errorsInCVLibrary
        lib_errs <- sl_fit$errorsInLibrary
        all_errs <- as.logical(cv_errs) | as.logical(lib_errs)
        n_failed <- sum(all_errs)
        n_total <- length(all_errs)
        if (n_failed > 0) {
          failed_names <- sl_fit$libraryNames[all_errs]
          ctx <- if (nzchar(context)) paste0(" [", context, "]") else ""
          if (n_failed >= n_total * 0.5) {
            warning(sprintf(
              "%d/%d learner(s) failed%s: %s. SL is relying heavily on survivors.",
              n_failed, n_total, ctx, paste(failed_names, collapse = ", ")),
              call. = FALSE)
          } else if (verbose) {
            .vmsg("  %d/%d learner(s) failed%s: %s",
                  n_failed, n_total, ctx, paste(failed_names, collapse = ", "))
          }
        }

        list(
          predictions = as.numeric(sl_fit$SL.predict),
          fit = sl_fit,
          method = "SuperLearner",
          sl_risk = sl_fit$cvRisk,
          sl_coef = sl_fit$coef
        )
      },
      error = function(e) {
        ctx <- if (nzchar(context)) paste0(" [", context, "]") else ""
        warning(sprintf("SuperLearner failed%s (n=%d, event_rate=%.3f): %s. Falling back to glm.",
                        ctx, length(Y), mean(Y), e$message), call. = FALSE)
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

#' Predict from a .safe_sl() fit result on new data
#'
#' @param fit_result List returned by \code{.safe_sl()}.
#' @param newdata data.frame of new covariates.
#' @return Numeric vector of predicted values.
#' @noRd
.predict_from_fit <- function(fit_result, newdata) {
  newdata <- as.data.frame(newdata)
  if (fit_result$method == "SuperLearner") {
    as.numeric(stats::predict(fit_result$fit, newdata = newdata)$pred)
  } else if (fit_result$method == "glm") {
    as.numeric(stats::predict(fit_result$fit, newdata = newdata,
                               type = "response"))
  } else {
    rep(mean(fit_result$predictions), nrow(newdata))
  }
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
