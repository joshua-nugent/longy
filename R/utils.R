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
#' \code{$obj}). Also detects longy_data objects that lost their class
#' attribute during serialization (e.g. \code{saveRDS}/\code{readRDS})
#' and restores the class automatically.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object.
#' @return A \code{longy_data} object.
#' @noRd
.as_longy_data <- function(obj) {
  if (inherits(obj, "longy_data")) {
    .repair_longy_data(obj)
    return(obj)
  }
  if (inherits(obj, "longy_result") && !is.null(obj$obj)) {
    .repair_longy_data(obj$obj)
    return(obj$obj)
  }
  if (inherits(obj, "longy_results") && length(obj) > 0 &&
      !is.null(obj[[1]]$obj)) {
    .repair_longy_data(obj[[1]]$obj)
    return(obj[[1]]$obj)
  }
  # Detect longy_data by structure when class is lost (e.g. after readRDS)
  if (is.list(obj) && .is_longy_data_structure(obj)) {
    class(obj) <- "longy_data"
    .repair_longy_data(obj)
    return(obj)
  }
  stop("Expected a longy_data, longy_result, or longy_results object.",
       call. = FALSE)
}

#' Repair a longy_data object after deserialization
#'
#' Fixes stale data.table internal selfref pointers that are invalidated
#' by \code{saveRDS}/\code{readRDS}. Without this, the first \code{:=}
#' operation on a deserialized data.table triggers a warning and copy.
#'
#' @param obj A longy_data object (modified by reference where possible).
#' @return Invisible NULL.
#' @noRd
.repair_longy_data <- function(obj) {
  # Repair the main data.table
  if (data.table::is.data.table(obj$data)) {
    data.table::setalloccol(obj$data)
  }
  # Repair prediction data.tables inside fits
  for (mtype in c("treatment", "censoring", "observation", "outcome")) {
    for (rname in names(obj$fits[[mtype]])) {
      fit <- obj$fits[[mtype]][[rname]]
      if (is.list(fit)) {
        if (data.table::is.data.table(fit$predictions)) {
          data.table::setalloccol(fit$predictions)
        }
        # Censoring has nested per-cause lists
        if (mtype == "censoring" && is.list(fit)) {
          for (cvar in names(fit)) {
            if (is.list(fit[[cvar]]) &&
                data.table::is.data.table(fit[[cvar]]$predictions)) {
              data.table::setalloccol(fit[[cvar]]$predictions)
            }
          }
        }
      }
    }
  }
  # Repair weight data.tables
  for (rname in names(obj$weights)) {
    wt <- obj$weights[[rname]]
    if (is.list(wt) && data.table::is.data.table(wt$weights_dt)) {
      data.table::setalloccol(wt$weights_dt)
    }
  }
  invisible(NULL)
}

#' Check if a list has the structural signature of a longy_data object
#'
#' Used by \code{.as_longy_data()} to detect objects that lost their class
#' attribute during serialization.
#'
#' @param obj A list to check.
#' @return Logical.
#' @noRd
.is_longy_data_structure <- function(obj) {
  required <- c("data", "nodes", "fits", "meta")
  all(required %in% names(obj)) &&
    is.list(obj$nodes) &&
    is.list(obj$fits) &&
    is.list(obj$meta) &&
    !is.null(obj$nodes$id) &&
    !is.null(obj$nodes$time)
}

#' Restore and Repair a longy_data Object
#'
#' Validates a longy_data object, restoring the S3 class if it was lost
#' during serialization (e.g. \code{saveRDS}/\code{readRDS}) and repairing
#' internal data.table pointers. Call this after loading a saved longy_data
#' object to ensure all methods (\code{plot}, \code{print}, diagnostics)
#' work correctly.
#'
#' @param obj A \code{longy_data} object, or a list with longy_data structure
#'   that lost its class attribute during serialization.
#'
#' @return A repaired \code{longy_data} object with class restored.
#'
#' @examples
#' \dontrun{
#' # Save and reload a longy_data object
#' saveRDS(obj, "my_results.rds")
#' obj2 <- as_longy_data(readRDS("my_results.rds"))
#' plot(obj2)
#' weight_diagnostics(obj2)
#' }
#'
#' @export
as_longy_data <- function(obj) {
  .as_longy_data(obj)
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
          n_clip <- sum(out$pred < q_lo | out$pred > q_hi)
          if (n_clip > 0) {
            message(sprintf("SL.xgboost.reg: %d/%d predictions clipped to [%s, %s].",
                            n_clip, length(out$pred), q_lo, q_hi))
          }
          out$pred <- pmin(pmax(out$pred, q_lo), q_hi)
          out
        }
        assign("SL.xgboost.reg", xgb_reg_fn, envir = sl_env)
        learners[learners == "SL.xgboost"] <- "SL.xgboost.reg"
      }

      if ("SL.glmnet" %in% learners) {
        y_is_binary <- all(Y %in% c(0, 1))
        if (!y_is_binary) {
          glmnet_reg_fn <- function(Y, X, newX, family, obsWeights, id, ...) {
            out <- SuperLearner::SL.glmnet(Y = Y, X = X, newX = newX,
                                            family = stats::gaussian(),
                                            obsWeights = obsWeights, id = id, ...)
            n_clip <- sum(out$pred < q_lo | out$pred > q_hi)
            if (n_clip > 0) {
              message(sprintf("SL.glmnet.reg: %d/%d predictions clipped to [%s, %s].",
                              n_clip, length(out$pred), q_lo, q_hi))
            }
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
          # SL.ranger with gaussian returns a regression fit, but
          # predict.SL.ranger checks family$family == "binomial" (passed by
          # predict.SuperLearner) and tries pred[, "1"] on a vector.
          # Fix: call ranger directly and use a custom predict method.
          ranger_reg_fn <- function(Y, X, newX, family, obsWeights, id, ...) {
            if (!requireNamespace("ranger", quietly = TRUE))
              stop("ranger package required", call. = FALSE)
            if (is.matrix(X)) X <- data.frame(X)
            if (is.matrix(newX)) newX <- data.frame(newX)
            fit_obj <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                                       case.weights = obsWeights,
                                       write.forest = TRUE,
                                       num.threads = 1, verbose = FALSE)
            pred <- stats::predict(fit_obj, data = newX)$predictions
            n_clip <- sum(pred < q_lo | pred > q_hi)
            if (n_clip > 0) {
              message(sprintf("SL.ranger.reg: %d/%d predictions clipped to [%s, %s].",
                              n_clip, length(pred), q_lo, q_hi))
            }
            pred <- pmin(pmax(pred, q_lo), q_hi)
            # Store bounds in fit so predict method is self-contained
            fit <- list(object = fit_obj, bounds = c(q_lo, q_hi))
            class(fit) <- "SL.ranger.reg"
            list(pred = pred, fit = fit)
          }
          # Self-contained predict method (reads bounds from fit object,
          # no closure dependencies for robust serialization/dispatch)
          ranger_reg_pred_fn <- function(object, newdata, family, ...) {
            if (!requireNamespace("ranger", quietly = TRUE))
              stop("ranger package required", call. = FALSE)
            if (is.matrix(newdata)) newdata <- data.frame(newdata)
            pred <- stats::predict(object$object, data = newdata,
                            num.threads = 1)$predictions
            lo <- object$bounds[1]
            hi <- object$bounds[2]
            pmin(pmax(pred, lo), hi)
          }
          assign("SL.ranger.reg", ranger_reg_fn, envir = sl_env)
          # Register predict method in multiple locations for robust S3 dispatch:
          # 1. sl_env (so ffSL workers can find it)
          assign("predict.SL.ranger.reg", ranger_reg_pred_fn, envir = sl_env)
          # 2. globalenv (fallback for direct R sessions)
          assign("predict.SL.ranger.reg", ranger_reg_pred_fn,
                 envir = globalenv())
          # 3. S3 methods table (most robust for UseMethod dispatch)
          tryCatch(
            registerS3method("predict", "SL.ranger.reg",
                             ranger_reg_pred_fn,
                             envir = asNamespace("stats")),
            error = function(e) NULL
          )
          learners[learners == "SL.ranger"] <- "SL.ranger.reg"
        }
      }
    }

    # Capture warnings during SL call to surface learner error messages
    sl_warnings <- list()
    fit <- tryCatch(
      withCallingHandlers(
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
          # Always print learner failure info with error details
          err_msgs <- vapply(sl_warnings, function(w) conditionMessage(w),
                             character(1))
          # Extract the actual error messages (from SL's "Error in algorithm X" warnings)
          algo_errs <- grep("Error in algorithm", err_msgs, value = TRUE)
          err_detail <- if (length(algo_errs) > 0) {
            paste0("\n    ", paste(unique(algo_errs), collapse = "\n    "))
          } else {
            ""
          }
          msg <- sprintf(
            "%d/%d learner(s) failed%s: %s%s",
            n_failed, n_total, ctx,
            paste(failed_names, collapse = ", "),
            err_detail)
          if (n_failed >= n_total * 0.5) {
            warning(paste0(msg, "\n    SL is relying heavily on survivors."),
                    call. = FALSE)
          } else {
            message(msg)
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
      warning = function(w) {
        sl_warnings[[length(sl_warnings) + 1L]] <<- w
        # Muffle SL's own algorithm error warnings (we report them better)
        # Let all other warnings through normally
        msg <- conditionMessage(w)
        if (grepl("Error in algorithm|Coefficients already 0", msg)) {
          invokeRestart("muffleWarning")
        }
      }),
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

#' Format covariate listing for verbose output
#'
#' Prints baseline, time-varying, and lagged covariates as a compact message.
#' Cleans up internal lag column names for readability.
#' @param baseline Character vector of baseline covariate names
#' @param timevarying Character vector of time-varying covariate names
#' @param lag_covs Character vector of lag covariate column names
#' @param indent Character. Indentation prefix.
#' @noRd
.vmsg_covariates <- function(baseline, timevarying, lag_covs,
                              indent = "    ") {
  # Clean lag names: .longy_lag_A_1 -> A_lag1
  clean_lags <- if (length(lag_covs) > 0) {
    gsub("^\\.longy_lag_(.+)_(\\d+)$", "\\1_lag\\2", lag_covs)
  } else {
    character(0)
  }

  parts <- character(0)
  if (length(baseline) > 0) {
    parts <- c(parts, sprintf("baseline: %s", paste(baseline, collapse = ", ")))
  }
  if (length(timevarying) > 0) {
    parts <- c(parts, sprintf("time-varying: %s",
                               paste(timevarying, collapse = ", ")))
  }
  if (length(clean_lags) > 0) {
    parts <- c(parts, sprintf("lagged: %s", paste(clean_lags, collapse = ", ")))
  } else {
    parts <- c(parts, "lagged: (none)")
  }
  .vmsg("%s%s", indent, paste(parts, collapse = "; "))
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
