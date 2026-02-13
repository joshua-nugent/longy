#' Create a longy_data Object for Longitudinal Causal Inference
#'
#' Constructs and validates a `longy_data` S3 object from a long-format dataset
#' (one row per subject-time). This is the entry point for the longy pipeline.
#'
#' @param data A data.frame or data.table in long format (one row per subject-time).
#' @param id Character. Column name for subject identifier.
#' @param time Character. Column name for time index (integer).
#' @param outcome Character. Column name for the outcome variable.
#' @param treatment Character. Column name for binary treatment (0/1).
#' @param censoring Character vector. Column name(s) for absorbing censoring
#'   indicators (1 = censored). Once censored, subject is gone. Can be NULL
#'   if no censoring. When multiple censoring sources are provided, each is
#'   modeled separately and their weights are multiplied. **Important**: if
#'   multiple censoring events can occur within the same interval, the data
#'   should be prepared so that only the first censoring event is recorded.
#'   The package does not model within-interval ordering of censoring sources.
#' @param observation Character. Column name for intermittent outcome
#'   measurement indicator (1 = observed). Subject can return after R=0.
#'   If NULL, assumes outcome is always observed when uncensored.
#' @param baseline Character vector. Column names for time-invariant covariates.
#' @param timevarying Character vector. Column names for time-varying covariates.
#' @param sampling_weights Character. Column name for external sampling/survey
#'   weights (e.g., to generalize from the study to a target population). Must
#'   be non-negative (zero allowed to exclude subjects) and constant within
#'   subject. These weights are used in nuisance model fitting (as observation
#'   weights in glm/SuperLearner), in computing marginal rates for stabilized
#'   weights, and multiply the final IPW weight. Follows ltmle's
#'   `observation.weights` convention. NULL if none.
#' @param outcome_type Character. One of `"binary"`, `"continuous"`, `"survival"`.
#' @param competing_risks Logical. Whether the outcome involves competing risks.
#' @param verbose Logical. Print progress messages.
#'
#' @return An S3 object of class `"longy_data"`.
#'
#' @examples
#' \dontrun{
#' obj <- longy_data(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C", observation = "R",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2")
#' )
#' }
#'
#' @export
longy_data <- function(data,
                       id,
                       time,
                       outcome,
                       treatment,
                       censoring = NULL,
                       observation = NULL,
                       baseline = character(0),
                       timevarying = character(0),
                       sampling_weights = NULL,
                       outcome_type = "binary",
                       competing_risks = FALSE,
                       verbose = TRUE) {

  # --- Convert to data.table ---
  dt <- data.table::as.data.table(data)

  # --- Validate column existence ---
  all_nodes <- c(id, time, outcome, treatment, censoring, observation,
                 baseline, timevarying, sampling_weights)
  missing_cols <- setdiff(all_nodes, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # --- Validate uniqueness of (id, time) ---
  if (anyDuplicated(dt, by = c(id, time)) > 0) {
    stop("Duplicate (id, time) pairs found. Each subject-time must be unique.",
         call. = FALSE)
  }

  # --- Validate treatment is binary {0, 1} ---
  a_vals <- dt[[treatment]]
  a_vals_nona <- a_vals[!is.na(a_vals)]
  if (!all(a_vals_nona %in% c(0L, 1L, 0, 1))) {
    stop("Treatment column must be binary {0, 1}.", call. = FALSE)
  }

  # --- Validate censoring columns are binary {0, 1} ---
  if (!is.null(censoring)) {
    for (cvar in censoring) {
      c_vals <- dt[[cvar]]
      c_vals_nona <- c_vals[!is.na(c_vals)]
      if (!all(c_vals_nona %in% c(0L, 1L, 0, 1))) {
        stop(sprintf("Censoring column '%s' must be binary {0, 1}.", cvar),
             call. = FALSE)
      }
    }
  }

  # --- Validate observation column is binary {0, 1} ---
  if (!is.null(observation)) {
    r_vals <- dt[[observation]]
    r_vals_nona <- r_vals[!is.na(r_vals)]
    if (!all(r_vals_nona %in% c(0L, 1L, 0, 1))) {
      stop("Observation column must be binary {0, 1}.", call. = FALSE)
    }
  }

  # --- Validate sampling weights ---
  if (!is.null(sampling_weights)) {
    sw_vals <- dt[[sampling_weights]]
    if (!is.numeric(sw_vals)) {
      stop("Sampling weights column must be numeric.", call. = FALSE)
    }
    if (any(sw_vals < 0, na.rm = TRUE)) {
      stop("Sampling weights must be non-negative.", call. = FALSE)
    }
    if (all(sw_vals == 0, na.rm = TRUE)) {
      stop("At least one sampling weight must be positive.", call. = FALSE)
    }
    # Must be constant within subject
    sw_unique <- dt[, list(nu = data.table::uniqueN(get(sampling_weights),
                                                     na.rm = TRUE)),
                    by = c(id)]
    if (any(sw_unique$nu > 1)) {
      bad_n <- sum(sw_unique$nu > 1)
      stop(sprintf(
        "Sampling weights vary within %d subject(s). Must be constant within subject.",
        bad_n), call. = FALSE)
    }
  }

  # --- Validate outcome_type ---
  outcome_type <- match.arg(outcome_type, c("binary", "continuous", "survival"))

  # --- Validate binary outcome values ---
  if (outcome_type == "binary") {
    y_vals <- dt[[outcome]]
    y_vals_nona <- y_vals[!is.na(y_vals)]
    if (length(y_vals_nona) > 0 && !all(y_vals_nona %in% c(0L, 1L, 0, 1))) {
      stop("Outcome must be binary {0, 1} when outcome_type = 'binary'.",
           call. = FALSE)
    }
  }

  # --- Validate survival outcome monotonicity ---
  if (outcome_type == "survival") {
    data.table::setkeyv(dt, c(id, time))
    dt[, .longy_check_surv := {
      y <- get(outcome)
      y_nona <- y[!is.na(y)]
      if (length(y_nona) < 2) TRUE
      else {
        first_event <- which(y_nona == 1)[1]
        if (is.na(first_event)) TRUE
        else all(y_nona[first_event:length(y_nona)] == 1)
      }
    }, by = c(id)]
    if (!all(dt$.longy_check_surv)) {
      bad_ids <- unique(dt[[id]][!dt$.longy_check_surv])
      dt[, .longy_check_surv := NULL]
      stop(sprintf(
        "Survival outcome must be absorbing (once Y=1, all subsequent Y=1). Violated for %d subject(s).",
        length(bad_ids)),
        call. = FALSE)
    }
    dt[, .longy_check_surv := NULL]
  }

  # --- Validate baseline covariates are constant within subject ---
  if (length(baseline) > 0) {
    data.table::setkeyv(dt, c(id, time))
    for (bvar in baseline) {
      n_unique <- dt[, list(nu = data.table::uniqueN(get(bvar), na.rm = TRUE)),
                     by = c(id)]
      if (any(n_unique$nu > 1)) {
        bad_n <- sum(n_unique$nu > 1)
        stop(sprintf(
          "Baseline covariate '%s' varies within %d subject(s). Baseline columns must be constant within subject.",
          bvar, bad_n),
          call. = FALSE)
        }
    }
  }

  # --- Convert factors/characters to dummies ---
  all_covariates <- c(baseline, timevarying)
  factor_cols <- character(0)
  for (col in all_covariates) {
    if (is.factor(dt[[col]]) || is.character(dt[[col]])) {
      factor_cols <- c(factor_cols, col)
    }
  }

  dummy_map <- list()  # old_name -> new_names
  if (length(factor_cols) > 0) {
    if (verbose) {
      .vmsg("Converting factor/character columns to dummies: %s",
            paste(factor_cols, collapse = ", "))
    }
    for (col in factor_cols) {
      vals <- dt[[col]]
      if (is.character(vals)) vals <- as.factor(vals)
      lvls <- levels(vals)
      new_names <- character(0)
      if (length(lvls) > 1) {
        for (i in seq_along(lvls)[-1]) {
          new_name <- paste0(col, "_", lvls[i])
          data.table::set(dt, j = new_name,
                          value = as.integer(vals == lvls[i]))
          new_names <- c(new_names, new_name)
        }
      }
      data.table::set(dt, j = col, value = NULL)
      dummy_map[[col]] <- new_names
    }

    # Update baseline/timevarying lists
    baseline <- .replace_with_dummies(baseline, dummy_map)
    timevarying <- .replace_with_dummies(timevarying, dummy_map)
  }

  # --- Ensure integer types for binary columns ---
  dt[, (treatment) := as.integer(get(treatment))]
  if (!is.null(censoring)) {
    for (cvar in censoring) {
      dt[, (cvar) := as.integer(get(cvar))]
    }
  }
  if (!is.null(observation)) {
    dt[, (observation) := as.integer(get(observation))]
  }

  # --- Set key ---
  data.table::setkeyv(dt, c(id, time))

  # --- Compute metadata ---
  time_vals <- sort(unique(dt[[time]]))
  ids <- unique(dt[[id]])

  nodes <- list(
    id = id,
    time = time,
    outcome = outcome,
    treatment = treatment,
    censoring = censoring,
    observation = observation,
    sampling_weights = sampling_weights,
    baseline = baseline,
    timevarying = timevarying,
    outcome_type = outcome_type,
    competing_risks = competing_risks
  )

  meta <- list(
    n_subjects = length(ids),
    n_obs = nrow(dt),
    n_times = length(time_vals),
    time_values = time_vals,
    max_time = max(time_vals),
    min_time = min(time_vals)
  )

  obj <- list(
    data = dt,
    nodes = nodes,
    regimes = list(),
    fits = list(treatment = NULL, censoring = list(), observation = NULL),
    weights = NULL,
    crossfit = list(enabled = FALSE, n_folds = NULL, fold_id = NULL),
    meta = meta
  )

  class(obj) <- "longy_data"

  if (verbose) {
    .vmsg("longy_data: %d subjects, %d observations, %d time points (%d-%d)",
          meta$n_subjects, meta$n_obs, meta$n_times,
          meta$min_time, meta$max_time)
  }

  obj
}

#' Replace factor column names with their dummy names
#' @noRd
.replace_with_dummies <- function(vars, dummy_map) {
  out <- character(0)
  for (v in vars) {
    if (v %in% names(dummy_map)) {
      out <- c(out, dummy_map[[v]])
    } else {
      out <- c(out, v)
    }
  }
  out
}

#' Set Up Cross-Fitting for a longy_data Object
#'
#' Assigns cross-fitting folds at the subject level.
#'
#' @param obj A `longy_data` object.
#' @param n_folds Integer. Number of folds.
#' @param fold_column Character. Name of an existing column with fold assignments.
#'   If provided, `n_folds` and `seed` are ignored.
#' @param seed Integer. Random seed for fold assignment.
#'
#' @return Modified `longy_data` object with crossfit information.
#' @export
set_crossfit <- function(obj, n_folds = 5L, fold_column = NULL, seed = NULL) {
  stopifnot(inherits(obj, "longy_data"))

  id_col <- obj$nodes$id

  if (!is.null(fold_column)) {
    # Use existing column
    if (!fold_column %in% names(obj$data)) {
      stop(sprintf("Fold column '%s' not found in data.", fold_column),
           call. = FALSE)
    }
    fold_dt <- unique(obj$data[, c(id_col, fold_column), with = FALSE])
    obj$crossfit <- list(
      enabled = TRUE,
      n_folds = data.table::uniqueN(fold_dt[[fold_column]]),
      fold_id = fold_column
    )
  } else {
    # Create folds at subject level
    ids <- unique(obj$data[[id_col]])
    n <- length(ids)
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(rep(seq_len(n_folds), length.out = n))
    fold_dt <- data.table::data.table(V1 = ids, .longy_fold = folds)
    data.table::setnames(fold_dt, "V1", id_col)

    obj$data <- merge(obj$data, fold_dt, by = id_col, all.x = TRUE)
    data.table::setkeyv(obj$data, c(id_col, obj$nodes$time))

    obj$crossfit <- list(
      enabled = TRUE,
      n_folds = n_folds,
      fold_id = ".longy_fold"
    )
  }

  obj
}

#' @export
print.longy_data <- function(x, ...) {
  cat("longy_data object\n")
  cat(sprintf("  Subjects:    %d\n", x$meta$n_subjects))
  cat(sprintf("  Time points: %d (%d to %d)\n",
              x$meta$n_times, x$meta$min_time, x$meta$max_time))
  cat(sprintf("  Observations: %d\n", x$meta$n_obs))
  cat(sprintf("  Outcome:     %s (%s)\n",
              x$nodes$outcome, x$nodes$outcome_type))
  cat(sprintf("  Treatment:   %s\n", x$nodes$treatment))
  if (!is.null(x$nodes$censoring)) {
    cat(sprintf("  Censoring:   %s\n",
                paste(x$nodes$censoring, collapse = ", ")))
  }
  if (!is.null(x$nodes$observation)) {
    cat(sprintf("  Observation: %s\n", x$nodes$observation))
  }
  if (!is.null(x$nodes$sampling_weights)) {
    cat(sprintf("  Sampling wt: %s\n", x$nodes$sampling_weights))
  }
  if (length(x$regimes) > 0) {
    cat(sprintf("  Regimes:     %s\n",
                paste(names(x$regimes), collapse = ", ")))
  }
  if (!is.null(x$weights)) {
    cat("  Weights:     computed\n")
  }
  invisible(x)
}

#' @export
summary.longy_data <- function(object, ...) {
  cat("=== longy_data summary ===\n\n")
  print(object)

  cat("\nNode columns:\n")
  cat(sprintf("  Baseline:     %s\n",
              if (length(object$nodes$baseline) > 0)
                paste(object$nodes$baseline, collapse = ", ")
              else "(none)"))
  cat(sprintf("  Time-varying: %s\n",
              if (length(object$nodes$timevarying) > 0)
                paste(object$nodes$timevarying, collapse = ", ")
              else "(none)"))

  # Treatment distribution
  a <- object$data[[object$nodes$treatment]]
  cat(sprintf("\nTreatment (%s): mean = %.3f\n",
              object$nodes$treatment, mean(a, na.rm = TRUE)))

  # Censoring summary
  if (!is.null(object$nodes$censoring)) {
    for (cvar in object$nodes$censoring) {
      c_vals <- object$data[[cvar]]
      cat(sprintf("Censoring (%s): %.1f%% censored overall\n",
                  cvar, 100 * mean(c_vals, na.rm = TRUE)))
    }
  }

  # Observation summary
  if (!is.null(object$nodes$observation)) {
    r_vals <- object$data[[object$nodes$observation]]
    cat(sprintf("Observation (%s): %.1f%% observed overall\n",
                object$nodes$observation, 100 * mean(r_vals, na.rm = TRUE)))
  }

  # Cross-fitting
  if (object$crossfit$enabled) {
    cat(sprintf("\nCross-fitting: %d folds (column: %s)\n",
                object$crossfit$n_folds, object$crossfit$fold_id))
  }

  invisible(object)
}
