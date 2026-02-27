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
#' @param censoring Character (length 1). Column name for a character/factor
#'   censoring status column. The column must contain \code{"uncensored"} for
#'   rows where the subject is not censored; any other value (e.g.
#'   \code{"censored"}, \code{"death"}, \code{"ltfu"}) is treated as a
#'   distinct censoring cause. Censoring is absorbing: once a subject has a
#'   non-\code{"uncensored"} value, they should have no further rows (or all
#'   subsequent rows should also be non-\code{"uncensored"}).
#'
#'   Each non-\code{"uncensored"} level is automatically decomposed into an
#'   internal binary indicator column (named \code{.cens_<cause>}) and modeled
#'   separately by \code{\link{fit_censoring}}. The internal column names are
#'   stored in \code{nodes$censoring}; the original column name and cause
#'   labels are in \code{nodes$censoring_col} and \code{nodes$censoring_levels}.
#'
#'   NAs are not allowed in the censoring column. If all values are
#'   \code{"uncensored"}, censoring is treated as NULL (no censoring).
#'   NULL if no censoring column exists.
#' @param observation Character. Column name for intermittent outcome
#'   measurement indicator (1 = observed). Subject can return after R=0.
#'   If NULL, auto-detected from NA values in the outcome column among
#'   uncensored rows. If NAs are found, an observation column is created
#'   automatically.
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
#' @param competing Character. Column name for a binary absorbing competing
#'   event indicator (1 = competing event occurred). Once D=1, must remain 1.
#'   Used with \code{outcome_type = "survival"} to estimate cause-specific
#'   cumulative incidence. Subjects who experience the competing event are
#'   excluded from the risk set with Q hard-coded to 0. NULL if no competing
#'   risks.
#' @param k Integer or \code{Inf}. Number of lagged time steps of covariate
#'   history to include as additional predictors in nuisance models. At each
#'   time index \code{i}, up to \code{min(k, i - 1)} lags are added for the
#'   treatment, outcome, and time-varying columns. \code{k = Inf} (default)
#'   uses all available history (matching ltmle's wide-format conditioning
#'   set). \code{k = 0} disables lag columns entirely.
#' @param verbose Logical. Print progress messages.
#'
#' @return An S3 object of class \code{"longy_data"}, a list with components:
#'   \describe{
#'     \item{data}{A \code{data.table} keyed on (id, time), including any
#'       internal columns created during validation (e.g. \code{.cens_*},
#'       \code{.obs}).}
#'     \item{nodes}{A list of column name mappings: \code{id}, \code{time},
#'       \code{outcome}, \code{treatment}, \code{censoring} (internal binary
#'       column names or NULL), \code{censoring_col} (original column name or
#'       NULL), \code{censoring_levels} (cause labels or NULL),
#'       \code{observation}, \code{sampling_weights}, \code{baseline},
#'       \code{timevarying}, \code{outcome_type}, \code{competing}.}
#'     \item{regimes}{Initially empty list; populated by
#'       \code{\link{define_regime}}.}
#'     \item{fits}{Initially empty; populated by \code{fit_treatment},
#'       \code{fit_censoring}, \code{fit_observation}.}
#'     \item{weights}{NULL until \code{\link{compute_weights}} is called.}
#'     \item{crossfit}{Cross-fitting configuration (enabled, n_folds,
#'       fold_id).}
#'     \item{meta}{Dataset metadata: n_subjects, n_obs, n_times, time_values,
#'       max_time, min_time.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Single censoring cause (C has values "uncensored" and "censored")
#' obj <- longy_data(
#'   data = sim_longy,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C", observation = "R",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2")
#' )
#' obj$nodes$censoring       # ".cens_censored"
#' obj$nodes$censoring_col   # "C"
#'
#' # Multiple censoring causes (C has values "uncensored", "death", "ltfu")
#' obj2 <- longy_data(
#'   data = my_data,
#'   id = "id", time = "time", outcome = "Y",
#'   treatment = "A", censoring = "C_status",
#'   baseline = c("W1", "W2"), timevarying = c("L1", "L2")
#' )
#' obj2$nodes$censoring       # c(".cens_death", ".cens_ltfu")
#' obj2$nodes$censoring_levels # c("death", "ltfu")
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
                       competing = NULL,
                       k = Inf,
                       verbose = TRUE) {

  # --- Convert to data.table ---
  dt <- data.table::as.data.table(data)

  # --- Validate column existence ---
  all_nodes <- c(id, time, outcome, treatment, censoring, observation,
                 baseline, timevarying, sampling_weights, competing)
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

  # --- Validate and decompose censoring column ---
  censoring_col <- NULL
  censoring_levels <- NULL
  if (!is.null(censoring)) {
    if (length(censoring) != 1L) {
      stop("censoring must be a single column name (or NULL).", call. = FALSE)
    }
    c_vals <- dt[[censoring]]
    if (!is.character(c_vals) && !is.factor(c_vals)) {
      stop(sprintf(
        "Censoring column '%s' must be character or factor (e.g. 'uncensored', 'censored').",
        censoring), call. = FALSE)
    }
    if (is.factor(c_vals)) c_vals <- as.character(c_vals)
    if (any(is.na(c_vals))) {
      stop(sprintf("Censoring column '%s' must not contain NAs.", censoring),
           call. = FALSE)
    }
    if (!"uncensored" %in% c_vals) {
      stop(sprintf(
        "Censoring column '%s' must contain 'uncensored' as a value.", censoring),
        call. = FALSE)
    }
    # Decompose into binary indicator columns
    censoring_col <- censoring
    all_levels <- sort(unique(c_vals))
    censoring_levels <- setdiff(all_levels, "uncensored")
    if (length(censoring_levels) == 0L) {
      # All rows are "uncensored" — treat as no censoring
      censoring_col <- NULL
      censoring_levels <- NULL
      censoring <- NULL
    } else {
      internal_cols <- character(length(censoring_levels))
      for (j in seq_along(censoring_levels)) {
        lvl <- censoring_levels[j]
        col_name <- paste0(".cens_", lvl)
        data.table::set(dt, j = col_name,
                        value = as.integer(c_vals == lvl))
        internal_cols[j] <- col_name
      }
      censoring <- internal_cols
    }
  }

  # --- Validate competing event column ---
  if (!is.null(competing)) {
    # Must be used with survival outcome
    if (outcome_type != "survival") {
      stop("Competing event column can only be used with outcome_type = 'survival'.",
           call. = FALSE)
    }
    # Must be binary {0, 1}
    d_vals <- dt[[competing]]
    d_vals_nona <- d_vals[!is.na(d_vals)]
    if (!all(d_vals_nona %in% c(0L, 1L, 0, 1))) {
      stop(sprintf("Competing event column '%s' must be binary {0, 1}.", competing),
           call. = FALSE)
    }
    # Must be absorbing: once D=1, all subsequent D=1
    data.table::setkeyv(dt, c(id, time))
    dt[, .longy_check_comp := {
      d <- get(competing)
      d_nona <- d[!is.na(d)]
      if (length(d_nona) < 2) TRUE
      else {
        first_comp <- which(d_nona == 1)[1]
        if (is.na(first_comp)) TRUE
        else all(d_nona[first_comp:length(d_nona)] == 1)
      }
    }, by = c(id)]
    if (!all(dt$.longy_check_comp)) {
      bad_ids <- unique(dt[[id]][!dt$.longy_check_comp])
      dt[, .longy_check_comp := NULL]
      stop(sprintf(
        "Competing event must be absorbing (once D=1, all subsequent D=1). Violated for %d subject(s).",
        length(bad_ids)),
        call. = FALSE)
    }
    dt[, .longy_check_comp := NULL]
    # Mutual exclusivity: no subject has both Y=1 and D=1 at the same time
    y_vals <- dt[[outcome]]
    d_vals <- dt[[competing]]
    both <- !is.na(y_vals) & !is.na(d_vals) & y_vals == 1 & d_vals == 1
    if (any(both)) {
      n_both <- sum(both)
      stop(sprintf(
        "Primary event (Y=1) and competing event (D=1) occur simultaneously in %d row(s). Events must be mutually exclusive.",
        n_both),
        call. = FALSE)
    }
  }

  # --- Survival LTCF: carry forward Y=1 after first event ---
  # For survival outcomes, NAs in Y after a Y=1 event represent absorption
  # (subject had the event and was removed), not intermittent missingness.
  # Carry forward Y=1 so that: (1) auto-detection of .obs doesn't fire on
  # these rows, (2) lag columns correctly show Y=1 for absorbed subjects,
  # (3) the data matches ltmle's expected format (LTCF).
  if (outcome_type == "survival") {
    data.table::setkeyv(dt, c(id, time))
    y_vals <- dt[[outcome]]
    # Find first event time per subject
    event_mask <- !is.na(y_vals) & as.numeric(y_vals) == 1
    if (any(event_mask)) {
      time_col_name <- time  # avoid data.table scope confusion
      event_dt <- dt[event_mask, list(.longy_first_ev = min(get(time_col_name))),
                      by = c(id)]
      dt[event_dt, .longy_first_ev := i..longy_first_ev, on = id]
      # Rows after event: NA outcome -> 1
      post_event_na <- !is.na(dt$.longy_first_ev) &
                       dt[[time_col_name]] > dt$.longy_first_ev &
                       is.na(dt[[outcome]])
      n_ltcf <- sum(post_event_na)
      if (n_ltcf > 0) {
        data.table::set(dt, which(post_event_na), outcome, 1L)
        if (verbose) {
          .vmsg("Survival LTCF: set %d post-event NA outcome rows to 1 (absorbing state).",
                n_ltcf)
        }
      }
      dt[, .longy_first_ev := NULL]
    }
  }

  # --- Auto-detect observation from outcome NAs ---
  if (is.null(observation)) {
    # Check for NAs in outcome among uncensored rows
    if (!is.null(censoring_col)) {
      censored <- dt[[censoring_col]] != "uncensored"
      uncensored_na <- !censored & is.na(dt[[outcome]])
    } else {
      uncensored_na <- is.na(dt[[outcome]])
    }
    if (any(uncensored_na)) {
      observation <- ".obs"
      dt[, .obs := as.integer(!is.na(get(outcome)))]
      if (verbose) {
        n_miss <- sum(uncensored_na)
        .vmsg("Auto-detected intermittent missingness: %d uncensored rows with NA outcome. Created observation column '.obs'.",
              n_miss)
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
  # Internal .cens_* columns already created as integer during decomposition
  if (!is.null(observation)) {
    dt[, (observation) := as.integer(get(observation))]
  }
  if (!is.null(competing)) {
    dt[, (competing) := as.integer(get(competing))]
  }

  # --- Set key ---
  data.table::setkeyv(dt, c(id, time))

  # --- Add lag columns for covariate history ---
  # Columns to lag: time-varying covariates, treatment, and outcome
  lag_vars <- character(0)
  lag_k <- 0
  if (!is.null(k) && (is.infinite(k) || k > 0)) {
    cols_to_lag <- unique(c(timevarying, treatment, outcome))
    if (length(cols_to_lag) > 0) {
      time_vals_pre <- sort(unique(dt[[time]]))
      n_times_pre <- length(time_vals_pre)
      max_lags <- if (is.infinite(k)) n_times_pre - 1L else min(as.integer(k), n_times_pre - 1L)
      if (max_lags > 0) {
        # First, create a LOCF-filled copy of each column to lag from.
        # This ensures that unobserved values (e.g., Y=NA when R=0) are
        # carried forward from the last observed value rather than becoming
        # NA or 0 in the lag columns. Original columns are NOT modified.
        locf_cols <- character(length(cols_to_lag))
        names(locf_cols) <- cols_to_lag
        for (col in cols_to_lag) {
          locf_name <- paste0(".longy_locf_", col)
          locf_cols[col] <- locf_name
          dt[, (locf_name) := get(col)]
          # LOCF within subject: fill NAs with last observed value
          dt[, (locf_name) := data.table::nafill(get(locf_name), type = "locf"),
             by = c(id)]
          # Any remaining NAs (no prior value at all) become 0
          data.table::setnafill(dt, fill = 0, cols = locf_name)
        }

        # Now create lag columns from the LOCF-filled copies
        for (col in cols_to_lag) {
          locf_name <- locf_cols[col]
          for (j in seq_len(max_lags)) {
            lag_name <- paste0(".longy_lag_", col, "_", j)
            dt[, (lag_name) := data.table::shift(get(locf_name), n = j,
                                                  type = "lag"),
               by = c(id)]
            # Structural NAs from shift (first j time points) become 0
            data.table::setnafill(dt, fill = 0, cols = lag_name)
          }
        }

        # Remove temporary LOCF columns
        locf_to_remove <- unname(locf_cols)
        dt[, (locf_to_remove) := NULL]
        lag_vars <- cols_to_lag
        lag_k <- if (is.infinite(k)) Inf else as.integer(k)
        if (verbose) {
          .vmsg("Added lag columns: %d variable(s) x %d max lag(s) = %d columns",
                length(cols_to_lag), max_lags, length(cols_to_lag) * max_lags)
        }
      }
    }
  }

  # --- Compute metadata ---
  time_vals <- sort(unique(dt[[time]]))
  ids <- unique(dt[[id]])

  nodes <- list(
    id = id,
    time = time,
    outcome = outcome,
    treatment = treatment,
    censoring = censoring,
    censoring_col = censoring_col,
    censoring_levels = censoring_levels,
    observation = observation,
    sampling_weights = sampling_weights,
    baseline = baseline,
    timevarying = timevarying,
    outcome_type = outcome_type,
    competing = competing,
    lag_vars = lag_vars,
    lag_k = lag_k
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
    fits = list(treatment = list(), censoring = list(), observation = list(), outcome = list()),
    weights = list(),
    results = list(),
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
  obj <- .as_longy_data(obj)

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
    if (!is.null(x$nodes$censoring_col)) {
      cat(sprintf("  Censoring:   %s (causes: %s)\n",
                  x$nodes$censoring_col,
                  paste(x$nodes$censoring_levels, collapse = ", ")))
    } else {
      cat(sprintf("  Censoring:   %s\n",
                  paste(x$nodes$censoring, collapse = ", ")))
    }
  }
  if (!is.null(x$nodes$observation)) {
    cat(sprintf("  Observation: %s\n", x$nodes$observation))
  }
  if (!is.null(x$nodes$competing)) {
    cat(sprintf("  Competing:   %s\n", x$nodes$competing))
  }
  if (!is.null(x$nodes$sampling_weights)) {
    cat(sprintf("  Sampling wt: %s\n", x$nodes$sampling_weights))
  }
  if (length(x$regimes) > 0) {
    cat(sprintf("  Regimes:     %s\n",
                paste(names(x$regimes), collapse = ", ")))
  }
  if (length(x$weights) > 0) {
    cat(sprintf("  Weights:     computed (%s)\n",
                paste(names(x$weights), collapse = ", ")))
  }
  # Show fitted models
  fitted_parts <- character(0)
  for (mtype in c("treatment", "censoring", "observation", "outcome")) {
    rnames <- names(x$fits[[mtype]])
    rnames <- rnames[vapply(x$fits[[mtype]], function(f) length(f) > 0, logical(1))]
    if (length(rnames) > 0) {
      fitted_parts <- c(fitted_parts,
                        sprintf("%s (%s)", mtype, paste(rnames, collapse = ", ")))
    }
  }
  if (length(fitted_parts) > 0) {
    cat(sprintf("  Fits:        %s\n", paste(fitted_parts, collapse = "; ")))
  }
  # Show results
  if (length(x$results) > 0) {
    cat(sprintf("  Results:     %s\n",
                paste(names(x$results), collapse = ", ")))
  }
  invisible(x)
}

#' @export
summary.longy_data <- function(object, ...) {
  nodes <- object$nodes
  meta <- object$meta

  # Always show metadata header (longy_data summary)
  cat("longy_data summary\n")
  cat(sprintf("  Outcome type: %s\n", nodes$outcome_type))
  cat(sprintf("  Subjects:     %d\n", meta$n_subjects))
  cat(sprintf("  Time points:  %d (%s)\n", meta$n_times,
              paste(meta$time_values, collapse = ", ")))
  cat(sprintf("  Treatment:    %s\n", nodes$treatment))
  if (!is.null(nodes$censoring_col)) {
    causes <- nodes$censoring_levels
    cat(sprintf("  Censoring:    %s (causes: %s)\n",
                nodes$censoring_col, paste(causes, collapse = ", ")))
  }
  if (!is.null(nodes$observation)) {
    cat(sprintf("  Observation:  %s\n", nodes$observation))
  }
  if (object$crossfit$enabled) {
    cat(sprintf("  Cross-fit:    %d folds\n", object$crossfit$n_folds))
  }

  # Show regimes if defined
  if (length(object$regimes) > 0) {
    reg_names <- names(object$regimes)
    cat(sprintf("  Regimes:      %s\n", paste(reg_names, collapse = ", ")))
  }

  # Show fitted models
  fitted_parts <- character(0)
  for (mtype in c("treatment", "censoring", "observation", "outcome")) {
    if (length(object$fits[[mtype]]) > 0) {
      rnames <- names(object$fits[[mtype]])
      fitted_parts <- c(fitted_parts,
                        sprintf("%s (%s)", mtype, paste(rnames, collapse = ", ")))
    }
  }
  if (length(fitted_parts) > 0) {
    cat(sprintf("  Fits:         %s\n", paste(fitted_parts, collapse = "; ")))
  }

  if (length(object$results) == 0) {
    cat("\n  No results yet. Run an estimator to see estimates.\n")
    return(invisible(object))
  }

  # Combine all results into a single table (mirrors plot.longy_data logic)
  all_est <- lapply(names(object$results), function(rname) {
    res <- object$results[[rname]]
    est <- as.data.frame(res$estimates)
    est_type <- if (!is.null(res$estimator)) res$estimator else "ipw"
    est$estimator <- switch(est_type,
                            gcomp = "G-comp", tmle = "TMLE",
                            unadjusted = "Unadjusted", "IPW")
    est$regime <- res$regime
    est
  })
  combined <- data.table::rbindlist(all_est, fill = TRUE)

  estimators <- unique(combined$estimator)
  regimes <- unique(combined$regime)
  times <- sort(unique(combined$time))
  has_ci <- "ci_lower" %in% names(combined) && "ci_upper" %in% names(combined)
  has_se <- "se" %in% names(combined)

  cat(sprintf("\nEstimator(s): %s | Regime(s): %s\n",
              paste(estimators, collapse = ", "),
              paste(regimes, collapse = ", ")))

  # Print one table per estimator (like plot facets)
  for (est_lab in estimators) {
    cat(sprintf("\n--- %s ---\n", est_lab))
    sub <- combined[combined$estimator == est_lab, ]

    if (length(regimes) == 1) {
      # Single regime: simple table
      rg_sub <- sub[sub$regime == regimes[1], ]
      .print_summary_table(rg_sub, has_se, has_ci)
    } else {
      # Multiple regimes: wide format with regime columns side by side
      .print_summary_table_wide(sub, regimes, times, has_se, has_ci)
    }
  }

  invisible(object)
}

#' Print a single-regime summary table
#' @noRd
.print_summary_table <- function(dt, has_se, has_ci) {
  # Check if SE/CI actually have non-NA values for this subset
  real_se <- has_se && any(!is.na(dt$se))
  real_ci <- has_ci && any(!is.na(dt$ci_lower))

  if (real_ci && real_se) {
    cat(sprintf("  %5s %10s %8s   %s\n", "time", "estimate", "se", "95% CI"))
    cat(sprintf("  %s\n", paste(rep("-", 48), collapse = "")))
    for (i in seq_len(nrow(dt))) {
      cat(sprintf("  %5d %10.4f %8.4f   [%7.4f, %7.4f]\n",
                  dt$time[i], dt$estimate[i], dt$se[i],
                  dt$ci_lower[i], dt$ci_upper[i]))
    }
  } else if (real_se) {
    cat(sprintf("  %5s %10s %8s\n", "time", "estimate", "se"))
    cat(sprintf("  %s\n", paste(rep("-", 28), collapse = "")))
    for (i in seq_len(nrow(dt))) {
      cat(sprintf("  %5d %10.4f %8.4f\n", dt$time[i], dt$estimate[i], dt$se[i]))
    }
  } else {
    cat(sprintf("  %5s %10s\n", "time", "estimate"))
    cat(sprintf("  %s\n", paste(rep("-", 18), collapse = "")))
    for (i in seq_len(nrow(dt))) {
      cat(sprintf("  %5d %10.4f\n", dt$time[i], dt$estimate[i]))
    }
  }
}

#' Print a multi-regime wide summary table
#' @noRd
.print_summary_table_wide <- function(dt, regimes, times, has_se, has_ci) {
  # Check if SE/CI actually have non-NA values for this subset
  real_se <- has_se && any(!is.na(dt$se))
  real_ci <- has_ci && any(!is.na(dt$ci_lower))

  # Build header: time | regime1 est (CI) | regime2 est (CI) | ...
  if (real_ci && real_se) {
    col_w <- 30  # width per regime column
    hdr <- sprintf("  %5s", "time")
    for (rg in regimes) hdr <- paste0(hdr, sprintf("  %*s", col_w, rg))
    cat(hdr, "\n")
    cat(sprintf("  %s\n", paste(rep("-", 5 + (col_w + 2) * length(regimes)),
                                collapse = "")))
    for (tt in times) {
      line <- sprintf("  %5d", tt)
      for (rg in regimes) {
        row <- dt[dt$time == tt & dt$regime == rg, ]
        if (nrow(row) == 1) {
          cell <- sprintf("%7.4f (%6.4f) [%7.4f, %7.4f]",
                          row$estimate, row$se, row$ci_lower, row$ci_upper)
          # Truncate/pad to col_w
          line <- paste0(line, sprintf("  %*s", col_w, cell))
        } else {
          line <- paste0(line, sprintf("  %*s", col_w, ""))
        }
      }
      cat(line, "\n")
    }
  } else {
    col_w <- 12
    hdr <- sprintf("  %5s", "time")
    for (rg in regimes) hdr <- paste0(hdr, sprintf("  %*s", col_w, rg))
    cat(hdr, "\n")
    cat(sprintf("  %s\n", paste(rep("-", 5 + (col_w + 2) * length(regimes)),
                                collapse = "")))
    for (tt in times) {
      line <- sprintf("  %5d", tt)
      for (rg in regimes) {
        row <- dt[dt$time == tt & dt$regime == rg, ]
        if (nrow(row) == 1) {
          line <- paste0(line, sprintf("  %*s", col_w,
                                       sprintf("%.4f", row$estimate)))
        } else {
          line <- paste0(line, sprintf("  %*s", col_w, ""))
        }
      }
      cat(line, "\n")
    }
  }
}
