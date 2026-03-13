#' Define an Intervention Regime
#'
#' Adds a treatment regime (static, dynamic, or stochastic) to a `longy_data` object.
#'
#' @param obj A `longy_data` object.
#' @param name Character. Name for this regime (must be unique).
#' @param static Treatment regime specification (one of three forms):
#'   \describe{
#'     \item{Scalar 0/1}{Constant regime — same treatment at every time point.
#'       E.g., \code{static = 1L} (always treat).}
#'     \item{Named numeric vector}{Time-varying regime — names are time values
#'       (as character), values are 0/1. E.g.,
#'       \code{static = c("1" = 1, "2" = 1, "3" = 0)}.
#'       All time values present in the data must have a defined value.}
#'     \item{Function}{A function \code{f(time_values)} receiving a numeric
#'       vector of time points, returning a 0/1 integer vector of the same
#'       length. E.g., \code{static = function(t) as.integer(t <= 3)}.}
#'   }
#' @param shifted Character. Column name in \code{obj$data} containing
#'   pre-computed counterfactual treatment values (0/1) for each subject-time
#'   row. This is the most flexible option — compute the regime however you
#'   want, store it in a column, and point longy at it. Internally handled
#'   like a dynamic regime (consistency tracked per subject).
#' @param dynamic A function `f(data_row)` returning 0 or 1. For dynamic regimes.
#' @param stochastic A function `f(data_row)` returning P(A=1). For stochastic regimes.
#' @param description Character. Optional human-readable description.
#'
#' @return Modified `longy_data` object with the regime added.
#'
#' @details
#' Exactly one of \code{static}, \code{dynamic}, or \code{stochastic} must be
#' provided.
#'
#' \strong{Static regimes} assign treatment based only on time — the same
#' sequence for every subject. They can be constant (always/never treat) or
#' time-varying (treat for the first K periods, then stop). Three input forms
#' are supported:
#' \itemize{
#'   \item \strong{Scalar}: \code{static = 1L} — treat at every time point.
#'   \item \strong{Named vector}: \code{static = c("0" = 1, "1" = 1, "2" = 0)}
#'     — names must match the time values in your data (as character strings).
#'     Every time value must be covered.
#'   \item \strong{Function}: \code{static = function(t) as.integer(t <= 5)}
#'     — receives a numeric vector of time points, returns 0/1 integer vector.
#'     Convenient when there are many time points.
#' }
#'
#' \strong{Shifted (pre-computed) regimes} are the most flexible option. The
#' user pre-computes the counterfactual treatment assignment as a 0/1 column
#' in the data, then points longy at it. This is useful for complex rules that
#' depend on subject history, multiple covariates, or external logic. Internally
#' handled like a dynamic regime.
#'
#' \strong{Dynamic regimes} assign treatment based on a subject's covariates
#' at each time point. The function receives a single-row data.table and must
#' return 0 or 1.
#'
#' \strong{Stochastic regimes} specify a probability of treatment. The function
#' receives a single-row data.table and must return a value between 0 and 1.
#'
#' @examples
#' \dontrun{
#' # --- Static: constant ---
#' # Always treat / never treat
#' obj <- obj |>
#'   define_regime(name = "always", static = 1L) |>
#'   define_regime(name = "never",  static = 0L)
#'
#' # --- Static: time-varying (named vector) ---
#' # Treat for 6 months, then stop (time values 0-11)
#' obj <- define_regime(obj, name = "six_months",
#'   static = c("0"=1, "1"=1, "2"=1, "3"=1, "4"=1, "5"=1,
#'              "6"=0, "7"=0, "8"=0, "9"=0, "10"=0, "11"=0))
#'
#' # --- Static: time-varying (function) ---
#' # Same regime, but cleaner when there are many time points
#' obj <- define_regime(obj, name = "six_months_fn",
#'   static = function(t) as.integer(t <= 5))
#'
#' # Treat only at even-numbered visits
#' obj <- define_regime(obj, name = "alternating",
#'   static = function(t) as.integer(t %% 2 == 0))
#'
#' # --- Shifted (pre-computed) ---
#' # Compute the regime yourself, store in a column, point longy at it.
#' # Useful for complex rules that are hard to express as a function.
#' my_data$A_cf <- ifelse(my_data$time <= 5 & my_data$L1 > 0, 1L, 0L)
#' obj <- define_regime(obj, name = "custom_rule", shifted = "A_cf")
#'
#' # --- Dynamic ---
#' # Treat if a covariate exceeds a threshold
#' obj <- define_regime(obj, name = "treat_high_bp",
#'   dynamic = function(row) as.integer(row[["blood_pressure"]] > 140))
#'
#' # --- Stochastic ---
#' # 70% probability of treatment at every time point
#' obj <- define_regime(obj, name = "prob70",
#'   stochastic = function(row) 0.7)
#' }
#'
#' @export
define_regime <- function(obj, name, static = NULL, shifted = NULL,
                          dynamic = NULL, stochastic = NULL,
                          description = NULL) {
  obj <- .as_longy_data(obj)

  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("'name' must be a non-empty character string.", call. = FALSE)
  }

  if (name %in% names(obj$regimes)) {
    stop(sprintf("Regime '%s' already exists. Use a different name.", name),
         call. = FALSE)
  }

  # Exactly one of static/shifted/dynamic/stochastic

  n_specified <- sum(!is.null(static), !is.null(shifted),
                     !is.null(dynamic), !is.null(stochastic))
  if (n_specified != 1) {
    stop("Exactly one of 'static', 'shifted', 'dynamic', or 'stochastic' must be provided.",
         call. = FALSE)
  }

  if (!is.null(static)) {
    if (is.function(static)) {
      # Function f(time_values) → 0/1 vector
      regime <- list(
        name = name,
        type = "static",
        value = static,
        time_varying = TRUE,
        description = description %||% "Static time-varying regime (function)"
      )
    } else if (is.numeric(static) && length(static) == 1 && is.null(names(static))) {
      # Scalar 0/1
      if (!static %in% c(0L, 1L, 0, 1))
        stop("'static' scalar must be 0 or 1.", call. = FALSE)
      regime <- list(
        name = name,
        type = "static",
        value = as.integer(static),
        time_varying = FALSE,
        description = description %||% sprintf("Static regime: A = %d at all times", static)
      )
    } else if (is.numeric(static) && length(static) > 0 && !is.null(names(static))) {
      # Named vector keyed by time
      if (!all(static %in% c(0, 1)))
        stop("All values in 'static' vector must be 0 or 1.", call. = FALSE)
      if (any(names(static) == "" | is.na(names(static))))
        stop("All elements of 'static' vector must be named (names = time values).",
             call. = FALSE)
      regime <- list(
        name = name,
        type = "static",
        value = setNames(as.integer(static), names(static)),
        time_varying = TRUE,
        description = description %||% "Static time-varying regime (vector)"
      )
    } else {
      stop("'static' must be 0/1 (scalar), a named numeric vector of 0/1, or a function.",
           call. = FALSE)
    }
  } else if (!is.null(shifted)) {
    if (!is.character(shifted) || length(shifted) != 1) {
      stop("'shifted' must be a single column name (character string).", call. = FALSE)
    }
    if (!shifted %in% names(obj$data)) {
      stop(sprintf("Shifted column '%s' not found in data.", shifted), call. = FALSE)
    }
    col_vals <- obj$data[[shifted]]
    if (!is.numeric(col_vals) && !is.integer(col_vals)) {
      stop(sprintf("Shifted column '%s' must be numeric (0/1 values).", shifted),
           call. = FALSE)
    }
    bad_vals <- unique(col_vals[!col_vals %in% c(0L, 1L, 0, 1) & !is.na(col_vals)])
    if (length(bad_vals) > 0) {
      stop(sprintf("Shifted column '%s' contains non-0/1 values: %s",
                   shifted, paste(head(bad_vals, 5), collapse = ", ")),
           call. = FALSE)
    }
    regime <- list(
      name = name,
      type = "shifted",
      value = shifted,
      description = description %||% sprintf("Pre-computed regime from column '%s'", shifted)
    )
  } else if (!is.null(dynamic)) {
    if (!is.function(dynamic)) {
      stop("'dynamic' must be a function.", call. = FALSE)
    }
    regime <- list(
      name = name,
      type = "dynamic",
      value = dynamic,
      description = description %||% "Dynamic regime"
    )
  } else {
    if (!is.function(stochastic)) {
      stop("'stochastic' must be a function.", call. = FALSE)
    }
    regime <- list(
      name = name,
      type = "stochastic",
      value = stochastic,
      description = description %||% "Stochastic regime"
    )
  }

  obj$regimes[[name]] <- regime
  obj
}

#' Resolve static regime values for given time points
#'
#' Handles all three static forms: scalar, named vector, function.
#'
#' @param value Regime value (scalar, named vector, or function)
#' @param time_values Numeric vector of time points to resolve
#' @return Integer vector of 0/1, same length as \code{time_values}
#' @noRd
.resolve_static_at_time <- function(value, time_values) {
  if (is.function(value)) {
    result <- value(time_values)
  } else if (length(value) == 1 && is.null(names(value))) {
    # Scalar — same value at all times
    return(rep(as.integer(value), length(time_values)))
  } else {
    # Named vector — look up by time
    result <- value[as.character(time_values)]
    if (anyNA(result)) {
      bad <- unique(time_values[is.na(result)])
      stop(sprintf("Static regime has no value defined for time(s): %s",
                   paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  result <- as.integer(result)
  if (!all(result %in% c(0L, 1L))) {
    stop("Static regime values must all be 0 or 1.", call. = FALSE)
  }
  result
}

#' Evaluate a regime to produce counterfactual treatment values
#'
#' @param regime A regime list (element of obj$regimes)
#' @param data A data.table of subject-time rows
#' @param time_col Character. Name of the time column in \code{data}. Required
#'   for time-varying static regimes; ignored for dynamic/stochastic.
#' @return Numeric vector: 0/1 for static/dynamic, probabilities for stochastic
#' @noRd
.evaluate_regime <- function(regime, data, time_col = NULL) {
  n <- nrow(data)

  if (regime$type == "static") {
    if (isTRUE(regime$time_varying)) {
      if (is.null(time_col))
        stop("time_col required for time-varying static regimes.", call. = FALSE)
      .resolve_static_at_time(regime$value, data[[time_col]])
    } else {
      rep(regime$value, n)
    }
  } else if (regime$type == "shifted") {
    col <- regime$value
    if (!col %in% names(data))
      stop(sprintf("Shifted column '%s' not found in data.", col), call. = FALSE)
    as.integer(data[[col]])
  } else if (regime$type == "dynamic") {
    vapply(seq_len(n), function(i) {
      tryCatch(
        regime$value(data[i, ]),
        error = function(e) stop(sprintf(
          "Regime '%s' evaluation failed at row %d: %s",
          regime$name, i, e$message
        ), call. = FALSE)
      )
    }, numeric(1))
  } else if (regime$type == "stochastic") {
    vapply(seq_len(n), function(i) {
      tryCatch(
        regime$value(data[i, ]),
        error = function(e) stop(sprintf(
          "Regime '%s' evaluation failed at row %d: %s",
          regime$name, i, e$message
        ), call. = FALSE)
      )
    }, numeric(1))
  } else {
    stop(sprintf("Unknown regime type: %s", regime$type), call. = FALSE)
  }
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
