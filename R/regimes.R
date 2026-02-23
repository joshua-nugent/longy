#' Define an Intervention Regime
#'
#' Adds a treatment regime (static, dynamic, or stochastic) to a `longy_data` object.
#'
#' @param obj A `longy_data` object.
#' @param name Character. Name for this regime (must be unique).
#' @param static Integer (0 or 1). For static regimes: always treat (1) or never treat (0).
#' @param dynamic A function `f(data_row)` returning 0 or 1. For dynamic regimes.
#' @param stochastic A function `f(data_row)` returning P(A=1). For stochastic regimes.
#' @param description Character. Optional human-readable description.
#'
#' @return Modified `longy_data` object with the regime added.
#'
#' @details
#' Exactly one of `static`, `dynamic`, or `stochastic` must be provided.
#'
#' @examples
#' \dontrun{
#' obj <- obj |>
#'   define_regime(name = "always_treat", static = 1L) |>
#'   define_regime(name = "never_treat", static = 0L)
#' }
#'
#' @export
define_regime <- function(obj, name, static = NULL, dynamic = NULL,
                          stochastic = NULL, description = NULL) {
  obj <- .as_longy_data(obj)

  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("'name' must be a non-empty character string.", call. = FALSE)
  }

  if (name %in% names(obj$regimes)) {
    stop(sprintf("Regime '%s' already exists. Use a different name.", name),
         call. = FALSE)
  }

  # Exactly one of static/dynamic/stochastic

  n_specified <- sum(!is.null(static), !is.null(dynamic), !is.null(stochastic))
  if (n_specified != 1) {
    stop("Exactly one of 'static', 'dynamic', or 'stochastic' must be provided.",
         call. = FALSE)
  }

  if (!is.null(static)) {
    if (!static %in% c(0L, 1L, 0, 1)) {
      stop("'static' must be 0 or 1.", call. = FALSE)
    }
    regime <- list(
      name = name,
      type = "static",
      value = as.integer(static),
      description = description %||% sprintf("Static regime: A = %d at all times", static)
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

#' Evaluate a regime to produce counterfactual treatment values
#'
#' @param regime A regime list (element of obj$regimes)
#' @param data A data.table of subject-time rows
#' @return Numeric vector: 0/1 for static/dynamic, probabilities for stochastic
#' @noRd
.evaluate_regime <- function(regime, data) {
  n <- nrow(data)

  if (regime$type == "static") {
    rep(regime$value, n)
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
