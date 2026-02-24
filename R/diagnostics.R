#' Weight Diagnostics
#'
#' Summarizes weight distributions by time point, including mean, median,
#' maximum, quantiles, and effective sample size.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with weights computed. For \code{longy_results}, uses the first
#'   element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param by_time Logical. If TRUE (default), summarize by time point.
#'
#' @return A data.table with weight summary statistics.
#' @export
weight_diagnostics <- function(obj, regime = NULL, by_time = TRUE) {
  obj <- .as_longy_data(obj)
  if (length(obj$weights) == 0) {
    stop("Weights not computed. Run compute_weights() first.", call. = FALSE)
  }
  if (is.null(regime)) regime <- names(obj$weights)[1]
  if (is.null(obj$weights[[regime]])) {
    stop(sprintf("No weights found for regime '%s'.", regime), call. = FALSE)
  }

  w_dt <- obj$weights[[regime]]$weights_dt

  if (by_time) {
    diag <- w_dt[, list(
      n = .N,
      mean_weight = round(mean(.final_weight), 2),
      median_weight = round(stats::median(.final_weight), 2),
      sd_weight = round(stats::sd(.final_weight), 2),
      min_weight = round(min(.final_weight), 2),
      max_weight = round(max(.final_weight), 2),
      p01 = round(stats::quantile(.final_weight, 0.01), 2),
      p05 = round(stats::quantile(.final_weight, 0.05), 2),
      p95 = round(stats::quantile(.final_weight, 0.95), 2),
      p99 = round(stats::quantile(.final_weight, 0.99), 2),
      ess = round(.ess(.final_weight), 1),
      cv = round(stats::sd(.final_weight) / mean(.final_weight), 2),
      mean_sw_a = round(mean(.sw_a), 2),
      mean_sw_c = round(mean(.sw_c), 2),
      mean_csw_ac = round(mean(.csw_ac), 2),
      mean_sw_r = round(mean(.sw_r), 2)
    ), by = .time]
    data.table::setnames(diag, ".time", "time")
    data.table::setkey(diag, time)
  } else {
    diag <- w_dt[, list(
      n = .N,
      mean_weight = round(mean(.final_weight), 2),
      median_weight = round(stats::median(.final_weight), 2),
      sd_weight = round(stats::sd(.final_weight), 2),
      min_weight = round(min(.final_weight), 2),
      max_weight = round(max(.final_weight), 2),
      p01 = round(stats::quantile(.final_weight, 0.01), 2),
      p99 = round(stats::quantile(.final_weight, 0.99), 2),
      ess = round(.ess(.final_weight), 1),
      cv = round(stats::sd(.final_weight) / mean(.final_weight), 2)
    )]
  }

  diag
}

#' Positivity Diagnostics
#'
#' Identifies observations with extreme propensity scores that may indicate
#' positivity violations.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with treatment model fit. For \code{longy_results}, uses the first
#'   element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param threshold Numeric. Flag observations with predicted probabilities
#'   below this value (or above 1 - threshold). Default 0.025.
#'
#' @return A data.table with flagged observations and their propensity scores.
#' @export
positivity_diagnostics <- function(obj, regime = NULL, threshold = 0.025) {
  obj <- .as_longy_data(obj)
  if (length(obj$fits$treatment) == 0) {
    stop("Treatment model not fit. Run fit_treatment() first.", call. = FALSE)
  }
  if (is.null(regime)) regime <- names(obj$fits$treatment)[1]
  if (is.null(obj$fits$treatment[[regime]])) {
    stop(sprintf("No treatment fit found for regime '%s'.", regime), call. = FALSE)
  }

  gA <- obj$fits$treatment[[regime]]$predictions
  id_col <- obj$nodes$id

  # Flag extreme propensity scores
  flagged <- gA[gA$.p_a < threshold | gA$.p_a > (1 - threshold), ]

  if (nrow(flagged) == 0) {
    message("No positivity violations detected at threshold ", threshold)
    return(invisible(data.table::data.table()))
  }

  # Summary by time
  summary_dt <- flagged[, list(
    n_flagged = .N,
    min_p_a = min(.p_a),
    max_p_a = max(.p_a),
    mean_p_a = mean(.p_a)
  ), by = .time]
  data.table::setnames(summary_dt, ".time", "time")
  data.table::setkey(summary_dt, time)

  message(sprintf("Positivity: %d observations flagged (threshold = %.3f) across %d time points",
                  nrow(flagged), threshold, nrow(summary_dt)))

  summary_dt
}

#' SuperLearner Diagnostics
#'
#' Extracts SuperLearner fit information from all fitted nuisance models,
#' including method used (SL, glm, marginal), CV risk, and learner coefficients.
#'
#' @param obj A \code{longy_data}, \code{longy_result}, or \code{longy_results}
#'   object with at least one model fit. For \code{longy_results}, uses the
#'   first element.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param model Character. Which model(s) to include: \code{"all"} (default),
#'   \code{"treatment"}, \code{"censoring"}, \code{"observation"}, or
#'   \code{"outcome"}.
#'
#' @return A data.table with columns: \code{model} (model type),
#'   \code{submodel} (censoring cause or NA), \code{time} (time point),
#'   \code{target_time} (target time for outcome models, or NA),
#'   \code{method} (fitting method used), \code{sl_risk} (CV risk of the
#'   SuperLearner, or NA), \code{sl_coef} (list-column of learner weights,
#'   or NULL).
#' @export
sl_diagnostics <- function(obj, regime = NULL, model = "all") {
  obj <- .as_longy_data(obj)

  model <- match.arg(model, c("all", "treatment", "censoring",
                               "observation", "outcome"))
  include <- if (model == "all") {
    c("treatment", "censoring", "observation", "outcome")
  } else {
    model
  }

  # Resolve regime: default to first available
  if (is.null(regime)) {
    available <- names(obj$fits$treatment)
    if (length(available) == 0) available <- names(obj$fits$outcome)
    if (length(available) == 0) available <- names(obj$fits$observation)
    regime <- available[1]
  }

  rows <- list()

  # If no regime resolved (no fits at all), skip to empty result
  if (is.na(regime) || is.null(regime)) {
    message("No model fits found.")
    return(invisible(data.table::data.table(
      model = character(), submodel = character(), time = integer(),
      target_time = integer(),
      method = character(), sl_risk = numeric(), sl_coef = list()
    )))
  }

  # Treatment
  trt_fit <- obj$fits$treatment[[regime]]
  if ("treatment" %in% include && !is.null(trt_fit$sl_info)) {
    for (entry in trt_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "treatment",
        submodel = NA_character_,
        time = entry$time,
        target_time = NA_integer_,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  # Censoring (one sub-list per cause)
  cens_fits <- obj$fits$censoring[[regime]]
  if ("censoring" %in% include && length(cens_fits) > 0) {
    for (cvar in names(cens_fits)) {
      sl_info <- cens_fits[[cvar]]$sl_info
      if (is.null(sl_info)) next
      for (entry in sl_info) {
        rows[[length(rows) + 1L]] <- list(
          model = "censoring",
          submodel = cvar,
          time = entry$time,
          target_time = NA_integer_,
          method = entry$method,
          sl_risk = entry$sl_risk %||% NA_real_,
          sl_coef = list(entry$sl_coef)
        )
      }
    }
  }

  # Observation
  obs_fit <- obj$fits$observation[[regime]]
  if ("observation" %in% include && !is.null(obs_fit$sl_info)) {
    for (entry in obs_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "observation",
        submodel = NA_character_,
        time = entry$time,
        target_time = NA_integer_,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  # Outcome
  out_fit <- obj$fits$outcome[[regime]]
  if ("outcome" %in% include && !is.null(out_fit$sl_info)) {
    for (entry in out_fit$sl_info) {
      rows[[length(rows) + 1L]] <- list(
        model = "outcome",
        submodel = NA_character_,
        time = entry$time,
        target_time = entry$target_time %||% NA_integer_,
        method = entry$method,
        sl_risk = entry$sl_risk %||% NA_real_,
        sl_coef = list(entry$sl_coef)
      )
    }
  }

  if (length(rows) == 0) {
    message("No model fits found.")
    return(invisible(data.table::data.table(
      model = character(), submodel = character(), time = integer(),
      target_time = integer(),
      method = character(), sl_risk = numeric(), sl_coef = list()
    )))
  }

  data.table::rbindlist(rows)
}

#' Plot SuperLearner Diagnostics
#'
#' Visualizes SuperLearner fit information across time points and nuisance
#' models. Requires ggplot2.
#'
#' @param obj A \code{longy_data} object, or a data.table from
#'   \code{sl_diagnostics()}.
#' @param regime Character. Regime name (passed to \code{sl_diagnostics()} if
#'   \code{obj} is a \code{longy_data}).
#' @param model Character. Which model(s) to include: \code{"all"} (default),
#'   \code{"treatment"}, \code{"censoring"}, \code{"observation"}, or
#'   \code{"outcome"}.
#' @param type Character. Plot type:
#'   \describe{
#'     \item{\code{"weights"}}{Stacked bar chart of ensemble weights by time,
#'       faceted by nuisance model. Shows which learners dominate and how the
#'       ensemble composition changes over time.}
#'     \item{\code{"risk"}}{Per-learner CV risk across time, faceted by nuisance
#'       model. Highlights where model fit degrades.}
#'     \item{\code{"heatmap"}}{Heatmap of learner weights (time x learner),
#'       faceted by nuisance model. Scales better than stacked bars with many
#'       learners or time points.}
#'     \item{\code{"method"}}{Tile plot showing which fitting method
#'       (SuperLearner, glm, marginal) was used at each time point.}
#'   }
#' @param drop_zero Logical. If TRUE (default), omit learners with zero weight
#'   from weight-based plots (\code{"weights"} and \code{"heatmap"}).
#'
#' @return A ggplot2 object, or NULL invisibly if there is nothing to plot.
#' @export
plot_sl_diagnostics <- function(obj, regime = NULL, model = "all",
                                type = c("weights", "risk", "heatmap", "method"),
                                drop_zero = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_sl_diagnostics().", call. = FALSE)
  }

  type <- match.arg(type)

  # Accept longy_data or pre-computed data.table
  if (inherits(obj, c("longy_data", "longy_result"))) {
    diag <- sl_diagnostics(obj, regime = regime, model = model)
  } else if (is.data.frame(obj)) {
    diag <- data.table::as.data.table(data.table::copy(obj))
  } else {
    stop("obj must be a longy_data object or sl_diagnostics() output.",
         call. = FALSE)
  }

  if (nrow(diag) == 0) {
    message("No SuperLearner diagnostics to plot.")
    return(invisible(NULL))
  }

  # Create human-readable facet labels
  diag[, facet_label := .make_sl_facet_label(model, submodel, target_time),
       by = seq_len(nrow(diag))]

  switch(type,
    weights = .plot_sl_weights(diag, drop_zero),
    risk    = .plot_sl_risk(diag),
    heatmap = .plot_sl_heatmap(diag, drop_zero),
    method  = .plot_sl_method(diag)
  )
}


# -- Internal helpers for plot_sl_diagnostics --------------------------------

#' Build a human-readable facet label
#' @noRd
.make_sl_facet_label <- function(model, submodel, target_time) {
  label <- model
  if (!is.na(submodel)) {
    label <- paste0(label, ": ", sub("^\\.cens_", "", submodel))
  }
  if (!is.na(target_time)) {
    label <- paste0(label, " (t=", target_time, ")")
  }
  label
}

#' Unnest sl_coef list-column into long format
#' @return data.table with facet_label, time, learner, weight; or NULL
#' @noRd
.unnest_sl_coef <- function(diag) {
  has_coef <- vapply(diag$sl_coef, function(x) {
    !is.null(x) && length(x) > 0
  }, logical(1))
  diag_sl <- diag[has_coef, ]
  if (nrow(diag_sl) == 0) return(NULL)

  data.table::rbindlist(lapply(seq_len(nrow(diag_sl)), function(i) {
    coefs <- diag_sl$sl_coef[[i]]
    data.table::data.table(
      facet_label = diag_sl$facet_label[i],
      time = diag_sl$time[i],
      learner = names(coefs),
      weight = as.numeric(coefs)
    )
  }))
}

#' Unnest sl_risk into long format (per-learner CV risk)
#' @return data.table with facet_label, time, learner, cv_risk; or NULL
#' @noRd
.unnest_sl_risk <- function(diag) {
  # sl_risk may be a plain numeric column (all NA) or a list-column
  if (!is.list(diag$sl_risk)) {
    if (all(is.na(diag$sl_risk))) return(NULL)
    diag_sl <- diag[!is.na(sl_risk)]
    return(data.table::data.table(
      facet_label = diag_sl$facet_label,
      time = diag_sl$time,
      learner = "SuperLearner",
      cv_risk = diag_sl$sl_risk
    ))
  }

  has_risk <- vapply(diag$sl_risk, function(x) {
    is.numeric(x) && length(x) > 0 && !all(is.na(x))
  }, logical(1))
  diag_sl <- diag[has_risk, ]
  if (nrow(diag_sl) == 0) return(NULL)

  data.table::rbindlist(lapply(seq_len(nrow(diag_sl)), function(i) {
    risks <- diag_sl$sl_risk[[i]]
    if (is.null(names(risks))) {
      data.table::data.table(
        facet_label = diag_sl$facet_label[i],
        time = diag_sl$time[i],
        learner = "SuperLearner",
        cv_risk = as.numeric(risks)
      )
    } else {
      data.table::data.table(
        facet_label = diag_sl$facet_label[i],
        time = diag_sl$time[i],
        learner = names(risks),
        cv_risk = as.numeric(risks)
      )
    }
  }))
}

#' Stacked bar chart of ensemble weights
#' @noRd
.plot_sl_weights <- function(diag, drop_zero) {
  long <- .unnest_sl_coef(diag)
  if (is.null(long) || nrow(long) == 0) {
    message("No SuperLearner fits found (all glm/marginal).")
    return(invisible(NULL))
  }

  if (drop_zero) long <- long[weight > 0]
  if (nrow(long) == 0) {
    message("All learner weights are zero.")
    return(invisible(NULL))
  }

  ggplot2::ggplot(long,
    ggplot2::aes(x = factor(time), y = weight, fill = learner)) +
    ggplot2::geom_col(position = "stack", width = 0.8) +
    ggplot2::facet_wrap(~facet_label, scales = "free_x") +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    ggplot2::labs(x = "Time", y = "Ensemble Weight", fill = "Learner",
                  title = "SuperLearner Ensemble Weights") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

#' Per-learner CV risk across time
#' @noRd
.plot_sl_risk <- function(diag) {
  long <- .unnest_sl_risk(diag)
  if (is.null(long) || nrow(long) == 0) {
    message("No CV risk information available.")
    return(invisible(NULL))
  }

  ggplot2::ggplot(long,
    ggplot2::aes(x = time, y = cv_risk, colour = learner)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~facet_label, scales = "free") +
    ggplot2::labs(x = "Time", y = "CV Risk", colour = "Learner",
                  title = "Per-Learner Cross-Validated Risk") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "bottom")
}

#' Heatmap of learner weights
#' @noRd
.plot_sl_heatmap <- function(diag, drop_zero) {
  long <- .unnest_sl_coef(diag)
  if (is.null(long) || nrow(long) == 0) {
    message("No SuperLearner fits found (all glm/marginal).")
    return(invisible(NULL))
  }

  if (drop_zero) {
    # Keep learners that have non-zero weight at ANY time, then show all their
    # time points (including zeros) for a complete picture
    active <- unique(long[weight > 0, learner])
    long <- long[learner %in% active]
  }
  if (nrow(long) == 0) {
    message("All learner weights are zero.")
    return(invisible(NULL))
  }

  ggplot2::ggplot(long,
    ggplot2::aes(x = factor(time), y = learner, fill = weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::facet_wrap(~facet_label, scales = "free") +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2166AC",
                                 limits = c(0, 1)) +
    ggplot2::labs(x = "Time", y = "Learner", fill = "Weight",
                  title = "SuperLearner Weight Heatmap") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Tile plot of fitting method by model and time
#' @noRd
.plot_sl_method <- function(diag) {
  method_colours <- c(
    "SuperLearner" = "#2166AC",
    "glm"          = "#F4A582",
    "marginal"     = "#D6604D"
  )

  ggplot2::ggplot(diag,
    ggplot2::aes(x = factor(time), y = facet_label, fill = method)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = method_colours, drop = FALSE) +
    ggplot2::labs(x = "Time", y = "", fill = "Method",
                  title = "Fitting Method by Model and Time") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
