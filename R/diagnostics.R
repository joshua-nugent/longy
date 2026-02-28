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

#' Prediction Diagnostics
#'
#' Summarizes the distribution of predicted probabilities from each fitted
#' nuisance model (treatment, censoring, observation) by time point. Useful
#' for spotting near-positivity violations and understanding which models are
#' driving extreme weights.
#'
#' @param obj A \code{longy_data} object with at least one nuisance model fit.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available regime.
#' @param model Character. Which model(s) to include: \code{"all"} (default),
#'   \code{"treatment"}, \code{"censoring"}, or \code{"observation"}.
#'
#' @return A data.table with columns: \code{model}, \code{time}, \code{n_risk},
#'   \code{marginal}, \code{min}, \code{p05}, \code{median}, \code{mean},
#'   \code{p95}, \code{max}.
#' @export
prediction_diagnostics <- function(obj, regime = NULL, model = "all") {
  obj <- .as_longy_data(obj)

  model <- match.arg(model, c("all", "treatment", "censoring", "observation"))
  include <- if (model == "all") {
    c("treatment", "censoring", "observation")
  } else {
    model
  }

  # Resolve regime
  if (is.null(regime)) {
    available <- names(obj$fits$treatment)
    if (length(available) == 0) available <- names(obj$fits$censoring)
    if (length(available) == 0) available <- names(obj$fits$observation)
    regime <- available[1]
  }
  if (is.na(regime) || is.null(regime)) {
    stop("No model fits found. Run fit_treatment/fit_censoring/fit_observation first.",
         call. = FALSE)
  }

  rows <- list()

  # Treatment
  trt_fit <- obj$fits$treatment[[regime]]
  if ("treatment" %in% include && !is.null(trt_fit)) {
    preds <- trt_fit$predictions
    rows[[length(rows) + 1L]] <- preds[, list(
      model = "treatment",
      n_risk = .N,
      marginal = round(.marg_a[1], 4),
      min = round(min(.p_a), 4),
      p05 = round(stats::quantile(.p_a, 0.05), 4),
      median = round(stats::median(.p_a), 4),
      mean = round(mean(.p_a), 4),
      p95 = round(stats::quantile(.p_a, 0.95), 4),
      max = round(max(.p_a), 4)
    ), by = .time]
  }

  # Censoring (per cause)
  cens_fits <- obj$fits$censoring[[regime]]
  if ("censoring" %in% include && length(cens_fits) > 0) {
    for (cvar in names(cens_fits)) {
      preds <- cens_fits[[cvar]]$predictions
      label <- paste0("censoring: ", sub("^\\.cens_", "", cvar))
      rows[[length(rows) + 1L]] <- preds[, list(
        model = label,
        n_risk = .N,
        marginal = round(.marg_c[1], 4),
        min = round(min(.p_c), 4),
        p05 = round(stats::quantile(.p_c, 0.05), 4),
        median = round(stats::median(.p_c), 4),
        mean = round(mean(.p_c), 4),
        p95 = round(stats::quantile(.p_c, 0.95), 4),
        max = round(max(.p_c), 4)
      ), by = .time]
    }
  }

  # Observation
  obs_fit <- obj$fits$observation[[regime]]
  if ("observation" %in% include && !is.null(obs_fit)) {
    preds <- obs_fit$predictions
    rows[[length(rows) + 1L]] <- preds[, list(
      model = "observation",
      n_risk = .N,
      marginal = round(.marg_r[1], 4),
      min = round(min(.p_r), 4),
      p05 = round(stats::quantile(.p_r, 0.05), 4),
      median = round(stats::median(.p_r), 4),
      mean = round(mean(.p_r), 4),
      p95 = round(stats::quantile(.p_r, 0.95), 4),
      max = round(max(.p_r), 4)
    ), by = .time]
  }

  if (length(rows) == 0) {
    message("No model fits found for the requested model type(s).")
    return(invisible(data.table::data.table(
      model = character(), time = integer(), n_risk = integer(),
      marginal = numeric(), min = numeric(), p05 = numeric(),
      median = numeric(), mean = numeric(), p95 = numeric(),
      max = numeric()
    )))
  }

  result <- data.table::rbindlist(rows)
  data.table::setnames(result, ".time", "time")
  data.table::setkey(result, model, time)
  result
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
        sl_risk = list(entry$sl_risk %||% NA_real_),
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
          sl_risk = list(entry$sl_risk %||% NA_real_),
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
        sl_risk = list(entry$sl_risk %||% NA_real_),
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
        sl_risk = list(entry$sl_risk %||% NA_real_),
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

  # Accept longy_data, longy_result, or pre-computed data.table
  if (is.data.frame(obj) && !.is_longy_data_structure(obj)) {
    diag <- data.table::as.data.table(data.table::copy(obj))
  } else {
    obj <- tryCatch(.as_longy_data(obj), error = function(e) NULL)
    if (is.null(obj)) {
      stop("obj must be a longy_data object or sl_diagnostics() output.",
           call. = FALSE)
    }
    diag <- sl_diagnostics(obj, regime = regime, model = model)
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
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
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
    ggplot2::aes(x = time, y = cv_risk, colour = learner, group = learner)) +
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


# ============================================================================
# Influence / Estimate Decomposition Diagnostics
# ============================================================================

#' Influence Diagnostics: Decompose Adjusted vs Unadjusted Estimates
#'
#' Creates a per-subject-time data.table that merges outcomes, IPW weights,
#' G-comp counterfactual predictions, and propensity scores. Useful for
#' understanding why adjusted estimates differ from unadjusted means.
#'
#' For each time point, subjects are grouped by weight quintile. The summary
#' shows how the outcome distribution varies across weight groups, revealing
#' whether high-weight subjects have systematically different outcomes.
#'
#' @param obj A \code{longy_data} object with at least one estimator run.
#' @param regime Character. Regime name. If NULL, uses the first available.
#' @param times Numeric vector. Time points to include. If NULL, all available.
#' @param n_quantiles Integer. Number of weight groups for the summary
#'   (default 5 = quintiles).
#'
#' @return A list of class \code{"longy_influence_diag"} with elements:
#'   \describe{
#'     \item{subject_dt}{Per-subject-time data.table with columns: id, time,
#'       outcome, weight, Q_hat, p_a, marg_a, weight_group}
#'     \item{summary_dt}{Summary by weight group and time: n, mean_outcome,
#'       mean_weight, mean_p_a, weighted_contribution, mean_Q_hat}
#'     \item{time_summary}{Per-time-point comparison of unadjusted mean,
#'       IPW estimate, G-comp estimate, and weight-outcome correlation}
#'     \item{regime}{Regime name}
#'   }
#'
#' @export
influence_diagnostics <- function(obj, regime = NULL, times = NULL,
                                   n_quantiles = 5L) {
  obj <- .as_longy_data(obj)

  # Resolve regime
  if (is.null(regime)) {
    available <- names(obj$weights)
    if (length(available) == 0) available <- names(obj$fits$treatment)
    if (length(available) == 0) available <- names(obj$fits$outcome)
    regime <- available[1]
  }
  if (is.null(regime) || is.na(regime)) {
    stop("No fitted regime found. Run fit_treatment() or estimate_*() first.",
         call. = FALSE)
  }

  nodes <- obj$nodes
  id_col <- nodes$id
  time_col <- nodes$time
  y_col <- nodes$outcome
  dt <- obj$data

  # --- Determine available components ---
  has_weights <- !is.null(obj$weights[[regime]]) &&
    length(obj$weights[[regime]]) > 0
  has_treatment <- !is.null(obj$fits$treatment[[regime]])
  has_outcome <- !is.null(obj$fits$outcome[[regime]]) &&
    !is.null(obj$fits$outcome[[regime]]$predictions)

  if (!has_weights && !has_outcome) {
    stop("Need at least weights (IPW) or outcome model (G-comp) fitted.",
         call. = FALSE)
  }

  # --- Build per-subject-time data ---
  # Start with the weight data (has regime followers at each time)
  if (has_weights) {
    w_dt <- obj$weights[[regime]]$weights_dt
    base <- w_dt[, c(id_col, ".time", ".sw_a", ".sw_c", ".csw_ac",
                      ".sw_r", ".final_weight"), with = FALSE]
  } else {
    # If no weights, build from all subjects at each time
    all_times <- obj$meta$time_values
    base <- dt[, c(id_col, time_col), with = FALSE]
    data.table::setnames(base, time_col, ".time")
    base[, .final_weight := 1]
  }

  # Resolve time points
  available_times <- sort(unique(base$.time))
  if (is.null(times)) {
    times <- available_times
  } else {
    times <- intersect(times, available_times)
  }
  base <- base[base$.time %in% times, ]

  # Merge outcome from raw data
  y_dt <- dt[, c(id_col, time_col, y_col), with = FALSE]
  data.table::setnames(y_dt, c(time_col, y_col), c(".time", ".outcome"))
  base <- merge(base, y_dt, by = c(id_col, ".time"), all.x = TRUE)

  # Merge treatment propensity
  if (has_treatment) {
    trt_preds <- obj$fits$treatment[[regime]]$predictions
    p_dt <- trt_preds[, c(id_col, ".time", ".p_a", ".marg_a"), with = FALSE]
    base <- merge(base, p_dt, by = c(id_col, ".time"), all.x = TRUE)
  }

  # Merge G-comp predictions (Q_hat from earliest backward step)
  if (has_outcome) {
    q_preds <- obj$fits$outcome[[regime]]$predictions
    # Q_hat is stored per target_time; id column uses original id name
    q_dt <- q_preds[, c(id_col, ".target_time", ".Q_hat"), with = FALSE]
    data.table::setnames(q_dt, ".target_time", ".time")
    base <- merge(base, q_dt, by = c(id_col, ".time"), all.x = TRUE)
  }

  # --- Compute influence metrics ---
  # Weight quintiles (within each time point)
  if (has_weights) {
    base[, .weight_group := {
      if (all(is.na(.final_weight))) {
        rep(NA_integer_, .N)
      } else {
        brks <- stats::quantile(.final_weight, probs = seq(0, 1, length.out = n_quantiles + 1L),
                                 na.rm = TRUE)
        # Handle ties in quantile boundaries
        brks <- unique(brks)
        if (length(brks) < 2) {
          rep(1L, .N)
        } else {
          as.integer(cut(.final_weight, breaks = brks, include.lowest = TRUE,
                          labels = FALSE))
        }
      }
    }, by = .time]
  } else {
    base[, .weight_group := 1L]
  }

  # --- Compute g-comp estimates from FULL predictions (all subjects) ---
  # G-comp uses mean(Q_hat) over ALL subjects, not just regime followers.
  # The base data only has followers (from weights), so we compute separately.
  gcomp_by_time <- NULL
  if (has_outcome) {
    q_preds_full <- obj$fits$outcome[[regime]]$predictions
    gcomp_by_time <- q_preds_full[, list(
      .gcomp_est = mean(.Q_hat, na.rm = TRUE),
      .gcomp_n = .N
    ), by = .target_time]
    data.table::setnames(gcomp_by_time, ".target_time", ".time")
  }

  # --- Compute population unadjusted mean (all subjects, not just followers) ---
  pop_unadj <- dt[, list(
    .pop_unadj = mean(as.numeric(get(y_col)), na.rm = TRUE),
    .pop_n = sum(!is.na(get(y_col)))
  ), by = c(time_col)]
  data.table::setnames(pop_unadj, time_col, ".time")

  # --- Time-level summary ---
  time_summary <- base[!is.na(.outcome), {
    follower_mean <- mean(.outcome, na.rm = TRUE)
    wm <- if (has_weights && any(!is.na(.final_weight))) {
      stats::weighted.mean(.outcome, .final_weight, na.rm = TRUE)
    } else {
      NA_real_
    }
    wt_y_cor <- if (has_weights && sum(!is.na(.final_weight) & !is.na(.outcome)) > 2) {
      tryCatch(stats::cor(.final_weight, .outcome, use = "complete.obs"),
               error = function(e) NA_real_)
    } else {
      NA_real_
    }
    list(
      n_followers = .N,
      follower_mean = round(follower_mean, 2),
      ipw_estimate = round(wm, 2),
      weight_outcome_cor = round(wt_y_cor, 3),
      mean_weight = round(mean(.final_weight, na.rm = TRUE), 3),
      max_weight = round(max(.final_weight, na.rm = TRUE), 2),
      ess = if (has_weights) round(.ess(.final_weight), 1) else as.numeric(.N)
    )
  }, by = .time]

  # Merge population unadjusted mean
  time_summary <- merge(time_summary, pop_unadj, by = ".time", all.x = TRUE)
  time_summary[, pop_unadj_mean := round(.pop_unadj, 2)]
  time_summary[, .pop_unadj := NULL]

  # Merge g-comp estimates from full predictions
  if (!is.null(gcomp_by_time)) {
    time_summary <- merge(time_summary, gcomp_by_time, by = ".time", all.x = TRUE)
    time_summary[, gcomp_estimate := round(.gcomp_est, 2)]
    time_summary[, c(".gcomp_est", ".gcomp_n") := NULL]
  } else {
    time_summary[, gcomp_estimate := NA_real_]
  }

  data.table::setkey(time_summary, .time)

  # --- Weight-group summary ---
  group_summary <- base[!is.na(.outcome) & !is.na(.weight_group), {
    out <- list(
      n = .N,
      mean_outcome = round(mean(.outcome, na.rm = TRUE), 2),
      mean_weight = round(mean(.final_weight, na.rm = TRUE), 3)
    )
    if (has_treatment && ".p_a" %in% names(base)) {
      out$mean_p_a <- round(mean(.p_a, na.rm = TRUE), 4)
    }
    if (has_outcome && ".Q_hat" %in% names(base)) {
      out$mean_Q_hat <- round(mean(.Q_hat, na.rm = TRUE), 2)
    }
    # Contribution: what fraction of the total weighted sum comes from this group
    total_wsum <- sum(.final_weight * .outcome, na.rm = TRUE)
    out$weighted_sum <- round(total_wsum, 1)
    out
  }, by = list(.time, .weight_group)]
  data.table::setkey(group_summary, .time, .weight_group)

  # Add contribution percentage (within each time point)
  group_summary[, contribution_pct := {
    total <- sum(weighted_sum, na.rm = TRUE)
    if (total != 0) round(100 * weighted_sum / total, 1) else NA_real_
  }, by = .time]

  result <- list(
    subject_dt = base,
    summary_dt = group_summary,
    time_summary = time_summary,
    regime = regime,
    n_quantiles = n_quantiles
  )
  class(result) <- "longy_influence_diag"
  result
}


#' @export
print.longy_influence_diag <- function(x, ...) {
  cat(sprintf("Influence Diagnostics -- regime: %s\n", x$regime))
  cat(strrep("=", 55), "\n\n")

  ts <- x$time_summary
  for (i in seq_len(nrow(ts))) {
    row <- ts[i]
    cat(sprintf("Time %s (n_pop=%d, n_followers=%d, ESS=%.0f):\n",
                as.character(row$.time), row$.pop_n, row$n_followers, row$ess))
    cat(sprintf("  Population mean (all):   %8.2f\n", row$pop_unadj_mean))
    cat(sprintf("  Follower mean (unwt'd):  %8.2f\n", row$follower_mean))
    if (!is.na(row$ipw_estimate))
      cat(sprintf("  IPW estimate (wt'd):     %8.2f\n", row$ipw_estimate))
    if (!is.na(row$gcomp_estimate))
      cat(sprintf("  G-comp estimate:         %8.2f\n", row$gcomp_estimate))
    if (!is.na(row$weight_outcome_cor))
      cat(sprintf("  Weight-outcome corr:     %8.3f  %s\n",
                  row$weight_outcome_cor,
                  if (row$weight_outcome_cor < -0.1) "(high-weight obs have LOWER outcomes)"
                  else if (row$weight_outcome_cor > 0.1) "(high-weight obs have HIGHER outcomes)"
                  else "(weak association)"))
    cat(sprintf("  Mean weight: %.3f | Max weight: %.2f\n\n",
                row$mean_weight, row$max_weight))

    # Weight group breakdown
    gs <- x$summary_dt[.time == row$.time]
    if (nrow(gs) > 0) {
      cat("  By weight group:\n")
      header <- sprintf("  %5s  %5s  %10s  %10s",
                        "group", "n", "mean_Y", "mean_wt")
      if ("mean_p_a" %in% names(gs)) header <- paste0(header, sprintf("  %8s", "mean_pA"))
      if ("mean_Q_hat" %in% names(gs)) header <- paste0(header, sprintf("  %10s", "mean_Qhat"))
      header <- paste0(header, sprintf("  %8s", "contrib%"))
      cat(header, "\n")
      for (j in seq_len(nrow(gs))) {
        g <- gs[j]
        line <- sprintf("  %5s  %5d  %10.2f  %10.3f",
                        paste0("Q", g$.weight_group), g$n,
                        g$mean_outcome, g$mean_weight)
        if ("mean_p_a" %in% names(gs))
          line <- paste0(line, sprintf("  %8.4f", g$mean_p_a))
        if ("mean_Q_hat" %in% names(gs))
          line <- paste0(line, sprintf("  %10.2f", g$mean_Q_hat))
        line <- paste0(line, sprintf("  %7.1f%%", g$contribution_pct))
        cat(line, "\n")
      }
      cat("\n")
    }
  }
  invisible(x)
}


#' Plot Influence Diagnostics
#'
#' Multi-panel plot showing the relationship between weights, outcomes, and
#' counterfactual predictions. Helps diagnose why adjusted estimates differ
#' from unadjusted means.
#'
#' @param x A \code{longy_influence_diag} object from
#'   \code{influence_diagnostics()}.
#' @param type Character. Plot type:
#'   \describe{
#'     \item{\code{"scatter"}}{Weight vs outcome scatter with loess smooth,
#'       faceted by time. Shows the correlation driving the difference.}
#'     \item{\code{"quintile"}}{Bar chart of mean outcome by weight quintile.
#'       Shows subgroup patterns driving the IPW estimate.}
#'     \item{\code{"density"}}{Overlaid densities of observed Y vs G-comp
#'       counterfactual predictions Q_hat. Shows distributional shift.}
#'     \item{\code{"cumulative"}}{Cumulative contribution plot. Subjects sorted
#'       by absolute influence; shows how concentrated the estimate is.}
#'   }
#' @param ... Unused.
#'
#' @return A ggplot2 object.
#' @export
plot_influence_diagnostics <- function(x, type = c("scatter", "quintile",
                                                     "density", "cumulative"),
                                        ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)
  if (!inherits(x, "longy_influence_diag"))
    stop("x must be a longy_influence_diag object.", call. = FALSE)

  type <- match.arg(type)

  switch(type,
    scatter    = .plot_infl_scatter(x),
    quintile   = .plot_infl_quintile(x),
    density    = .plot_infl_density(x),
    cumulative = .plot_infl_cumulative(x)
  )
}


# --- Internal plot helpers --------------------------------------------------

#' Scatter: weight vs outcome
#' @noRd
.plot_infl_scatter <- function(x) {
  dt <- x$subject_dt[!is.na(.outcome) & !is.na(.final_weight)]
  if (nrow(dt) == 0) {
    message("No data to plot.")
    return(invisible(NULL))
  }

  dt[, time_label := paste("t =", .time)]

  ggplot2::ggplot(dt, ggplot2::aes(x = .final_weight, y = .outcome)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.8, colour = "#2166AC") +
    ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "#D6604D",
                          linewidth = 1) +
    ggplot2::facet_wrap(~time_label, scales = "free") +
    ggplot2::labs(
      x = "IPW Weight", y = "Outcome",
      title = sprintf("Weight vs Outcome -- %s", x$regime),
      subtitle = "Negative slope = high-weight subjects have lower outcomes, pulling IPW estimate down"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}


#' Bar chart: mean outcome by weight quintile
#' @noRd
.plot_infl_quintile <- function(x) {
  gs <- x$summary_dt
  if (nrow(gs) == 0) {
    message("No summary data to plot.")
    return(invisible(NULL))
  }

  gs[, time_label := paste("t =", .time)]
  gs[, group_label := paste0("Q", .weight_group)]

  # Two-bar comparison: mean outcome (unweighted) and contribution
  ggplot2::ggplot(gs, ggplot2::aes(x = group_label, y = mean_outcome)) +
    ggplot2::geom_col(fill = "#2166AC", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("w=%.1f", mean_weight)),
                        vjust = -0.3, size = 3, colour = "#666666") +
    ggplot2::facet_wrap(~time_label, scales = "free_y") +
    ggplot2::labs(
      x = "Weight Quintile (Q1=lowest, Q5=highest)",
      y = "Mean Outcome",
      title = sprintf("Mean Outcome by Weight Group -- %s", x$regime),
      subtitle = "Labels show mean weight per group. Declining pattern = confounding drives IPW below unadjusted."
    ) +
    ggplot2::theme_minimal(base_size = 13)
}


#' Density: observed Y vs Q_hat
#' @noRd
.plot_infl_density <- function(x) {
  dt <- x$subject_dt

  if (!".Q_hat" %in% names(dt) || all(is.na(dt$.Q_hat))) {
    message("No G-comp predictions (Q_hat) available. Run fit_outcome() first.")
    return(invisible(NULL))
  }

  # Stack observed and predicted for overlaid density
  obs <- dt[!is.na(.outcome), list(value = .outcome, source = "Observed Y",
                                     .time = .time)]
  pred <- dt[!is.na(.Q_hat), list(value = .Q_hat, source = "Counterfactual Q_hat",
                                    .time = .time)]
  stacked <- data.table::rbindlist(list(obs, pred))
  stacked[, time_label := paste("t =", .time)]

  ggplot2::ggplot(stacked, ggplot2::aes(x = value, fill = source, colour = source)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 0.7) +
    ggplot2::facet_wrap(~time_label, scales = "free") +
    ggplot2::scale_fill_manual(values = c("Observed Y" = "#2166AC",
                                           "Counterfactual Q_hat" = "#D6604D")) +
    ggplot2::scale_colour_manual(values = c("Observed Y" = "#2166AC",
                                              "Counterfactual Q_hat" = "#D6604D")) +
    ggplot2::labs(
      x = "Value", y = "Density",
      fill = "", colour = "",
      title = sprintf("Observed Outcome vs Counterfactual Prediction -- %s",
                      x$regime),
      subtitle = "Shift between distributions shows the confounding adjustment from the outcome model"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "bottom")
}


#' Cumulative contribution (Lorenz-like)
#' @noRd
.plot_infl_cumulative <- function(x) {
  dt <- x$subject_dt[!is.na(.outcome) & !is.na(.final_weight)]
  if (nrow(dt) == 0) {
    message("No data to plot.")
    return(invisible(NULL))
  }

  dt[, time_label := paste("t =", .time)]

  # Per-subject influence: weighted contribution relative to uniform
  cum_dt <- dt[, {
    total_w <- sum(.final_weight)
    contrib <- (.final_weight / total_w) * .outcome
    # Sort by absolute deviation from uniform contribution
    uniform <- .outcome / .N
    influence <- abs(contrib - uniform)
    ord <- order(influence, decreasing = TRUE)
    cum_pct <- cumsum(rep(1 / .N, .N))
    cum_estimate <- cumsum(contrib[ord])
    # Normalize cumulative estimate to final value
    final_val <- sum(contrib)
    list(
      pct_subjects = cum_pct * 100,
      cum_weighted_mean = cum_estimate / final_val * 100
    )
  }, by = time_label]

  ggplot2::ggplot(cum_dt, ggplot2::aes(x = pct_subjects,
                                         y = cum_weighted_mean)) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          colour = "grey50") +
    ggplot2::facet_wrap(~time_label) +
    ggplot2::labs(
      x = "% of Subjects (sorted by influence)",
      y = "% of IPW Estimate",
      title = sprintf("Cumulative Influence -- %s", x$regime),
      subtitle = "Deviation from diagonal = concentration of influence in few subjects"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

#' TMLE Diagnostics
#'
#' Summarizes TMLE fluctuation steps and EIF decomposition for each target
#' time point. Shows epsilon values, Q-model methods, risk set sizes, and
#' the relative contribution of the initial (outcome model) vs augmentation
#' (IPW correction) terms to the EIF variance.
#'
#' @param obj A \code{longy_data} object with TMLE results, or a
#'   \code{longy_result} from \code{estimate_tmle()}.
#' @param regime Character. Name of the regime. If NULL, uses the first
#'   available TMLE result.
#'
#' @return Invisibly returns a list with \code{steps} (data.table of per-step
#'   info) and \code{eif} (data.table of EIF decomposition). Also prints a
#'   formatted summary.
#' @export
tmle_diagnostics <- function(obj, regime = NULL) {
  obj <- .as_longy_data(obj)

  # Find TMLE result
  tmle_keys <- grep("_tmle$", names(obj$results), value = TRUE)
  if (length(tmle_keys) == 0) {
    stop("No TMLE results found. Run estimate_tmle() first.", call. = FALSE)
  }

  if (!is.null(regime)) {
    key <- paste0(regime, "_tmle")
    if (!key %in% tmle_keys) {
      stop(sprintf("No TMLE result for regime '%s'. Available: %s",
                   regime, paste(sub("_tmle$", "", tmle_keys), collapse = ", ")),
           call. = FALSE)
    }
  } else {
    key <- tmle_keys[1]
    regime <- sub("_tmle$", "", key)
  }

  res <- obj$results[[key]]
  tmle_info <- res$tmle_info

  if (is.null(tmle_info) || length(tmle_info) == 0) {
    cat("No TMLE diagnostic info stored. Re-run estimate_tmle() to generate diagnostics.\n")
    return(invisible(list()))
  }

  steps <- tmle_info$steps
  eif <- tmle_info$eif

  cat(sprintf("TMLE Diagnostics -- regime: %s\n", regime))
  cat(paste(rep("=", 60), collapse = ""), "\n")

  if (!is.null(steps) && nrow(steps) > 0) {
    target_times <- sort(unique(steps$target_time))

    for (tt in target_times) {
      tt_steps <- steps[steps$target_time == tt, ]
      tt_steps <- tt_steps[order(tt_steps$step_time, decreasing = TRUE), ]

      cat(sprintf("\nTarget time %d (%d backward steps):\n", tt, nrow(tt_steps)))

      # Summary line: total epsilon, methods
      total_eps <- sum(abs(tt_steps$epsilon))
      methods_used <- unique(tt_steps$method[tt_steps$method != "none"])
      cat(sprintf("  Sum |epsilon|: %.4f | Methods: %s\n",
                  total_eps, paste(methods_used, collapse = ", ")))

      # Per-step table
      cat(sprintf("  %6s %9s %8s %7s %7s %8s %8s %8s\n",
                  "step_t", "epsilon", "method", "n_risk", "n_fluc",
                  "mean_Qb", "mean_Q*", "min_g"))
      cat(sprintf("  %s\n", paste(rep("-", 72), collapse = "")))

      for (j in seq_len(nrow(tt_steps))) {
        s <- tt_steps[j, ]
        method_short <- switch(s$method,
                               SuperLearner = "SL",
                               "cross-fitted" = "CF",
                               marginal = "marg",
                               none = "none",
                               glm = "glm",
                               s$method)
        cat(sprintf("  %6d %9.4f %8s %7d %7d %8.4f %8.4f %8.4f\n",
                    s$step_time, s$epsilon, method_short,
                    s$n_at_risk, s$n_fluct,
                    ifelse(is.na(s$mean_Q_bar), 0, s$mean_Q_bar),
                    ifelse(is.na(s$mean_Q_star), 0, s$mean_Q_star),
                    ifelse(is.na(s$min_g_denom), 0, s$min_g_denom)))
      }

      # Flag large epsilons
      large_eps <- tt_steps[abs(tt_steps$epsilon) > 0.5, ]
      if (nrow(large_eps) > 0) {
        cat(sprintf("  WARNING: %d step(s) with |epsilon| > 0.5 (large fluctuation)\n",
                    nrow(large_eps)))
      }

      # Flag small fluctuation sets
      small_fluct <- tt_steps[tt_steps$n_fluct > 0 & tt_steps$n_fluct < 10, ]
      if (nrow(small_fluct) > 0) {
        cat(sprintf("  WARNING: %d step(s) with <10 subjects in fluctuation set\n",
                    nrow(small_fluct)))
      }

      # Flag near-positivity violations
      if (any(!is.na(tt_steps$min_g_denom) & tt_steps$min_g_denom < 0.02)) {
        cat("  WARNING: min g_denom < 0.02 at some steps (near-positivity violation)\n")
      }
    }
  }

  # EIF decomposition
  if (!is.null(eif) && nrow(eif) > 0) {
    cat(sprintf("\nEIF Variance Decomposition:\n"))
    cat(sprintf("  %6s %8s %8s %8s %7s %8s\n",
                "time", "SE(D)", "SE(init)", "SE(aug)", "pct_aug", "n_aug>0"))
    cat(sprintf("  %s\n", paste(rep("-", 54), collapse = "")))

    for (j in seq_len(nrow(eif))) {
      e <- eif[j, ]
      # Variance decomposition: what fraction of total variance comes from augmentation?
      var_total <- e$se_total^2
      var_aug <- e$se_augmentation^2
      pct_aug <- if (!is.na(var_total) && var_total > 0) {
        100 * var_aug / var_total
      } else {
        NA_real_
      }
      cat(sprintf("  %6d %8.4f %8.4f %8.4f %6.1f%% %5d/%d\n",
                  e$target_time, e$se_total, e$se_initial, e$se_augmentation,
                  ifelse(is.na(pct_aug), 0, pct_aug),
                  e$n_aug_nonzero, e$n_subjects))
    }

    cat("\n  SE(init) = SE of Q*_0 - psi (outcome model contribution)\n")
    cat("  SE(aug)  = SE of sum H_s*(Q*_{s+1} - Q*_s) (IPW correction)\n")
    cat("  High pct_aug = SE driven by weight variability, not outcome model\n")
  }

  cat("\n")
  invisible(list(steps = steps, eif = eif, regime = regime))
}
