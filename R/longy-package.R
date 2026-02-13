#' @keywords internal
#' @import data.table
"_PACKAGE"

# Suppress R CMD check NOTEs for non-standard evaluation symbols
utils::globalVariables(c(
  # data.table columns used in :=, .SD, etc.
  ".N", ".SD", ".time",
  ".longy_regime_consist", ".longy_cum_consist",
  ".longy_uncens", ".longy_cum_uncens",
  ".longy_consist_prev", ".longy_uncens_prev",
  ".longy_check_surv", ".longy_fold",
  ".sw_a", ".sw_c", ".sw_ac", ".csw_ac", ".sw_r", ".final_weight",
  ".marg_a", ".marg_c", ".marg_r",
  ".p_a", ".p_c", ".p_r",
  ".treatment", ".censored", ".observed",
  ".n_risk", ".method", ".id",
  # fit_outcome columns
  ".longy_Q", ".longy_regime_a", ".target_time",
  # .compute_cumulative_g columns
  ".g_a", ".g_c", ".g_point", ".g_cum",
  # Auto-detected observation column
  ".obs",
  # TMLE columns
  ".Q_star", ".Q_star_0", ".Q_star_next",
  ".H", ".aug", ".indicator", ".Y_T",
  ".tmp_order",
  # ggplot2 aes variables
  "time", "estimate", "ci_lower", "ci_upper", "regime", "estimator_label",
  # data.table special vars
  "nu", "nfold"
))

#' Simulated Longitudinal Dataset with Known Causal Effects
#'
#' A simulated dataset of 2000 subjects observed over 10 time points, with
#' time-varying confounders, informative censoring, intermittent outcome
#' measurement, and known true causal effects for validation.
#'
#' @format A data.frame with the following columns:
#' \describe{
#'   \item{id}{Subject identifier (integer 1-2000)}
#'   \item{time}{Time point (integer 0-9)}
#'   \item{W1}{Baseline covariate, continuous (N(0,1))}
#'   \item{W2}{Baseline covariate, binary (Bernoulli(0.3))}
#'   \item{L1}{Time-varying confounder, continuous (affected by past treatment)}
#'   \item{L2}{Time-varying confounder, binary}
#'   \item{A}{Treatment indicator (binary, 0/1)}
#'   \item{C}{Censoring indicator (binary, 1 = censored, absorbing)}
#'   \item{R}{Outcome observation indicator (binary, 1 = observed, intermittent)}
#'   \item{Y}{Outcome (binary, observed only when R=1 and C=0)}
#' }
#'
#' The dataset has an attribute `"true_effects"` containing a data.frame with
#' columns `time`, `EY1` (expected outcome under always-treat), and
#' `EY0` (expected outcome under never-treat) computed by
#' large-sample Monte Carlo simulation under the known DGP.
#'
#' @source Simulated via \code{data-raw/make_sim_longy.R}
"sim_longy"
