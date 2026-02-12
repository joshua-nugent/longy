#' @keywords internal
"_PACKAGE"

#' Simulated Longitudinal Dataset with Known Causal Effects
#'
#' A simulated dataset of 500 subjects observed over 10 time points, with
#' time-varying confounders, informative censoring, intermittent outcome
#' measurement, and known true causal effects for validation.
#'
#' @format A data.frame with the following columns:
#' \describe{
#'   \item{id}{Subject identifier (integer 1-500)}
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
#' columns `time`, `EY1` (E[Y_t(a=1)]), and `EY0` (E[Y_t(a=0)]) computed by
#' large-sample Monte Carlo simulation under the known DGP.
#'
#' @source Simulated via \code{data-raw/make_sim_longy.R}
"sim_longy"
