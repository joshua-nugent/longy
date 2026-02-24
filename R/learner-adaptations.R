#' Automatic Learner Adaptations for Continuous Pseudo-Outcomes
#'
#' @description
#' When longy fits outcome models for TMLE or G-computation, the pseudo-outcomes
#' at intermediate time steps are continuous values in \eqn{[0,1]} (backward ICE
#' predictions). These are modeled with \code{quasibinomial} family to keep
#' predictions naturally bounded. However, many SuperLearner wrappers do not
#' support \code{quasibinomial}, so longy automatically adapts certain learners
#' at fit time. This page documents those adaptations.
#'
#' @section quasibinomial to binomial swap:
#' SuperLearner and most of its wrappers only recognize \code{"gaussian"} and
#' \code{"binomial"} as family strings. Since \code{quasibinomial} uses the same
#' link and variance function as \code{binomial}, longy substitutes
#' \code{binomial()} when calling SuperLearner. The quasibinomial family is
#' still used for the TMLE fluctuation step (which uses \code{glm.fit}
#' directly).
#'
#' @section SL.glmnet adaptation (h2o.glm):
#' \code{glmnet} does not support \code{quasibinomial} family. Fitting with
#' \code{gaussian} family produces predictions outside \eqn{[0,1]} that, when
#' clipped and logit-transformed in the TMLE fluctuation step, create extreme
#' values and inflated standard errors. Alternative workarounds (logit-transform
#' Y, fit gaussian, back-transform) also proved unstable in practice.
#'
#' When \code{SL.glmnet} is requested and the outcome is quasibinomial, longy
#' replaces it with an internal \code{SL.h2o.glm} wrapper that uses
#' \code{h2o::h2o.glm} with \code{family = "quasibinomial"}. Unlike
#' \code{glmnet}, \code{h2o.glm} natively supports quasibinomial with
#' elastic-net penalties (\code{alpha}, \code{lambda}), so predictions are
#' naturally bounded in \eqn{(0,1)} without transformation hacks.
#'
#' The wrapper:
#' \enumerate{
#'   \item Initializes h2o if not already running (\code{nthreads = -1, max_mem = "2g"})
#'   \item Converts data to h2o frames
#'   \item Fits \code{h2o.glm} with \code{family = "quasibinomial"},
#'     \code{alpha = 0.5}, and \code{lambda_search = TRUE}
#'   \item Extracts predicted probabilities (the \code{p1} column)
#'   \item Clips predictions to \eqn{[0.005, 0.995]}
#'   \item Cleans up temporary h2o frames
#' }
#'
#' If the \pkg{h2o} package is not installed, \code{SL.glmnet} is dropped from
#' the learner library with a warning.
#'
#' @section SL.xgboost adaptation:
#' \code{xgboost} with \code{binary:logistic} objective crashes when Y is
#' continuous in \eqn{[0,1]} rather than strictly binary. longy replaces
#' \code{SL.xgboost} with an internal wrapper \code{SL.xgboost.reg} that fits
#' with \code{gaussian} family (i.e., \code{reg:squarederror} objective) and
#' clips predictions to \eqn{[0.005, 0.995]}. Because tree-based methods
#' naturally restrict predictions to the range of the training data, the
#' transformation approach used for glmnet is not needed here.
#'
#' @section When do these adaptations apply?:
#' These swaps only occur when the internal family is \code{quasibinomial},
#' which happens in two contexts:
#' \enumerate{
#'   \item \strong{TMLE}: All backward ICE Q-model steps use quasibinomial,
#'     since pseudo-outcomes are continuous \eqn{[0,1]} predictions from the
#'     previous step (or scaled observed outcomes at the target time).
#'   \item \strong{G-computation}: Backward ICE steps for continuous outcomes
#'     also use quasibinomial after scaling Y to \eqn{[0,1]}.
#' }
#'
#' For binary outcomes at the target time (where Y is strictly 0/1), the
#' quasibinomial-to-binomial swap is harmless since the families are equivalent
#' on binary data.
#'
#' @section Other learners:
#' Learners not listed above (e.g., \code{SL.glm}, \code{SL.mean},
#' \code{SL.earth}, \code{SL.ranger}, \code{SL.gam}) are passed through
#' with the \code{binomial} family swap only. If a learner fails internally,
#' SuperLearner's built-in error handling drops it from the ensemble and
#' longy emits a warning reporting which learners failed.
#'
#' @section Interaction with cross-fitting:
#' The same adaptations apply within each cross-fitting fold. The wrapper
#' functions are injected into the SuperLearner evaluation environment at
#' call time, so they are available to both standard \code{SuperLearner} and
#' the parallel \code{ffSL} implementation.
#'
#' @name learner_adaptations
#' @aliases SL.h2o.glm SL.xgboost.reg
NULL
