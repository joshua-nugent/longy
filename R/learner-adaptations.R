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
#' @section SL.glmnet adaptation:
#' \code{glmnet}'s \code{binomial} family rejects continuous \eqn{[0,1]} Y
#' (requires discrete 0/1 classes). longy replaces \code{SL.glmnet} with an
#' internal wrapper \code{SL.glmnet.reg} that fits with \code{gaussian} family
#' (penalized least squares on \eqn{[0,1]} pseudo-outcomes) and clips
#' predictions to \eqn{[0.005, 0.995]}. A message is printed to the console
#' when this swap occurs.
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
#' @section SL.ranger adaptation:
#' \code{ranger}'s formula interface rejects column names containing special
#' characters (e.g., dummy columns like \code{uacr_c_>=300}). longy
#' \strong{always} replaces \code{SL.ranger} with an internal wrapper
#' \code{SL.ranger.longy} that uses ranger's non-formula interface
#' (\code{x = X, y = Y}). This avoids formula-related errors regardless of
#' column names. When the outcome is continuous \eqn{[0,1]} (quasibinomial
#' context), predictions are clipped to \eqn{[0.005, 0.995]} as a safety
#' measure.
#'
#' @section When do these adaptations apply?:
#' The \code{SL.ranger} adaptation applies unconditionally (all model fits).
#' The \code{SL.xgboost} and \code{SL.glmnet} adaptations apply only when the
#' internal family is \code{quasibinomial}, which happens in two contexts:
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
#' \code{SL.earth}, \code{SL.gam}, \code{SL.nnet}) are passed through
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
#' @aliases SL.xgboost.reg SL.glmnet.reg SL.ranger.longy
NULL
