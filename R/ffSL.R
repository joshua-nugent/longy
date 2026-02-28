# Internal future-factorial SuperLearner
#
# Adapted from ffSL by Joshua Nugent (originally based on mcSuperLearner by
# Eric Polley). Distributes CV fold x algorithm combinations via
# future.apply::future_lapply for efficient parallelism on multi-core / HPC.
#
# NOT exported. Called by .safe_sl() when sl_fn = "ffSL".

#' Future Factorial SuperLearner (internal)
#'
#' Same interface and return type as \code{SuperLearner::SuperLearner()}.
#' Parallelizes across all fold x algorithm combinations using
#' \code{future.apply::future_lapply}.
#'
#' @inheritParams SuperLearner::SuperLearner
#' @return An object of class \code{"SuperLearner"}.
#' @noRd
.ffSL <- function(Y, X, newX = NULL, family = stats::gaussian(), SL.library,
                  method = 'method.NNLS', id = NULL, verbose = FALSE,
                  control = list(), cvControl = list(), obsWeights = NULL,
                  env = parent.frame()) {

  if (verbose) message("Running Future Factorial Super Learner")

  time_start <- proc.time()

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for ffSL. ",
         "Install with install.packages('future.apply').", call. = FALSE)
  }

  if (is.character(method)) {
    if (exists(method, mode = 'list')) {
      method <- get(method, mode = 'list')
    } else if (exists(method, mode = 'function')) {
      method <- get(method, mode = 'function')()
    } else if (exists(method, envir = asNamespace("SuperLearner"))) {
      method <- get(method, envir = asNamespace("SuperLearner"))
      if (is.function(method)) method <- method()
    }
  } else if (is.function(method)) {
    method <- method()
  }
  if (!is.list(method)) {
    stop("method is not in the appropriate format. Check out help('method.template')")
  }
  if (!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }

  control <- do.call(SuperLearner::SuperLearner.control, control)
  cvControl <- do.call(SuperLearner::SuperLearner.CV.control, cvControl)

  library <- utils::getFromNamespace(".createLibrary", "SuperLearner")(SL.library)
  utils::getFromNamespace(".check.SL.library", "SuperLearner")(
    library = c(unique(library$library$predAlgorithm), library$screenAlgorithm)
  )

  call <- match.call(expand.dots = TRUE)

  if (!inherits(X, 'data.frame') && verbose) {
    message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
  }

  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
  k <- nrow(library$library)
  kScreen <- length(library$screenAlgorithm)
  Z <- matrix(NA, N, k)
  libraryNames <- paste(library$library$predAlgorithm,
                         library$screenAlgorithm[library$library$rowScreen],
                         sep = "_")

  if (p < 2 & !identical(library$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }

  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k), envir = fitLibEnv)
  assign('libraryNames', libraryNames, envir = fitLibEnv)
  evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)

  errorsInCVLibrary <- rep(0, k)
  errorsInLibrary <- rep(0, k)

  if (is.null(newX)) {
    newX <- X
  }
  if (!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) {
    stop("the outcome Y must be a numeric vector")
  }

  if (is.character(family))
    family <- get(family, mode = "function", envir = env)
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (family$family != "binomial" & isTRUE("cvAUC" %in% method$require)) {
    stop("'method.AUC' is designed for the 'binomial' family only")
  }

  validRows <- utils::getFromNamespace("CVFolds", "SuperLearner")(
    N = N, id = id, Y = Y, cvControl = cvControl
  )

  if (is.null(id)) {
    id <- seq(N)
  }
  if (!identical(length(id), N)) {
    stop("id vector must have the same dimension as Y")
  }

  if (is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if (!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }

  # Single CV fold x algorithm task
  .singleCVTask <- function(task, Y, dataX, id, obsWeights, library,
                             kScreen, p, libraryNames, validRows,
                             family, verbose, env) {
    fold_idx <- task$fold_idx
    alg_idx <- task$alg_idx
    valid <- validRows[[fold_idx]]

    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- Y[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]

    screenIndex <- library$library$rowScreen[alg_idx]
    screen_fn_name <- library$screenAlgorithm[screenIndex]

    if (exists(screen_fn_name, envir = env)) {
      screen_fn <- get(screen_fn_name, envir = env)
    } else if (exists(screen_fn_name, envir = asNamespace("SuperLearner"))) {
      screen_fn <- get(screen_fn_name, envir = asNamespace("SuperLearner"))
    } else {
      stop(paste("Cannot find screening function:", screen_fn_name))
    }

    whichScreen <- tryCatch(
      do.call(screen_fn, list(
        Y = tempOutcome, X = tempLearn, family = family,
        id = tempId, obsWeights = tempObsWeights
      )),
      error = function(e) NULL
    )
    if (is.null(whichScreen)) {
      whichScreen <- rep(TRUE, p)
    }

    pred_fn_name <- library$library$predAlgorithm[alg_idx]

    if (exists(pred_fn_name, envir = env)) {
      pred_fn <- get(pred_fn_name, envir = env)
    } else if (exists(pred_fn_name, envir = asNamespace("SuperLearner"))) {
      pred_fn <- get(pred_fn_name, envir = asNamespace("SuperLearner"))
    } else {
      stop(paste("Cannot find function:", pred_fn_name))
    }

    testAlg <- tryCatch(
      do.call(pred_fn, list(
        Y = tempOutcome,
        X = tempLearn[, whichScreen, drop = FALSE],
        newX = tempValid[, whichScreen, drop = FALSE],
        family = family, id = tempId, obsWeights = tempObsWeights
      )),
      error = function(e) {
        message(sprintf("ffSL: %s failed in CV fold %d: %s",
                        libraryNames[alg_idx], fold_idx, e$message))
        NULL
      }
    )

    if (is.null(testAlg)) {
      pred <- rep(NA, nrow(tempValid))
    } else {
      pred <- testAlg$pred
    }

    if (verbose) message(paste("CV fold", fold_idx, "-", libraryNames[alg_idx]))

    list(fold_idx = fold_idx, alg_idx = alg_idx, valid = valid, pred = pred)
  }

  # Factorial design: all fold x algorithm combinations
  cvTasks <- expand.grid(fold_idx = seq_along(validRows), alg_idx = seq(k))
  cvTasksList <- split(cvTasks, seq(nrow(cvTasks)))
  cvTasksList <- cvTasksList[sample(length(cvTasksList))]

  time_train <- system.time({
    cvResults <- future.apply::future_lapply(
      cvTasksList, FUN = .singleCVTask,
      Y = Y, dataX = X, id = id,
      obsWeights = obsWeights, library = library,
      kScreen = kScreen, p = p,
      libraryNames = libraryNames,
      validRows = validRows,
      family = family,
      verbose = verbose,
      env = env,
      future.scheduling = Inf
    )

    for (result in cvResults) {
      Z[result$valid, result$alg_idx] <- result$pred
    }

    errorsInCVLibrary <- apply(Z, 2, function(x) any(is.na(x)))
    if (sum(errorsInCVLibrary) > 0) {
      n_cv_failed <- sum(errorsInCVLibrary)
      failed_cv_names <- libraryNames[as.logical(errorsInCVLibrary)]
      warning(sprintf(
        "ffSL: %d learner(s) failed in CV and will be dropped: %s.",
        n_cv_failed, paste(failed_cv_names, collapse = ", ")), call. = FALSE)
      Z[, as.logical(errorsInCVLibrary)] <- 0
    }
    if (all(Z == 0)) {
      stop("All algorithms dropped from library")
    }

    getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames,
                                   obsWeights = obsWeights, control = control,
                                   verbose = verbose,
                                   errorsInLibrary = errorsInCVLibrary)
    coef <- getCoef$coef
    names(coef) <- libraryNames
  })

  if (!("optimizer" %in% names(getCoef))) {
    getCoef["optimizer"] <- NA
  }

  m <- dim(newX)[1L]
  predY <- matrix(NA, nrow = m, ncol = k)

  .screenFun <- function(fun, list) {
    if (exists(fun, envir = env)) {
      screen_fn <- get(fun, envir = env)
    } else if (exists(fun, envir = asNamespace("SuperLearner"))) {
      screen_fn <- get(fun, envir = asNamespace("SuperLearner"))
    } else {
      stop(paste("Cannot find screening function:", fun))
    }

    out <- tryCatch(
      do.call(screen_fn, list),
      error = function(e) NULL
    )
    if (is.null(out)) {
      out <- rep(TRUE, ncol(list$X))
    }
    out
  }

  time_predict <- system.time({
    whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun,
                          list = list(Y = Y, X = X, family = family,
                                      id = id, obsWeights = obsWeights),
                          simplify = FALSE)
    whichScreen <- do.call(rbind, whichScreen)

    .predFun <- function(index, lib, Y, dataX, newX, whichScreen, family,
                          id, obsWeights, verbose, control, libraryNames, env) {
      out <- list(pred = NA, fitLibrary = NULL)
      pred_fn_name <- lib$predAlgorithm[index]

      if (exists(pred_fn_name, envir = env)) {
        pred_fn <- get(pred_fn_name, envir = env)
      } else if (exists(pred_fn_name, envir = asNamespace("SuperLearner"))) {
        pred_fn <- get(pred_fn_name, envir = asNamespace("SuperLearner"))
      } else {
        stop(paste("Cannot find function:", pred_fn_name))
      }

      testAlg <- tryCatch(
        do.call(pred_fn, list(
          Y = Y,
          X = dataX[, whichScreen[lib$rowScreen[index], ], drop = FALSE],
          newX = newX[, whichScreen[lib$rowScreen[index], ], drop = FALSE],
          family = family, id = id, obsWeights = obsWeights
        )),
        error = function(e) {
          message(sprintf("ffSL: %s failed on full data: %s",
                          libraryNames[index], e$message))
          NULL
        }
      )
      if (is.null(testAlg)) {
        out$pred <- rep.int(NA, times = nrow(newX))
      } else {
        out$pred <- testAlg$pred
        if (control$saveFitLibrary) {
          out$fitLibrary <- testAlg$fit
        }
      }
      if (verbose) {
        message(paste("full", libraryNames[index]))
      }
      invisible(out)
    }

    foo <- future.apply::future_lapply(
      seq(k), FUN = .predFun,
      lib = library$library, Y = Y, dataX = X, newX = newX,
      whichScreen = whichScreen, family = family, id = id,
      obsWeights = obsWeights, verbose = verbose, control = control,
      libraryNames = libraryNames, env = env,
      future.scheduling = Inf
    )
    predY <- do.call('cbind', lapply(foo, '[[', 'pred'))
    newFitLib <- lapply(foo, '[[', 'fitLibrary')
    # Replace NULL fitLibrary entries (failed learners) with SL.mean dummy
    # fits so predict.SuperLearner doesn't crash when re-predicting.
    # Must match SL.mean's fit structure: list(object = meanY) with class "SL.mean"
    # so predict.SL.mean dispatch works correctly.
    Y_mean <- mean(Y)
    for (j in seq_along(newFitLib)) {
      if (is.null(newFitLib[[j]])) {
        dummy_fit <- list(object = Y_mean)
        class(dummy_fit) <- "SL.mean"
        newFitLib[[j]] <- dummy_fit
      }
    }
    names(newFitLib) <- libraryNames
    assign('fitLibrary', newFitLib, envir = fitLibEnv)
    rm(foo)

    errorsInLibrary <- apply(predY, 2, function(xx) anyNA(xx))
    if (sum(errorsInLibrary) > 0) {
      n_pred_failed <- sum(errorsInLibrary)
      failed_pred_names <- libraryNames[as.logical(errorsInLibrary)]
      warning(sprintf(
        "ffSL: %d learner(s) failed in prediction: %s. Re-computing weights without them.",
        n_pred_failed, paste(failed_pred_names, collapse = ", ")), call. = FALSE)
      if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
        Z[, as.logical(errorsInLibrary)] <- 0
        if (all(Z == 0)) {
          stop("All algorithms dropped from library")
        }
        getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames,
                                       obsWeights = obsWeights, control = control,
                                       verbose = verbose,
                                       errorsInLibrary = errorsInLibrary)
        coef <- getCoef$coef
        names(coef) <- libraryNames
      } else {
        warning("coefficients already 0 for all failed algorithm(s)")
      }
    }

    getPred <- method$computePred(predY = predY, coef = coef, control = control)
  })

  colnames(predY) <- libraryNames

  if (sum(errorsInCVLibrary) > 0) {
    getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
  }

  time_end <- proc.time()

  times <- list(everything = time_end - time_start,
                train = time_train,
                predict = time_predict)

  out <- list(
    call = call,
    libraryNames = libraryNames,
    SL.library = library,
    SL.predict = getPred,
    coef = coef,
    library.predict = predY,
    Z = Z,
    cvRisk = getCoef$cvRisk,
    family = family,
    fitLibrary = get('fitLibrary', envir = fitLibEnv),
    id = id,
    varNames = varNames,
    validRows = validRows,
    method = method,
    whichScreen = whichScreen,
    control = control,
    errorsInCVLibrary = errorsInCVLibrary,
    errorsInLibrary = errorsInLibrary,
    obsWeights = obsWeights,
    metaOptimizer = getCoef$optimizer,
    env = env,
    times = times
  )
  class(out) <- c("SuperLearner")
  out
}
