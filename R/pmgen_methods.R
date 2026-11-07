# Built-in method implementations for pmgen_simulate()
#
# Each fit_*() function takes:
#   train_data  — data frame with predictor columns and outcome column
#   predictors  — character vector of predictor column names
#   outcome_col — name of the binary outcome column in train_data
#   method_args — named list of tuning overrides
#
# Each returns a *prediction closure*: function(new_data) -> numeric vector
# of predicted probabilities, one per row of new_data. The built-in closures
# also attach the linear predictor as attribute "lp", which is used for the
# calibration slope and intercept (exact even where probabilities round to
# 0 or 1).


# ---- Logistic regression (unpenalised) ----------------------------------------

#' @noRd
fit_logistic <- function(train_data, predictors, outcome_col, method_args = list()) {
  separation <- match.arg(method_args$separation %||% "drop", c("drop", "keep"))

  # Perfect-prediction handling: drop predictors (and rows) that predict the
  # outcome perfectly before fitting; they get coefficient 0 when predicting
  dropped <- character(0)
  if (separation == "drop") {
    pp         <- .drop_perfect_predictors(train_data, predictors, outcome_col)
    train_data <- pp$data
    dropped    <- pp$dropped
  }
  kept <- setdiff(predictors, dropped)

  rhs  <- if (length(kept) > 0) paste(kept, collapse = " + ") else "1"
  fmla <- stats::as.formula(paste(outcome_col, "~", rhs))
  fit <- stats::glm(fmla, data = train_data, family = stats::binomial(),
                    control = stats::glm.control(maxit = 100))

  # Full coefficient vector over all predictors: dropped and aliased
  # (NA) coefficients are set to 0
  fit_coefs <- stats::coef(fit)
  fit_coefs[is.na(fit_coefs)] <- 0
  coefs <- stats::setNames(rep(0, length(predictors) + 1L),
                           c("(Intercept)", predictors))
  coefs[names(fit_coefs)] <- fit_coefs

  # Optional van Houwelingen uniform shrinkage factor
  shrink <- isTRUE(method_args$shrinkage)
  if (shrink) {
    null_dev  <- fit$null.deviance
    resid_dev <- fit$deviance
    # number of predictor coefficients actually estimated
    df_pred   <- sum(!is.na(stats::coef(fit)[-1]))
    # heuristic: (chi^2 - df) / chi^2, floored at 0
    chi2  <- null_dev - resid_dev
    factor <- if (df_pred == 0) 1 else max(0, (chi2 - df_pred) / chi2)
    # shrink predictor coefficients (not intercept) toward zero
    coefs[-1] <- coefs[-1] * factor
  }

  pred_fn <- function(new_data) {
    lp <- coefs[1] + as.vector(as.matrix(new_data[, predictors, drop = FALSE]) %*% coefs[-1])
    structure(stats::plogis(lp), lp = lp)
  }
  attr(pred_fn, "dropped") <- dropped
  pred_fn
}


#' Drop predictors that predict a binary outcome perfectly
#'
#' For each two-valued predictor (in order), checks whether one of its levels
#' has a constant outcome. If so, that predictor and the rows at that level are
#' removed, and the check starts again on the remaining data, until no such
#' predictor remains. A predictor with more than two values that separates the
#' outcome completely (some cut-point splits it exactly) is an error, as is an
#' outcome left with no variation.
#'
#' @param data Data frame with predictor columns and outcome column.
#' @param predictors Character vector of predictor column names.
#' @param outcome_col Name of the binary outcome column in `data`.
#' @return A list with `data` (the remaining rows) and `dropped` (names of the
#'   dropped predictors, in the order they were dropped).
#' @noRd
.drop_perfect_predictors <- function(data, predictors, outcome_col) {
  dropped <- character(0)

  repeat {
    hit <- FALSE
    for (x in setdiff(predictors, dropped)) {
      lvls <- sort(unique(data[[x]]))
      if (length(lvls) != 2L) next
      for (lv in lvls) {
        at_lv <- data[[x]] == lv
        if (length(unique(data[[outcome_col]][at_lv])) == 1L) {
          data    <- data[!at_lv, , drop = FALSE]
          dropped <- c(dropped, x)
          hit     <- TRUE
          break
        }
      }
      if (hit) break
    }
    if (!hit) break
  }

  y <- data[[outcome_col]]
  if (length(unique(y)) < 2L) {
    stop("Outcome does not vary after removing perfectly predicted observations.")
  }

  # Complete separation on a predictor with more than two values
  for (x in setdiff(predictors, dropped)) {
    v <- data[[x]]
    if (length(unique(v)) <= 2L) next
    if (max(v[y == 0]) < min(v[y == 1]) || max(v[y == 1]) < min(v[y == 0])) {
      stop(paste0("Predictor '", x, "' predicts the outcome perfectly."))
    }
  }

  list(data = data, dropped = dropped)
}


# ---- Ridge regression (glmnet, alpha = 0) -------------------------------------

#' @noRd
fit_ridge <- function(train_data, predictors, outcome_col, method_args = list()) {
  X   <- as.matrix(train_data[, predictors, drop = FALSE])
  y   <- train_data[[outcome_col]]
  nfolds   <- method_args$nfolds   %||% 10
  type_msr <- method_args$type.measure %||% "deviance"
  lambda   <- method_args$lambda   %||% NULL

  if (is.null(lambda)) {
    cv_fit <- glmnet::cv.glmnet(X, y, alpha = 0, family = "binomial",
                                nfolds = nfolds, type.measure = type_msr)
    lambda_use <- cv_fit$lambda.min
    fit <- cv_fit$glmnet.fit
  } else {
    fit        <- glmnet::glmnet(X, y, alpha = 0, family = "binomial")
    lambda_use <- lambda
  }

  function(new_data) {
    X_new <- as.matrix(new_data[, predictors, drop = FALSE])
    lp <- as.vector(stats::predict(fit, newx = X_new, s = lambda_use, type = "link"))
    structure(stats::plogis(lp), lp = lp)
  }
}


# ---- Lasso regression (glmnet, alpha = 1) -------------------------------------

#' @noRd
fit_lasso <- function(train_data, predictors, outcome_col, method_args = list()) {
  X   <- as.matrix(train_data[, predictors, drop = FALSE])
  y   <- train_data[[outcome_col]]
  nfolds   <- method_args$nfolds   %||% 10
  type_msr <- method_args$type.measure %||% "deviance"
  lambda   <- method_args$lambda   %||% NULL

  if (is.null(lambda)) {
    cv_fit <- glmnet::cv.glmnet(X, y, alpha = 1, family = "binomial",
                                nfolds = nfolds, type.measure = type_msr)
    lambda_use <- cv_fit$lambda.min
    fit <- cv_fit$glmnet.fit
  } else {
    fit        <- glmnet::glmnet(X, y, alpha = 1, family = "binomial")
    lambda_use <- lambda
  }

  function(new_data) {
    X_new <- as.matrix(new_data[, predictors, drop = FALSE])
    lp <- as.vector(stats::predict(fit, newx = X_new, s = lambda_use, type = "link"))
    structure(stats::plogis(lp), lp = lp)
  }
}


# ---- Elastic net (glmnet, user-supplied alpha) ---------------------------------

#' @noRd
fit_elasticnet <- function(train_data, predictors, outcome_col, method_args = list()) {
  alpha <- method_args$alpha
  if (is.null(alpha)) {
    rlang::abort("method_args$alpha must be specified for elastic net (0 < alpha < 1).")
  }
  if (alpha <= 0 || alpha >= 1) {
    rlang::abort("method_args$alpha must be strictly between 0 and 1 for elastic net.")
  }

  X   <- as.matrix(train_data[, predictors, drop = FALSE])
  y   <- train_data[[outcome_col]]
  nfolds   <- method_args$nfolds   %||% 10
  type_msr <- method_args$type.measure %||% "deviance"
  lambda   <- method_args$lambda   %||% NULL

  if (is.null(lambda)) {
    cv_fit <- glmnet::cv.glmnet(X, y, alpha = alpha, family = "binomial",
                                nfolds = nfolds, type.measure = type_msr)
    lambda_use <- cv_fit$lambda.min
    fit <- cv_fit$glmnet.fit
  } else {
    fit        <- glmnet::glmnet(X, y, alpha = alpha, family = "binomial")
    lambda_use <- lambda
  }

  function(new_data) {
    X_new <- as.matrix(new_data[, predictors, drop = FALSE])
    lp <- as.vector(stats::predict(fit, newx = X_new, s = lambda_use, type = "link"))
    structure(stats::plogis(lp), lp = lp)
  }
}


# ---- Firth's penalised logistic regression ------------------------------------

#' @noRd
fit_firth <- function(train_data, predictors, outcome_col, method_args = list()) {
  fmla <- stats::as.formula(
    paste(outcome_col, "~", paste(predictors, collapse = " + "))
  )
  fit <- logistf::logistf(fmla, data = train_data, control = logistf::logistf.control(
    maxit  = method_args$maxit  %||% 100,
    maxstep = method_args$maxstep %||% 10
  ))

  coefs <- stats::coef(fit)

  function(new_data) {
    X_new <- as.matrix(new_data[, predictors, drop = FALSE])
    lp    <- coefs[1] + as.vector(X_new %*% coefs[-1])
    structure(stats::plogis(lp), lp = lp)
  }
}


# ---- Method dispatcher --------------------------------------------------------

#' Resolve a method specification to a fit function
#'
#' @param method Character string naming a built-in method, or a user-supplied
#'   function with signature `(train_data, predictors, outcome_col, method_args)`.
#' @return A function with that same signature returning a prediction closure.
#'
#' Built-in methods that need a suggested package check for it here, once,
#' because fit errors are caught per simulation: a missing package would
#' otherwise show up only as failed fits.
#' @noRd
resolve_method <- function(method) {
  if (is.function(method)) return(method)

  needs <- switch(method,
    ridge = , lasso = , elasticnet = "glmnet",
    firth = "logistf",
    NULL
  )
  if (!is.null(needs)) {
    rlang::check_installed(needs, reason = paste0("to use method = '", method, "'."))
  }

  switch(method,
    logistic   = fit_logistic,
    ridge      = fit_ridge,
    lasso      = fit_lasso,
    elasticnet = fit_elasticnet,
    firth      = fit_firth,
    rlang::abort(paste0(
      "Unknown method '", method, "'. ",
      "Built-in options: 'logistic', 'ridge', 'lasso', 'elasticnet', 'firth'. ",
      "Or supply a function with signature (train_data, predictors, outcome_col, method_args)."
    ))
  )
}
