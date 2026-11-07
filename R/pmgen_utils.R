# Internal utility functions

logit <- function(p) log(p / (1 - p))

invlogit <- function(x) 1 / (1 + exp(-x))

#' Compute the unit Fisher information matrix for a logistic model
#'
#' Computes the average per-observation Fisher information matrix from a design
#' matrix X and a linear predictor vector lp. The inverse of this matrix,
#' scaled by 1/N, gives the approximate sampling variance-covariance of the
#' coefficients at sample size N.
#'
#' @param X Numeric matrix. Design matrix including an intercept column (n x p).
#' @param lp Numeric vector. Linear predictor values (length n).
#' @return A p x p numeric matrix: the unit information matrix (uninverted).
#' @noRd
compute_unit_info_matrix <- function(X, lp) {
  # p*(1-p) = exp(lp) / (1 + exp(lp))^2 — logistic variance at each row
  w <- exp(lp) / (1 + exp(lp))^2

  p <- ncol(X)
  I_mat <- matrix(0, nrow = p, ncol = p)

  for (col in seq_len(p)) {
    for (row in col:p) {
      I_mat[row, col] <- mean(X[, row] * X[, col] * w)
      if (row != col) I_mat[col, row] <- I_mat[row, col]
    }
  }

  I_mat
}

#' Generate synthetic predictor data from distributional assumptions
#'
#' @param predictors Named list. Each element specifies a predictor distribution.
#'   Supported: `list(dist="normal", mean=, sd=)` and
#'   `list(dist="bernoulli", prob=)`.
#' @param n Integer. Number of rows to generate.
#' @param seed Integer or NULL.
#' @return A data frame with one column per predictor.
#' @noRd
generate_synth_data <- function(predictors, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  cols <- lapply(names(predictors), function(nm) {
    spec <- predictors[[nm]]
    if (is.null(spec$dist)) rlang::abort(paste0("Predictor '", nm, "' missing 'dist'."))
    switch(spec$dist,
      normal    = rnorm(n, mean = spec$mean %||% 0, sd = spec$sd %||% 1),
      bernoulli = rbinom(n, 1, prob = spec$prob %||% 0.5),
      rlang::abort(paste0("Unknown distribution '", spec$dist, "' for predictor '", nm, "'."))
    )
  })

  df <- as.data.frame(setNames(cols, names(predictors)))
  df
}

# Null-coalescing helper (base R doesn't have %||% outside rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Fit independent marginal distributions from a data frame.
# Binary columns (0/1 only) get Bernoulli; everything else gets Normal.
# Returns a named list suitable for generate_synth_data().
.fit_marginals <- function(df) {
  lapply(df, function(col) {
    vals <- unique(col[!is.na(col)])
    if (all(vals %in% c(0, 1))) {
      list(dist = "bernoulli", prob = mean(col, na.rm = TRUE))
    } else {
      list(dist = "normal", mean = mean(col, na.rm = TRUE), sd = sd(col, na.rm = TRUE))
    }
  })
}

# Is every non-NA value in a column 0 or 1? Same heuristic used by
# .fit_marginals() to decide Bernoulli vs. Normal, factored out so the
# "mvn" bisection path classifies columns identically.
.is_binary_col <- function(col) {
  vals <- unique(col[!is.na(col)])
  all(vals %in% c(0, 1))
}

# quickcstat via Wilcoxon-style case/control matching (from Pavlou, as used
# internally by .find_mu_sigma() for its Monte Carlo verification step).
.quickcstat <- function(y, pred, seed = 1) {
  set.seed(seed)
  casepred <- pred[y == 1]
  conpred  <- pred[y == 0]

  if (length(casepred) == 0 || length(conpred) == 0) {
    return(NA_real_)
  }

  if (length(conpred) > length(casepred)) {
    conpred  <- conpred[sample(length(conpred), length(casepred), replace = FALSE)]
    auc_true <- sum(casepred > conpred) / length(casepred)
  } else {
    casepred <- casepred[sample(length(casepred), length(conpred), replace = FALSE)]
    auc_true <- sum(casepred > conpred) / length(conpred)
  }

  auc_true
}

#' Find the mean and variance of a normal linear predictor giving a target
#' prevalence and c-statistic
#'
#' Assumes the linear predictor is normally distributed (as holds when
#' predictors are jointly multivariate normal and combined linearly). Solves
#' analytically/numerically for `(mu, sigma^2)` of that normal distribution
#' such that the implied outcome prevalence and c-statistic match the
#' targets, then verifies via a large Monte Carlo simulation.
#'
#' @param target_prev Target outcome prevalence.
#' @param target_c Target c-statistic (AUC).
#' @param seed Integer seed for the Monte Carlo verification draw.
#' @return Named numeric vector `c(mu, variance, <achieved c-stat>)`.
#' @noRd
.find_mu_sigma <- function(target_prev,
                            target_c,
                            min.opt = c(-10, 0),
                            max.opt = c(0.02, 5),
                            tol = 0.00001,
                            seed = 1) {

  pcfun <- function(x) {
    mean_lp     <- x[1]
    variance_lp <- x[2]

    f1 <- function(x) {
      stats::integrate(
        function(y) {
          stats::dnorm(x, mean = mean_lp, sd = sqrt(variance_lp)) *
            stats::dnorm(y, mean = mean_lp, sd = sqrt(variance_lp)) *
            (1 + exp(-x))^(-1) *
            (1 + exp(y))^(-1)
        },
        lower = -Inf,
        upper = x
      )$value
    }

    num <- stats::integrate(Vectorize(f1), lower = -Inf, upper = Inf)$value

    f2 <- function(x) {
      stats::integrate(
        function(y) {
          stats::dnorm(x, mean = mean_lp, sd = sqrt(variance_lp)) *
            stats::dnorm(y, mean = mean_lp, sd = sqrt(variance_lp)) *
            (1 + exp(-x))^(-1) *
            (1 + exp(y))^(-1)
        },
        lower = -Inf,
        upper = Inf
      )$value
    }

    denom <- stats::integrate(Vectorize(f2), lower = -Inf, upper = Inf)$value

    f3 <- function(x) {
      stats::dnorm(x, mean = mean_lp, sd = sqrt(variance_lp)) *
        (1 + exp(-x))^(-1)
    }

    c_hat <- num / denom

    prev <- stats::integrate(
      f3, lower = -Inf, upper = Inf, subdivisions = 1000L
    )$value

    abs(c_hat - target_c)^2 + abs(prev - target_prev)^2
  }

  if (target_c > 0.65) {
    out <- stats::optim(
      par = c(-2.65, 0.1),
      pcfun,
      c(min.opt, max.opt, tol = tol))$par
  } else {
    sigma_c <- sqrt(2) * stats::qnorm(target_c)
    mu <- 0.5 * (2 * target_prev - 1) * (sigma_c^2) +
      log(target_prev / (1 - target_prev))
    sigma <- sqrt((sigma_c^2) *
                    (1 + target_prev * (1 - target_prev) * (sigma_c^2)))
    out <- c(mu, sigma^2)
  }

  N    <- 2000000
  lp   <- stats::rnorm(N, mean = out[1], sd = sqrt(out[2]))
  p    <- (1 + exp(-lp))^(-1)
  y    <- stats::rbinom(N, 1, prob = p)
  c_achieved <- .quickcstat(y, lp, seed = seed)

  c(mu = out[1], variance = out[2], cstat = c_achieved)
}

#' Generate a synthetic population from a multivariate-normal approximation
#' to `joint_data`, preserving its correlation structure
#'
#' Draws `n` rows from a multivariate normal distribution parameterised by
#' the means, SDs, and correlation matrix estimated from `pop_data`, then
#' dichotomizes any binary (0/1) columns via a threshold on the
#' corresponding generated normal column chosen to match that column's
#' empirical proportion in `pop_data`. This is a standard Gaussian-copula
#' (NORTA) style approximation: it preserves the target marginal proportions
#' for binary columns, but the resulting correlation among dichotomized
#' columns is mildly attenuated relative to `corr_mat` (no correction is
#' applied for this).
#'
#' @param pop_data Data frame of real predictor data to approximate.
#' @param n Integer. Number of synthetic rows to generate.
#' @param seed Integer or NULL.
#' @return A list with `data` (the generated data frame, same columns/order
#'   as `pop_data`), `means`, `predictor_sds`, `corr_mat`, `sigma_x`
#'   (covariance matrix), and `is_binary` (named logical vector).
#' @noRd
generate_mvn_synth_data <- function(pop_data, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  predictors    <- names(pop_data)
  X             <- as.matrix(pop_data)
  means         <- colMeans(X)
  predictor_sds <- apply(X, 2, sd)
  corr_mat      <- stats::cor(X)
  is_binary     <- vapply(pop_data, .is_binary_col, logical(1))

  sd_mat  <- diag(predictor_sds, length(predictors))
  sigma_x <- sd_mat %*% corr_mat %*% sd_mat
  dimnames(sigma_x) <- list(predictors, predictors)

  z <- mvtnorm::rmvnorm(n = n, mean = means, sigma = sigma_x)
  colnames(z) <- predictors

  for (nm in predictors[is_binary]) {
    prop      <- unname(means[nm])
    threshold <- stats::qnorm(1 - prop, mean = means[nm], sd = predictor_sds[nm])
    z[, nm]   <- as.numeric(z[, nm] > threshold)
  }

  list(
    data          = as.data.frame(z),
    means         = means,
    predictor_sds = predictor_sds,
    corr_mat      = corr_mat,
    sigma_x       = sigma_x,
    is_binary     = is_binary
  )
}
