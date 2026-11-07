# ---- Perfect-prediction handling in fit_logistic() ----------------------------

# x1 continuous; b1 has a level (b1 = 1) where everyone has the event;
# b2 is fine on the full data
make_pp_data <- function() {
  set.seed(2)
  n  <- 60
  df <- data.frame(x1 = rnorm(n), b1 = rbinom(n, 1, 0.2), b2 = rbinom(n, 1, 0.5))
  df$y <- rbinom(n, 1, 0.5)
  df$y[df$b1 == 1] <- 1L
  df
}

test_that(".drop_perfect_predictors() drops the predictor and its rows", {
  df  <- make_pp_data()
  out <- .drop_perfect_predictors(df, c("x1", "b1", "b2"), "y")
  expect_equal(out$dropped, "b1")
  expect_equal(nrow(out$data), sum(df$b1 == 0))
  expect_true(all(out$data$b1 == 0))
})

test_that(".drop_perfect_predictors() repeats until no predictor is perfect", {
  df <- make_pp_data()
  # b2 = 1 has a constant outcome only once the b1 = 1 rows are removed
  df$y[df$b1 == 0 & df$b2 == 1] <- 0L
  expect_true(length(unique(df$y[df$b2 == 1])) == 2L)
  out <- .drop_perfect_predictors(df, c("x1", "b1", "b2"), "y")
  expect_equal(out$dropped, c("b1", "b2"))
  expect_true(all(out$data$b1 == 0 & out$data$b2 == 0))
})

test_that(".drop_perfect_predictors() leaves data without perfect prediction alone", {
  set.seed(3)
  df <- data.frame(x1 = rnorm(200), b1 = rbinom(200, 1, 0.5))
  df$y <- rbinom(200, 1, 0.5)
  out <- .drop_perfect_predictors(df, c("x1", "b1"), "y")
  expect_equal(out$dropped, character(0))
  expect_identical(out$data, df)
})

test_that(".drop_perfect_predictors() errors on complete separation by a continuous predictor", {
  df <- data.frame(x1 = c(-3, -2, -1, 1, 2, 3), b1 = c(0, 1, 0, 1, 0, 1),
                   y = c(0, 0, 0, 1, 1, 1))
  expect_error(.drop_perfect_predictors(df, c("x1", "b1"), "y"),
               "'x1' predicts the outcome perfectly")
  # overlap at one value is not complete separation
  df$x1[4] <- -1
  expect_silent(.drop_perfect_predictors(df, c("x1", "b1"), "y"))
})

test_that(".drop_perfect_predictors() errors if the outcome no longer varies", {
  # b1 separates the outcome completely: both levels have a constant outcome
  df <- data.frame(x1 = rnorm(8), b1 = rep(0:1, 4), y = rep(0:1, 4))
  expect_error(.drop_perfect_predictors(df, c("x1", "b1"), "y"),
               "Outcome does not vary")
})

test_that("fit_logistic() with separation = 'drop' matches glm() on the reduced data", {
  df <- make_pp_data()
  pred_fn <- fit_logistic(df, c("x1", "b1", "b2"), "y")
  expect_equal(attr(pred_fn, "dropped"), "b1")

  fit_ref <- glm(y ~ x1 + b2, family = binomial(), data = df[df$b1 == 0, ])
  new <- data.frame(x1 = c(-1, 0, 1), b1 = c(1, 0, 1), b2 = c(0, 1, 1))
  # b1 has coefficient 0, so its value does not affect predictions
  expect_equal(as.vector(pred_fn(new)),
               unname(predict(fit_ref, newdata = new, type = "response")))
})

test_that("fit_logistic() is the default and 'drop' is used when separation is unset", {
  df <- make_pp_data()
  new <- data.frame(x1 = 0.5, b1 = 1, b2 = 0)
  expect_equal(fit_logistic(df, c("x1", "b1", "b2"), "y")(new),
               fit_logistic(df, c("x1", "b1", "b2"), "y",
                            list(separation = "drop"))(new))
})

test_that("fit_logistic() with separation = 'keep' is plain glm()", {
  df <- make_pp_data()
  pred_fn <- suppressWarnings(
    fit_logistic(df, c("x1", "b1", "b2"), "y", list(separation = "keep"))
  )
  expect_equal(attr(pred_fn, "dropped"), character(0))
  fit_ref <- suppressWarnings(glm(y ~ x1 + b1 + b2, family = binomial(),
                                  data = df,
                                  control = glm.control(maxit = 100)))
  new <- data.frame(x1 = c(-1, 0, 1), b1 = c(1, 0, 1), b2 = c(0, 1, 1))
  expect_equal(as.vector(pred_fn(new)),
               unname(predict(fit_ref, newdata = new, type = "response")))
})

test_that("fit_logistic() errors on an unknown separation value", {
  df <- make_pp_data()
  expect_error(fit_logistic(df, c("x1", "b1", "b2"), "y",
                            list(separation = "firth")))
})

test_that("fit_logistic() shrinkage counts only the predictors kept", {
  df <- make_pp_data()
  pred_fn <- fit_logistic(df, c("x1", "b1", "b2"), "y",
                          list(shrinkage = TRUE))

  fit_ref <- glm(y ~ x1 + b2, family = binomial(), data = df[df$b1 == 0, ])
  chi2    <- fit_ref$null.deviance - fit_ref$deviance
  factor  <- max(0, (chi2 - 2) / chi2)
  b       <- coef(fit_ref)
  new     <- data.frame(x1 = c(-1, 0, 1), b1 = c(1, 0, 1), b2 = c(0, 1, 1))
  expect_equal(as.vector(pred_fn(new)),
               plogis(b[1] + factor * (b["x1"] * new$x1 + b["b2"] * new$b2)))
})

test_that("fit_logistic() fits an intercept-only model if every predictor is dropped", {
  df <- data.frame(b1 = c(1, 1, 0, 0, 0, 0), y = c(1, 1, 0, 1, 0, 1))
  pred_fn <- fit_logistic(df, "b1", "y", list(shrinkage = TRUE))
  expect_equal(attr(pred_fn, "dropped"), "b1")
  expect_equal(as.vector(pred_fn(data.frame(b1 = c(0, 1)))), c(0.5, 0.5))
})


# ---- Reporting in pmgen_simulate() ---------------------------------------------

suppressMessages({
  fx <- make_sim_fixtures()
})

test_that("pmgen_simulate() records perfect-prediction drops and failed fits", {
  # x3 is binary; at N = 15 one of its levels often has a constant outcome.
  # (glm() also warns about near-separation in some of these small fits.)
  res <- suppressWarnings(suppressMessages(
    pmgen_simulate(fx$ref, fx$tp, N = c(15L, 200L), n_sims = 40L,
                 method = "logistic", metrics = c("cstat", "cal_slope"),
                 seed = 5, verbose = FALSE)
  ))
  fi <- res$fit_info
  expect_equal(fi$N, c(15L, 200L))
  expect_equal(fi$n_sims, c(40L, 40L))
  expect_true(fi$n_dropped[1] > 0)
  expect_equal(fi$n_dropped[2], 0L)
  # failed fits (complete separation on x1 or x2) give NA metrics
  expect_equal(fi$n_failed, vapply(res$posterior, function(p) sum(is.na(p$cstat)),
                                   integer(1), USE.NAMES = FALSE))
  expect_equal(dim(res$dropped_predictors), c(2L, 3L))
  expect_equal(colnames(res$dropped_predictors), c("x1", "x2", "x3"))
  expect_equal(unname(res$dropped_predictors[, "x3"]), fi$n_dropped)
  expect_output(print(res), "Perfect prediction")
})

test_that("pmgen_simulate() separation = 'keep' drops nothing", {
  res <- suppressWarnings(suppressMessages(
    pmgen_simulate(fx$ref, fx$tp, N = 15L, n_sims = 40L, method = "logistic",
                 method_args = list(separation = "keep"),
                 metrics = "cstat", seed = 5, verbose = FALSE)
  ))
  expect_equal(res$fit_info$n_dropped, 0L)
})

test_that("pmgen_simulate() reports perfect prediction when verbose", {
  suppressWarnings(expect_message(
    pmgen_simulate(fx$ref, fx$tp, N = 15L, n_sims = 40L, method = "logistic",
                 metrics = "cstat", seed = 5, verbose = TRUE),
    "Perfect prediction in"
  ))
})

test_that("pmgen_simulate() checks method_args$separation", {
  expect_error(
    pmgen_simulate(fx$ref, fx$tp, N = 50L, n_sims = 2L,
                 method_args = list(separation = "omit"), verbose = FALSE),
    "must be \"drop\" or \"keep\""
  )
  skip_if_not_installed("glmnet")
  expect_warning(
    pmgen_simulate(fx$ref, fx$tp, N = 80L, n_sims = 2L, method = "lasso",
                 method_args = list(separation = "keep"),
                 metrics = "cstat", seed = 1, verbose = FALSE),
    "ignored"
  )
})

test_that("pmgen_simulate() warns about failed fits with the first error", {
  bad_method <- function(train_data, predictors, outcome_col, method_args) {
    stop("deliberate failure")
  }
  expect_warning(
    res <- pmgen_simulate(fx$ref, fx$tp, N = 50L, n_sims = 3L,
                        method = bad_method, metrics = "cstat",
                        seed = 1, verbose = FALSE),
    "3 of 3.*deliberate failure"
  )
  expect_equal(res$fit_info$n_failed, 3L)
  expect_true(all(is.na(res$posterior[["50"]]$cstat)))
})
