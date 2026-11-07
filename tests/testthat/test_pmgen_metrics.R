# ---- Calibration from the linear predictor -----------------------------------

# Target population where some linear predictors are large enough that
# plogis() rounds to exactly 1
make_saturated_tp <- function() {
  set.seed(4)
  n  <- 5000
  lp <- c(rnorm(n - 50), runif(50, 50, 60))
  list(data = data.frame(lp = lp, p_true = plogis(lp),
                         outcome = rbinom(n, 1, plogis(lp))))
}

# (glm() warns "fitted probabilities numerically 0 or 1" on this data)
test_that("calibration slope and intercept use the linear predictor when given", suppressWarnings({
  tp <- make_saturated_tp()
  lp <- 0.8 * tp$data$lp + 0.1
  p  <- plogis(lp)
  expect_true(any(p == 1))

  m <- compute_metrics_one_draw(p, tp, c("cal_slope", "cal_intercept"), lp = lp)
  y <- tp$data$outcome
  expect_equal(unname(m["cal_slope"]),
               unname(coef(glm(y ~ lp, family = binomial()))[2]))
  expect_equal(unname(m["cal_intercept"]),
               unname(coef(glm(y ~ 1, offset = lp, family = binomial()))[1]))
}))

test_that("without a linear predictor, probabilities of 0 or 1 are clipped", suppressWarnings({
  tp <- make_saturated_tp()
  lp <- 0.8 * tp$data$lp + 0.1
  p  <- plogis(lp)

  m_lp   <- compute_metrics_one_draw(p, tp, c("cal_slope", "cal_intercept"), lp = lp)
  m_clip <- compute_metrics_one_draw(p, tp, c("cal_slope", "cal_intercept"))
  expect_false(anyNA(m_clip))
  # clipping caps the linear predictor at qlogis(1 - 1e-15), about 34.5
  expect_equal(m_clip, m_lp, tolerance = 1e-3)

  # where nothing saturates, clipping changes nothing
  ok <- p > 1e-15 & p < 1 - 1e-15
  tp_ok <- list(data = tp$data[ok, ])
  expect_equal(compute_metrics_one_draw(p[ok], tp_ok, "cal_slope"),
               compute_metrics_one_draw(p[ok], tp_ok, "cal_slope", lp = lp[ok]))
}))

test_that("built-in logistic predictions carry the linear predictor", {
  df <- data.frame(x = c(-2, -1, 0, 1, 2, 0.5, -0.5, 1.5),
                   y = c(0, 0, 1, 0, 1, 1, 0, 1))
  pred_fn <- fit_logistic(df, "x", "y")
  p <- pred_fn(data.frame(x = c(-1, 0, 50)))
  fit <- glm(y ~ x, family = binomial(), data = df)
  expect_equal(attr(p, "lp"),
               unname(predict(fit, newdata = data.frame(x = c(-1, 0, 50)))))
  expect_equal(as.vector(p), plogis(attr(p, "lp")))
})

test_that("pmgen_simulate() works with a custom method that returns only probabilities", {
  fx <- suppressMessages(make_sim_fixtures())
  my_method <- function(train_data, predictors, outcome_col, method_args) {
    fit <- glm(reformulate(predictors, outcome_col), family = binomial(),
               data = train_data)
    function(new_data) unname(predict(fit, newdata = new_data, type = "response"))
  }
  res_custom <- suppressMessages(
    pmgen_simulate(fx$ref, fx$tp, N = 300L, n_sims = 5L, method = my_method,
                 metrics = c("cal_slope", "cal_intercept"), seed = 8,
                 verbose = FALSE)
  )
  res_builtin <- suppressMessages(
    pmgen_simulate(fx$ref, fx$tp, N = 300L, n_sims = 5L, method = "logistic",
                 metrics = c("cal_slope", "cal_intercept"), seed = 8,
                 verbose = FALSE)
  )
  # same fits; at N = 300 nothing saturates, so the two routes agree
  expect_equal(res_custom$posterior[["300"]], res_builtin$posterior[["300"]])
})
