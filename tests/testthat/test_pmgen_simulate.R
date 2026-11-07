suppressMessages({
  fx <- make_sim_fixtures()
})
ref <- fx$ref
tp  <- fx$tp

# ---- Target population partitioning -----------------------------------------

test_that("pmgen_target_population() produces a non-overlapping partition", {
  # joint_data has 10000 rows; target pop takes 2000 without replacement
  expect_equal(nrow(tp$data), 2000L)
  expect_equal(nrow(tp$dev_pool), 8000L)
  expect_equal(nrow(tp$data) + nrow(tp$dev_pool),
               nrow(tp$ref_model$joint_data))
})


test_that("pmgen_target_population() errors when n >= nrow(joint_data)", {
  expect_error(
    pmgen_target_population(ref, n = 10000L, seed = 1),
    "strictly less than"
  )
})

test_that("dev_pool_n caps the development pool size", {
  tp_capped <- pmgen_target_population(ref, n = 2000L,
                                     dev_pool_n = 3000L, seed = 1)
  expect_equal(nrow(tp_capped$dev_pool), 3000L)
  expect_equal(nrow(tp_capped$data), 2000L)
})

# ---- Smoke test: logistic regression -----------------------------------------

test_that("pmgen_simulate() returns a pmgen_result for logistic method", {
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = c(50L, 100L),
                 n_sims  = 10L,
                 method  = "logistic",
                 metrics = c("cstat", "cal_slope", "mape"),
                 seed    = 42,
                 verbose = FALSE)
  )

  expect_s3_class(result, "pmgen_result")
  expect_equal(result$N, c(50L, 100L))
  expect_equal(result$method, "logistic")
  expect_equal(result$metrics, c("cstat", "cal_slope", "mape"))
})

test_that("posterior dimensions are correct", {
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = c(80L),
                 n_sims  = 15L,
                 method  = "logistic",
                 metrics = c("cstat", "cal_slope"),
                 seed    = 7,
                 verbose = FALSE)
  )

  post <- result$posterior[["80"]]
  expect_true(is.data.frame(post))
  expect_equal(nrow(post), 15L)
  expect_equal(sort(names(post)), sort(c("cstat", "cal_slope")))
})

test_that("summary.pmgen_result() returns a data frame with correct columns", {
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = c(60L, 120L),
                 n_sims  = 10L,
                 method  = "logistic",
                 metrics = c("cstat", "mape"),
                 seed    = 3,
                 verbose = FALSE)
  )

  s <- summary(result)
  expect_true(is.data.frame(s))
  expect_true("N" %in% names(s))
  expect_true("metric" %in% names(s))
  expect_equal(sort(unique(s$N)), c(60L, 120L))
})


# ---- Built-in methods: ridge, lasso, elasticnet, firth -----------------------

test_that("pmgen_simulate() works with ridge method", {
  skip_if_not_installed("glmnet")
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 80L, n_sims = 5L, method = "ridge",
                 metrics = c("cstat"), seed = 10, verbose = FALSE)
  )
  expect_s3_class(result, "pmgen_result")
  expect_equal(result$method, "ridge")
  expect_false(all(is.na(result$posterior[["80"]]$cstat)))
})

test_that("pmgen_simulate() works with lasso method", {
  skip_if_not_installed("glmnet")
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 80L, n_sims = 5L, method = "lasso",
                 metrics = c("cstat"), seed = 11, verbose = FALSE)
  )
  expect_s3_class(result, "pmgen_result")
  expect_equal(result$method, "lasso")
})

test_that("pmgen_simulate() works with elasticnet method", {
  skip_if_not_installed("glmnet")
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 80L, n_sims = 5L, method = "elasticnet",
                 method_args = list(alpha = 0.5),
                 metrics = c("cstat"), seed = 12, verbose = FALSE)
  )
  expect_s3_class(result, "pmgen_result")
  expect_equal(result$method, "elasticnet")
})

test_that("pmgen_simulate() works with firth method", {
  skip_if_not_installed("logistf")
  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 80L, n_sims = 5L, method = "firth",
                 metrics = c("cstat"), seed = 13, verbose = FALSE)
  )
  expect_s3_class(result, "pmgen_result")
  expect_equal(result$method, "firth")
})


# ---- Custom method interface --------------------------------------------------

test_that("pmgen_simulate() accepts a custom method function", {
  my_method <- function(train_data, predictors, outcome_col, method_args) {
    fmla <- stats::as.formula(
      paste(outcome_col, "~", paste(predictors, collapse = "+"))
    )
    fit <- stats::glm(fmla, data = train_data, family = stats::binomial())
    function(new_data) stats::predict(fit, newdata = new_data, type = "response")
  }

  result <- suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 60L, n_sims = 5L, method = my_method,
                 metrics = c("cstat"), seed = 20, verbose = FALSE)
  )
  expect_s3_class(result, "pmgen_result")
  expect_equal(result$method, "custom")
})


# ---- Error handling ----------------------------------------------------------

test_that("pmgen_simulate() errors on unknown method string", {
  expect_error(
    suppressMessages(
      pmgen_simulate(ref, tp, N = 50L, n_sims = 2L,
                   method = "neuralnet", verbose = FALSE)
    ),
    "Unknown method"
  )
})

test_that("pmgen_simulate() errors on invalid metric name", {
  expect_error(
    suppressMessages(
      pmgen_simulate(ref, tp, N = 50L, n_sims = 2L,
                   metrics = c("cstat", "banana"), verbose = FALSE)
    ),
    "Unknown metric"
  )
})

test_that("pmgen_simulate() errors on bad threshold", {
  expect_error(
    suppressMessages(
      pmgen_simulate(ref, tp, N = 50L, n_sims = 2L,
                   threshold = 1.5, verbose = FALSE)
    ),
    "threshold"
  )
})

test_that("pmgen_simulate() warns when N exceeds dev pool size", {
  # dev pool has ~8000 rows; N = 1e6 should trigger the warning
  expect_warning(
    suppressMessages(
      pmgen_simulate(ref, tp,
                   N = 1e6L, n_sims = 2L,
                   method = "logistic",
                   metrics = "cstat",
                   seed = 1, verbose = FALSE)
    ),
    "with replacement"
  )
})


# ---- Reproducibility ----------------------------------------------------------

test_that("pmgen_simulate() is reproducible with same seed", {
  run <- function() suppressMessages(
    pmgen_simulate(ref, tp,
                 N = 60L, n_sims = 8L, method = "logistic",
                 metrics = c("cstat"), seed = 99, verbose = FALSE)
  )
  r1 <- run()
  r2 <- run()
  expect_equal(r1$posterior[["60"]]$cstat,
               r2$posterior[["60"]]$cstat)
})


# ---- Integration / regression test against paper results ---------------------

test_that("DVT full workflow at N=400 reproduces expected performance range", {
  skip_on_cran()

  data(dvt_data, package = "pmsampsize")
  predictors <- c("male", "OCuse", "malignancy", "surgery", "ab_leg_trauma",
                  "vein_distension", "calf_diff", "ab_ddimer",
                  "pregnancy", "immobile", "clotting")

  ref_dvt <- suppressMessages(
    pmgen_reference(
      joint_data  = dvt_data,
      predictors  = predictors,
      prevalence  = 0.21,
      cstatistic  = 0.78,
      tolerance   = 0.005,
      seed        = 42
    )
  )

  tp_dvt <- suppressMessages(
    pmgen_target_population(ref_dvt, n = 5000L, seed = 42)
  )

  res <- suppressMessages(
    pmgen_simulate(
      ref_model  = ref_dvt,
      target_pop = tp_dvt,
      N          = 400L,
      n_sims     = 500L,
      method     = "logistic",
      metrics    = c("cstat", "cal_slope"),
      seed       = 42,
      verbose    = FALSE
    )
  )

  post <- res$posterior[["400"]]

  mean_cal_slope <- mean(post$cal_slope, na.rm = TRUE)
  mean_cstat_deg <- mean(res$summaries[["400"]]$degradation_mean[
    res$summaries[["400"]]$metric == "cstat"
  ])

  expect_true(
    mean_cal_slope >= 0.75 && mean_cal_slope <= 1.15,
    info = sprintf("Mean cal slope %.3f outside [0.75, 1.15]", mean_cal_slope)
  )

  expect_true(
    mean_cstat_deg < 0.05,
    info = sprintf("Mean cstat degradation %.3f >= 0.05", mean_cstat_deg)
  )
})
