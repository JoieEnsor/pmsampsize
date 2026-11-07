sim <- function(fx, N, metrics, seed, method = "logistic", method_args = list()) {
  suppressMessages(pmgen_simulate(
    fx$ref, fx$tp, N = N, n_sims = 20L, method = method,
    method_args = method_args, metrics = metrics, seed = seed,
    verbose = FALSE
  ))
}

test_that("pmgen_compare errors on non-pmgen_result input", {
  fx <- make_sim_fixtures()
  r1 <- sim(fx, 200L, "cstat", 1)
  expect_error(pmgen_compare(list(), r1), "pmgen_result")
  expect_error(pmgen_compare(r1, list()), "pmgen_result")
})

test_that("pmgen_compare errors when results share no N values", {
  fx <- make_sim_fixtures()
  r1 <- sim(fx, 200L, "cstat", 1)
  r2 <- sim(fx, 500L, "cstat", 2)
  expect_error(pmgen_compare(r1, r2), "share no N values")
})

test_that("pmgen_compare errors when results share no metrics", {
  fx <- make_sim_fixtures()
  r1 <- sim(fx, 200L, "cstat", 1)
  r2 <- sim(fx, 200L, "mape", 2)
  expect_error(pmgen_compare(r1, r2), "share no metrics")
})

test_that("pmgen_compare returns correct structure and values", {
  fx <- make_sim_fixtures()
  r1 <- sim(fx, c(200L, 500L), c("cstat", "mape"), 1)
  r2 <- sim(fx, c(200L, 500L, 800L), c("cstat", "mape", "cal_slope"), 2,
            method_args = list(shrinkage = TRUE))
  r2$method <- "logistic_shrunk"

  comp <- pmgen_compare(r1, r2, stat = "mean")

  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 4L)   # 2 shared N x 2 shared metrics
  expect_named(comp, c("N", "metric", "logistic", "logistic_shrunk",
                       "difference"))
  expect_setequal(comp$N, c(200L, 500L))

  s1 <- r1$summaries[["500"]]
  s2 <- r2$summaries[["500"]]
  row <- comp[comp$N == 500L & comp$metric == "mape", ]
  expect_equal(row$logistic, s1$mean[s1$metric == "mape"])
  expect_equal(row$logistic_shrunk, s2$mean[s2$metric == "mape"])
  expect_equal(row$difference, row$logistic - row$logistic_shrunk)
})
