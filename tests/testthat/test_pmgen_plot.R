# plot() needs ggplot2 and scales (Suggests)
skip_if_not_installed("ggplot2")
skip_if_not_installed("scales")

test_that("plot degradation returns a ggplot object", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 300L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("cstat", "cal_slope", "mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  p <- plot(result, type = "degradation")
  expect_s3_class(p, "ggplot")
})

test_that("plot assurance returns a ggplot object", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 300L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("cstat", "cal_slope", "mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  p <- plot(result, type = "assurance")
  expect_s3_class(p, "ggplot")
})

test_that("instability plot errors without store_predictions", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = 200L,
    n_sims  = 30L,
    method  = "logistic",
    metrics = c("cstat"),
    seed    = 1L,
    verbose = FALSE
  ))

  expect_error(plot(result, type = "instability"),
               "store_predictions = TRUE")
})

test_that("calibration plot errors without store_predictions", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = 200L,
    n_sims  = 30L,
    method  = "logistic",
    metrics = c("cstat"),
    seed    = 1L,
    verbose = FALSE
  ))

  expect_error(plot(result, type = "calibration"),
               "store_predictions = TRUE")
})

test_that("instability plot returns ggplot when store_predictions = TRUE", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N                 = 200L,
    n_sims            = 50L,
    method            = "logistic",
    metrics           = c("cstat"),
    store_predictions = TRUE,
    seed              = 1L,
    verbose           = FALSE
  ))

  expect_false(is.null(result$pred_matrices))
  expect_false(is.null(result$p_true))
  expect_equal(ncol(result$pred_matrices[["200"]]), 50L)
  expect_equal(length(result$p_true), nrow(fix$tp$data))

  p <- suppressMessages(plot(result, type = "instability", N = 200L))
  expect_s3_class(p, "ggplot")
})

test_that("calibration plot returns ggplot when store_predictions = TRUE", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N                 = 200L,
    n_sims            = 50L,
    method            = "logistic",
    metrics           = c("cstat"),
    store_predictions = TRUE,
    seed              = 1L,
    verbose           = FALSE
  ))

  p <- suppressMessages(plot(result, type = "calibration", N = 200L))
  expect_s3_class(p, "ggplot")
})

test_that(".resolve_n defaults to first N and informs", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N                 = c(100L, 300L),
    n_sims            = 30L,
    method            = "logistic",
    metrics           = c("cstat"),
    store_predictions = TRUE,
    seed              = 1L,
    verbose           = FALSE
  ))

  expect_message(
    plot(result, type = "instability"),
    "N = 100"
  )
})

test_that("plot.pmgen_result errors on unknown type", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = 100L,
    n_sims  = 20L,
    method  = "logistic",
    metrics = c("cstat"),
    seed    = 1L,
    verbose = FALSE
  ))

  expect_error(plot(result, type = "bananas"))
})
