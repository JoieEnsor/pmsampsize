# pmgen_simulate() runs its simulations with furrr::future_map() when furrr is
# installed and with lapply() otherwise. Each simulation sets its own seed,
# drawn before the loop, so all routes must give identical results.

run_sim <- function(fx) {
  suppressMessages(pmgen_simulate(
    fx$ref, fx$tp, N = c(60L, 200L), n_sims = 12L, method = "logistic",
    metrics = c("cstat", "cal_slope", "mape"), store_predictions = TRUE,
    seed = 2024, verbose = FALSE
  ))
}

expect_same_results <- function(a, b) {
  expect_identical(a$posterior, b$posterior)
  expect_identical(a$summaries, b$summaries)
  expect_identical(a$pred_matrices, b$pred_matrices)
  expect_identical(a$fit_info, b$fit_info)
}

test_that("lapply() fallback and furrr (sequential plan) give identical results", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  fx <- make_sim_fixtures()

  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  res_furrr <- run_sim(fx)
  local_mocked_bindings(.use_furrr = function() FALSE)
  res_lapply <- run_sim(fx)

  expect_same_results(res_furrr, res_lapply)
})

test_that("parallel (multisession) and sequential runs give identical results", {
  skip_on_cran()
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  fx <- make_sim_fixtures()

  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  # Workers load the *installed* pmsampsize, which may predate pmgen (e.g.
  # under devtools::test() with an older version installed)
  has_pmgen <- future::value(future::future(
    requireNamespace("pmsampsize", quietly = TRUE) &&
      exists("pmgen_simulate", envir = asNamespace("pmsampsize"))
  ))
  skip_if_not(has_pmgen, "installed pmsampsize has no pmgen functions")

  res_par <- run_sim(fx)
  future::plan(future::sequential)
  res_seq <- run_sim(fx)

  expect_same_results(res_par, res_seq)
})
