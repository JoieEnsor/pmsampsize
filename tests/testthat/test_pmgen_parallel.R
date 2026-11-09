# pmgen_simulate() runs its simulations with furrr::future_map() when furrr is
# installed and with lapply() otherwise. Each simulation sets its own seed and
# the session's RNG kinds, so all routes must give identical results.

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

# Workers load the *installed* pmsampsize, which may predate pmgen (e.g.
# under devtools::test() with an older version installed)
skip_if_workers_lack_pmgen <- function() {
  has_pmgen <- future::value(future::future(
    requireNamespace("pmsampsize", quietly = TRUE) &&
      exists("pmgen_simulate", envir = asNamespace("pmsampsize")) &&
      exists(".set_seed_kinds", envir = asNamespace("pmsampsize"))
  ))
  skip_if_not(has_pmgen,
              "installed pmsampsize lacks pmgen or its current seeding")
}

test_that("lapply() fallback and furrr (sequential plan) give identical results", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  fx <- make_sim_fixtures()

  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  res_furrr <- run_sim(fx)
  rng_furrr <- .Random.seed
  local_mocked_bindings(.use_furrr = function() FALSE)
  res_lapply <- run_sim(fx)
  rng_lapply <- .Random.seed

  expect_same_results(res_furrr, res_lapply)
  # Both routes leave the session's stream in the same place
  expect_identical(rng_furrr, rng_lapply)
})

test_that("parallel (multisession) and sequential runs give identical results", {
  skip_on_cran()
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  fx <- make_sim_fixtures()

  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  skip_if_workers_lack_pmgen()

  res_par <- run_sim(fx)
  future::plan(future::sequential)
  res_seq <- run_sim(fx)

  expect_same_results(res_par, res_seq)
})

test_that("multisession workers use the session's non-default RNG kinds", {
  skip_on_cran()
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  fx <- make_sim_fixtures()

  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  skip_if_workers_lack_pmgen()

  restore_rng <- set_rng_kinds("L'Ecuyer-CMRG", "Box-Muller", "Rounding")
  on.exit(restore_rng(), add = TRUE)
  res_par <- run_sim(fx)
  future::plan(future::sequential)
  res_seq <- run_sim(fx)

  expect_same_results(res_par, res_seq)
  expect_identical(RNGkind(), c("L'Ecuyer-CMRG", "Box-Muller", "Rounding"))
})
