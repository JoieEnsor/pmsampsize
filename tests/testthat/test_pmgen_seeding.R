# pmgen_simulate() seeds each N's simulations from `seed`, the N value and
# n_sims only (see .pmgen_sim_seeds()).

sim_seeded <- function(fx, N, n_sims = 12L, seed = 2024) {
  suppressMessages(pmgen_simulate(
    fx$ref, fx$tp, N = N, n_sims = n_sims, method = "logistic",
    metrics = c("cstat", "cal_slope", "mape"), seed = seed, verbose = FALSE
  ))
}

test_that("an N value's results do not depend on the other N values", {
  fx <- make_sim_fixtures()
  alone  <- sim_seeded(fx, N = 200L)
  after  <- sim_seeded(fx, N = c(60L, 200L))
  before <- sim_seeded(fx, N = c(200L, 120L, 80L))

  expect_identical(after$posterior[["200"]], alone$posterior[["200"]])
  expect_identical(before$posterior[["200"]], alone$posterior[["200"]])
  expect_identical(after$summaries[["200"]], alone$summaries[["200"]])
})

test_that("different N values get different simulation seeds", {
  seeds <- lapply(c(60L, 61L, 200L), function(n) .pmgen_sim_seeds(123L, n, 50L))
  expect_false(any(seeds[[1]] %in% seeds[[2]]))
  expect_false(any(seeds[[1]] %in% seeds[[3]]))
})

test_that("increasing n_sims keeps the earlier simulations", {
  fx <- make_sim_fixtures()
  short <- sim_seeded(fx, N = 60L, n_sims = 6L)
  long  <- sim_seeded(fx, N = 60L, n_sims = 12L)

  expect_identical(unname(as.matrix(long$posterior[["60"]][1:6, ])),
                   unname(as.matrix(short$posterior[["60"]])))
  expect_identical(.pmgen_sim_seeds(99L, 60L, 1000L)[1:500],
                   .pmgen_sim_seeds(99L, 60L, 500L))
})

test_that("seed = NULL uses the session's stream", {
  fx <- make_sim_fixtures()
  set.seed(2024)
  from_session <- sim_seeded(fx, N = 60L, seed = NULL)
  from_arg     <- sim_seeded(fx, N = 60L, seed = 2024)
  expect_identical(from_session$posterior, from_arg$posterior)
})

test_that("the session's RNG stream moves only by set.seed(seed) and one draw", {
  fx <- make_sim_fixtures()
  sim_seeded(fx, N = c(60L, 200L))
  after_sim <- .Random.seed
  set.seed(2024)
  sample.int(.Machine$integer.max, 1L)
  expect_identical(after_sim, .Random.seed)
})

test_that("the session's RNG kinds are used and left unchanged", {
  fx <- make_sim_fixtures()
  default_kinds <- sim_seeded(fx, N = 60L)

  restore_rng <- set_rng_kinds("L'Ecuyer-CMRG", "Box-Muller", "Rounding")
  on.exit(restore_rng(), add = TRUE)
  other_kinds <- sim_seeded(fx, N = 60L)

  expect_false(identical(other_kinds$posterior, default_kinds$posterior))
  expect_identical(RNGkind(), c("L'Ecuyer-CMRG", "Box-Muller", "Rounding"))
  # Reproducible under the same kinds
  expect_identical(sim_seeded(fx, N = 60L)$posterior, other_kinds$posterior)
})
