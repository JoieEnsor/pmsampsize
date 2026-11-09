
# ---- Shared lightweight fixtures for Phase 2 tests ---------------------------
# These are created once per test session (not per test) by storing them in the
# helper environment. Tests that need them call make_sim_fixtures().

.sim_fixtures <- NULL

make_sim_fixtures <- function() {
  if (!is.null(.sim_fixtures)) return(.sim_fixtures)

  set.seed(1)
  n <- 10000
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rbinom(n, 1, 0.4)
  )

  ref <- pmgen_reference(
    joint_data = df,
    predictors = c("x1", "x2", "x3"),
    prevalence = 0.30,
    cstatistic = 0.70,
    tolerance  = 0.01,
    max_iter   = 60,
    seed       = 1
  )

  # n = 2000 leaves 8000 rows as the development pool
  tp <- pmgen_target_population(ref, n = 2000L, seed = 1)

  .sim_fixtures <<- list(ref = ref, tp = tp, df = df)
  .sim_fixtures
}

# Switch the session to other RNG kinds. Returns a function that puts the
# session's RNG state (stream and kinds) back as it was; use it with on.exit()
set_rng_kinds <- function(kind, normal.kind, sample.kind) {
  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  suppressWarnings(RNGkind(kind, normal.kind, sample.kind))  # "Rounding" warns
  function() {
    if (is.null(old_seed)) {
      suppressWarnings(RNGkind("default", "default", "default"))
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }
}
