test_that("summary.pmgen_result returns a data frame with correct structure", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 200L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("cstat", "cal_slope", "mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  s <- summary(result)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("N", "metric", "mean", "sd", "q025", "q975",
                    "degradation_mean") %in% names(s)))
  # One row per (N, metric)
  expect_equal(nrow(s), length(result$N) * length(result$metrics))
  # Both N values present
  expect_setequal(unique(s$N), c(100L, 200L))
  # All three metrics present
  expect_setequal(unique(s$metric), c("cstat", "cal_slope", "mape"))
  # Mean cstat is a probability
  cstat_rows <- s[s$metric == "cstat", ]
  expect_true(all(cstat_rows$mean > 0 & cstat_rows$mean < 1))
  # Plain row names (quantile() names such as "2.5%" must not leak in)
  expect_identical(rownames(s), as.character(seq_len(nrow(s))))
  expect_identical(rownames(result$summaries[[1]]),
                   as.character(seq_len(nrow(result$summaries[[1]]))))
})

test_that("pmgen_min_n errors on missing pred_matrices or bad input", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 200L, 500L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("cstat", "cal_slope", "mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  expect_error(pmgen_min_n("not a result"),
               "pmgen_result")
  expect_error(pmgen_min_n(result, criteria = list()),
               "non-empty")
  expect_error(pmgen_min_n(result, criteria = list(list(max_mean = 0.1))),
               "named")
  expect_error(pmgen_min_n(result, criteria = list(nonexistent = list(max_mean = 0.1))),
               "not in result")
})

test_that("pmgen_min_n with max_mean criterion returns correct structure", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 300L, 1000L),
    n_sims  = 100L,
    method  = "logistic",
    metrics = c("cstat", "mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  mn <- pmgen_min_n(result, criteria = list(
    mape = list(max_mean = 0.99)   # very lenient — should always pass
  ))

  expect_s3_class(mn, "pmgen_min_n")
  expect_true(!is.na(mn$min_n))
  expect_s3_class(mn$table, "data.frame")
  expect_true("N" %in% names(mn$table))
  expect_true("all_met" %in% names(mn$table))
  expect_true("mape" %in% names(mn$table))
  # min_n must be one of the evaluated N values
  expect_true(mn$min_n %in% result$N)
})

test_that("pmgen_min_n returns NA when no N meets criteria", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(50L, 100L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("mape"),
    seed    = 1L,
    verbose = FALSE
  ))

  mn <- pmgen_min_n(result, criteria = list(
    mape = list(max_mean = 0.0)   # impossible threshold
  ))

  expect_true(is.na(mn$min_n))
  expect_false(any(mn$table$all_met))
})

test_that("pmgen_min_n range+min_prob criterion works", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(200L, 1000L),
    n_sims  = 100L,
    method  = "logistic",
    metrics = c("cal_slope"),
    seed    = 1L,
    verbose = FALSE
  ))

  # Very wide range and low probability — should pass at both N
  mn <- pmgen_min_n(result, criteria = list(
    cal_slope = list(range = c(-10, 10), min_prob = 0.01)
  ))
  expect_false(is.na(mn$min_n))
  expect_equal(mn$min_n, 200L)

  # Impossible range — should fail at all N
  mn2 <- pmgen_min_n(result, criteria = list(
    cal_slope = list(range = c(100, 200), min_prob = 0.99)
  ))
  expect_true(is.na(mn2$min_n))
})

test_that("print.pmgen_min_n runs without error", {
  fix <- make_sim_fixtures()
  result <- suppressMessages(pmgen_simulate(
    fix$ref, fix$tp,
    N       = c(100L, 500L),
    n_sims  = 50L,
    method  = "logistic",
    metrics = c("mape"),
    seed    = 1L,
    verbose = FALSE
  ))
  mn <- pmgen_min_n(result, criteria = list(mape = list(max_mean = 0.99)))
  expect_output(print(mn), "pmgen minimum sample size")
})
