test_that("pmgen_reference() recovers target prevalence and c-stat from distributional assumptions", {
  # Minimal reproducible test: independent standard-normal predictors,
  # standard Bernoulli predictors. Tolerances are loose because the bisection
  # is run on a modest n_bisect_pop for speed.
  ref <- pmgen_reference(
    joint_data = list(
      x1 = list(dist = "normal",    mean = 0, sd = 1),
      x2 = list(dist = "normal",    mean = 0, sd = 1),
      x3 = list(dist = "bernoulli", prob  = 0.5)
    ),
    prevalence  = 0.20,
    cstatistic  = 0.72,
    tolerance   = 0.005,
    n_synth     = 50000,
    max_iter    = 60,
    seed        = 42
  )

  expect_s3_class(ref, "pmgen_ref_model")
  expect_equal(ref$prevalence_achieved, 0.20, tolerance = 0.02)
  expect_equal(ref$cstat_achieved,      0.72, tolerance = 0.02)
  expect_named(ref$coef, c("x1", "x2", "x3"))
  expect_true(is.matrix(ref$unit_info_matrix))
  expect_equal(dim(ref$unit_info_matrix), c(4L, 4L))  # 3 predictors + intercept
  expect_true(is.matrix(ref$inv_unit_info))
})

test_that("pmgen_reference() recovers pre-eclampsia-like parameters", {
  # 10-predictor scenario matching the paper (Section 3).
  # Predictors are 5 standard normal + 5 Bernoulli(0.5).
  # Target: prevalence = 0.68, c-stat = 0.76.
  # We use equal relative weights (coef = NULL) as a conservative default.
  skip_on_cran()  # slow: 1e5 synth population

  ref <- pmgen_reference(
    joint_data = list(
      x1  = list(dist = "normal",    mean = 0, sd = 1),
      x2  = list(dist = "normal",    mean = 0, sd = 1),
      x3  = list(dist = "normal",    mean = 0, sd = 1),
      x4  = list(dist = "normal",    mean = 0, sd = 1),
      x5  = list(dist = "normal",    mean = 0, sd = 1),
      x6  = list(dist = "bernoulli", prob  = 0.5),
      x7  = list(dist = "bernoulli", prob  = 0.5),
      x8  = list(dist = "bernoulli", prob  = 0.5),
      x9  = list(dist = "bernoulli", prob  = 0.5),
      x10 = list(dist = "bernoulli", prob  = 0.5)
    ),
    prevalence  = 0.68,
    cstatistic  = 0.76,
    tolerance   = 0.005,
    n_synth     = 100000,
    max_iter    = 80,
    seed        = 42
  )

  expect_s3_class(ref, "pmgen_ref_model")
  expect_equal(ref$prevalence_achieved, 0.68, tolerance = 0.02)
  expect_equal(ref$cstat_achieved,      0.76, tolerance = 0.02)
  expect_equal(dim(ref$unit_info_matrix), c(11L, 11L))
})

test_that("pmgen_reference() works from a supplied data frame", {
  set.seed(1)
  n <- 5000
  df <- data.frame(
    age = rnorm(n, 50, 10),
    bp  = rnorm(n,  0,  1),
    sex = rbinom(n, 1, 0.5)
  )

  ref <- pmgen_reference(
    joint_data  = df,
    predictors  = c("age", "bp", "sex"),
    prevalence  = 0.30,
    cstatistic  = 0.70,
    tolerance   = 0.01,
    n_synth     = 5000,
    max_iter    = 80,
    seed        = 7
  )

  expect_s3_class(ref, "pmgen_ref_model")
  expect_equal(ref$prevalence_achieved, 0.30, tolerance = 0.03)
  expect_equal(ref$cstat_achieved,      0.70, tolerance = 0.03)
})

test_that("pmgen_reference() errors on bad inputs", {
  expect_error(
    pmgen_reference(
      joint_data = list(x1 = list(dist = "normal", mean = 0, sd = 1)),
      prevalence = 1.5,
      cstatistic = 0.7
    ),
    "prevalence"
  )
  expect_error(
    pmgen_reference(
      joint_data = list(x1 = list(dist = "normal", mean = 0, sd = 1)),
      prevalence = 0.3,
      cstatistic = 0.4  # below 0.5
    ),
    "cstatistic"
  )
})

test_that("unit information matrix is symmetric and positive definite", {
  ref <- pmgen_reference(
    joint_data = list(
      x1 = list(dist = "normal",    mean = 0, sd = 1),
      x2 = list(dist = "bernoulli", prob  = 0.5)
    ),
    prevalence  = 0.25,
    cstatistic  = 0.65,
    tolerance   = 0.01,
    n_synth     = 20000,
    max_iter    = 60,
    seed        = 99
  )

  I <- ref$unit_info_matrix
  expect_true(isSymmetric(I, tol = 1e-10))
  eigvals <- eigen(I, only.values = TRUE)$values
  expect_true(all(eigvals > 0))
})
