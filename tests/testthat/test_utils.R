test_that(".fast_auc matches the Mann-Whitney statistic", {
  set.seed(1)
  y <- rbinom(500, 1, 0.3)
  p <- plogis(rnorm(500, mean = y))
  w <- suppressWarnings(wilcox.test(p[y == 1], p[y == 0], exact = FALSE))$statistic
  expect_equal(.fast_auc(p, y), unname(w) / (sum(y == 1) * sum(y == 0)))
})

test_that(".fast_auc counts ties as 1/2", {
  expect_equal(.fast_auc(c(0.5, 0.5), c(1, 0)), 0.5)
  expect_equal(.fast_auc(c(0.2, 0.8, 0.8, 0.4), c(0, 1, 0, 1)), 0.625)
})

test_that(".fast_auc returns NA when an outcome group is empty", {
  expect_identical(.fast_auc(c(0.1, 0.9), c(1, 1)), NA_real_)
  expect_identical(.fast_auc(c(0.1, 0.9), c(0, 0)), NA_real_)
})
