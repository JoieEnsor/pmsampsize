test_that("Test 4.1", {
  ss <- c(905, 4501, 4252, 4501)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), oeciwidth = 1, simobs = 1000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})

test_that("Test 4.2", {
  ss <- c(905, 1641, 4252, 4252)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpbeta = c(0.5,0.5), oeciwidth = 1, simobs = 1000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})

# pmvalsampsize 0.1.0 (lpcstat = -4.7 as a start value, used as is) gave
# c(905, 16175, 4252, 16175); 2.0.0 calculates the non-events mean exactly
test_that("Test 4.3", {
  ss <- c(905, 16270, 4252, 16270)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpcstat = TRUE, oeciwidth = 1, seed = 1234, simobs = 10000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})

test_that("lpcstat non-events mean is solved exactly", {
  expect_equal(.lpcstat_m2(prevalence = 0.018, cstatistic = 0.8), -4.7075, tolerance = 1e-4)
  expect_equal(.lpcstat_m2(prevalence = 0.2, cstatistic = 0.75), -1.8412, tolerance = 1e-4)
  expect_equal(.lpcstat_m2(prevalence = 0.5, cstatistic = 0.9), -1.6424, tolerance = 1e-4)
})

test_that("lpcstat returns for prevalence = 0.2 (looped forever in 0.1.0)", {
  test <- suppressMessages(pmvalsampsize(type = "b", prevalence = 0.2, cstatistic = 0.75, lpcstat = TRUE, simobs = 10000))
  expect_s3_class(test, "pmvalsampsize")
  expect_true(test$simulated_data_cstat > 0.7 & test$simulated_data_cstat < 0.8)
})

test_that("lpcstat gives one message with the non-events mean and prevalence", {
  expect_message(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpcstat = TRUE, simobs = 1000),
                 "mean in non-event group=-4.70755.*target prevalence=0.018")
})


