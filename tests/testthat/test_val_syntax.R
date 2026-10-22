# general checks
#======================================
test_that("Test 1.1.1", {
  errormessage <- 'argument "type" is missing, with no default'
  expect_error(pmvalsampsize(prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), oeciwidth = 1), errormessage)
})

test_that("Test 1.1.2", {
  errormessage <- 'prevalence must be specified for binary sample size'
  expect_error(pmvalsampsize(type = "b", cstatistic = 0.8, lpnormal = c(-5,2.5)), errormessage)
})

test_that("Test 1.1.3", {
  errormessage <- 'cstatistic must be specified for binary outcome models'
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpnormal = c(-5,2.5)), errormessage)
})


# checks for binary sample size
#======================================
test_that("Test 1.3.1", {
  errormessage <- "An LP distribution must be specified"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8), errormessage)
})

test_that("Test 1.3.2", {
  errormessage <- "Only one LP distribution option can be specified"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpnormal = c(-5,2.5), lpbeta = c(0.5,0.5), cstatistic = 0.8), errormessage)
})

test_that("Test 1.3.3", {
  errormessage <- "Only one LP distribution option can be specified"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpnormal = c(-5,2.5), lpcstat = TRUE, cstatistic = 0.8), errormessage)
})

test_that("Test 1.3.4", {
  errormessage <- "lpcstat must be TRUE or FALSE"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpcstat = "yes", cstatistic = 0.8), errormessage)
})

test_that("Test 1.3.5", {
  errormessage <- "cstatistic must be greater than 0.5 and less than 1 when lpcstat = TRUE"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpcstat = TRUE, cstatistic = 0.5), errormessage)
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpcstat = TRUE, cstatistic = 0.4), errormessage)
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpcstat = TRUE, cstatistic = 1), errormessage)
})


# retired arguments and legacy inputs (pmvalsampsize 0.1.0)
#======================================
test_that("Test 1.4.1", {
  for (arg in c("tolerance", "increment", "trace")) {
    args <- list(type = "b", prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), simobs = 1000)
    args[[arg]] <- if (arg == "trace") TRUE else 0.1
    expect_warning(suppressMessages(do.call(pmvalsampsize, args)),
                   paste(arg, "was retired in pmsampsize 2.0.0 and is ignored"))
  }
})

test_that("Test 1.4.2", {
  expect_warning(
    test <- suppressMessages(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpcstat = -4.7, oeciwidth = 1, seed = 1234, simobs = 10000)),
    "a numeric lpcstat is no longer used and is treated as lpcstat = TRUE")
  expected <- suppressMessages(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpcstat = TRUE, oeciwidth = 1, seed = 1234, simobs = 10000))
  expect_identical(test, expected)
})

test_that("Test 1.4.3", {
  # lpcstat = NA (the 0.1.0 default) means the option is not used
  expect_no_warning(suppressMessages(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), lpcstat = NA, simobs = 1000)))
})




