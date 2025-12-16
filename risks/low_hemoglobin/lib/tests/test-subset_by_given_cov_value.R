source("risk_hemoglobin/lib/subset_by_given_cov_value.R")

dat <- data.table::data.table(
  values = 1:4,
  cov_rev = c(0, 1, 0, 1),
  cov_con_socio = c(1, 1, 0, 0)
)

#' Ensure that subsetting by cov_con_socio works as expected
testthat::test_that("Subsetting to 'cov_con_socio = 0' works as expected", {
  expected <- dat[dat$cov_con_socio == 0, ]
  result <- subset_by_given_cov_value(dat, cov_con_socio = 0)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Subsetting to 'cov_con_socio = 1' works as expected", {
  expected <- dat[dat$cov_con_socio == 1, ]
  result <- subset_by_given_cov_value(dat, cov_con_socio = 1)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Correct message is printed for 'cov_con_socio = 0'", {
  testthat::expect_message(subset_by_given_cov_value(dat, cov_con_socio = 0),
                           "Subsetting by cov_con_socio = 0")
})

testthat::test_that("Correct message is printed for 'cov_con_socio = 1'", {
  testthat::expect_message(subset_by_given_cov_value(dat, cov_con_socio = 1),
                           "Subsetting by cov_con_socio = 1")
})


#' Ensure that subsetting by cov_rev works as expected
testthat::test_that("Subsetting to 'cov_rev = 0' works as expected", {
  expected <- dat[dat$cov_rev == 0, ]
  result <- subset_by_given_cov_value(dat, cov_rev = 0)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Subsetting to 'cov_rev = 1' works as expected", {
  expected <- dat[dat$cov_rev == 1, ]
  result <- subset_by_given_cov_value(dat, cov_rev = 1)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Correct message is printed for cov_rev = 0", {
  testthat::expect_message(subset_by_given_cov_value(dat, cov_rev = 0),
                           "Subsetting by cov_rev = 0")
})

testthat::test_that("Correct message is printed for 'cov_rev = 1'", {
  testthat::expect_message(subset_by_given_cov_value(dat, cov_rev = 1),
                           "Subsetting by cov_rev = 1")
})


#' Ensure that not subsetting works as expected
testthat::test_that("Passing '0,1' to either cov returns original data", {
  result <- subset_by_given_cov_value(dat, cov_con_socio = "0,1")
  testthat::expect_equal(result, dat)
  
  result <- subset_by_given_cov_value(dat, cov_rev = "0,1")
  testthat::expect_equal(result, dat)
})

testthat::test_that("Passing NULL or an empty string returns original data", {
  # For cov_con_socio
  result <- subset_by_given_cov_value(dat, cov_con_socio = NULL)
  testthat::expect_equal(result, dat)
  
  result <- subset_by_given_cov_value(dat, cov_con_socio = "")
  testthat::expect_equal(result, dat)
  
  # For cov_rev
  result <- subset_by_given_cov_value(dat, cov_rev = NULL)
  testthat::expect_equal(result, dat)
  
  result <- subset_by_given_cov_value(dat, cov_rev = "")
  testthat::expect_equal(result, dat)
  
})

testthat::test_that("Correct message prints when no subsetting is indicated", {
  testthat::expect_message(subset_by_given_cov_value(dat),
                           "No subsetting by cov_con_socio indicated")
  testthat::expect_message(subset_by_given_cov_value(dat),
                           "No subsetting by cov_rev indicated")
})
