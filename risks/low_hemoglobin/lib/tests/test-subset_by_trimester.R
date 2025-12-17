source("risk_hemoglobin/lib/subset_by_trimester.R")

dat <- data.table::data.table(
  values = 1:3,
  exp_timepoint = c("first_trimester", "second_trimester", "third_trimester")
)

testthat::test_that("Correct message and data returned when no trimester is specified", {
  testthat::expect_message(subset_by_trimester(dat, NA), "No trimester specified, using all relevant observations")
})

testthat::test_that("Data is correctly subset for the first trimester", {
  expected <- dat[dat$exp_timepoint == "first_trimester", ]
  result <- subset_by_trimester(dat, 1)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Data is correctly subset for the second trimester", {
  expected <- dat[dat$exp_timepoint == "second_trimester", ]
  result <- subset_by_trimester(dat, 2)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Data is correctly subset for the third trimester", {
  expected <- dat[dat$exp_timepoint == "third_trimester", ]
  result <- subset_by_trimester(dat, 3)
  testthat::expect_equal(result, expected)
})

testthat::test_that("Invalid trimester value results in error", {
  testthat::expect_error(subset_by_trimester(dat, 4), "Invalid value provided for `trimester`")
})
