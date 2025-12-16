source("risk_hemoglobin/lib/drop_cols.R")

dat <- data.table::data.table(
  col_1 = 1:5,
  col_2 = 1:5,
  col_3 = 1:5
)

testthat::test_that("Drops a single column", {
  result <- drop_cols(dat, covs_to_drop = "col_1")
  testthat::expect_false("col_1" %in% names(result))
})

testthat::test_that("Drops multiple columns with comma separator", {
  result <- drop_cols(dat, covs_to_drop = "col_1, col_2")
  testthat::expect_false(all(c("col_1", "col_2") %in% names(result)))
})

testthat::test_that("Drops multiple columns with space separator", {
  result <- drop_cols(dat, covs_to_drop = "col_1 col_2")
  testthat::expect_false(all(c("col_1", "col_2") %in% names(result)))
})

testthat::test_that("Drops columns with mixed comma and space separators", {
  result <- drop_cols(dat, covs_to_drop = "col_1, col_2 col_3")
  testthat::expect_false(all(c("col_1", "col_2", "col_3") %in% names(result)))
})

testthat::test_that("Returns message when passed non-existent column names", {
  testthat::expect_message(drop_cols(dat, covs_to_drop = "col_4, col_5"),
                           "Columns not found and not dropped: col_4, col_5")
})

testthat::test_that("Throws error when covs_to_drop is NULL", {
  testthat::expect_error(drop_cols(dat, covs_to_drop = NULL))
})
