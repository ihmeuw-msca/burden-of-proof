source("risk_hemoglobin/lib/invert_and_swap_risks.R")

dat <- data.table::data.table(
  alt_risk_upper = c(140, 120, 140),
  alt_risk_lower = c(110, 80, 120),
  ref_risk_upper = c(130, 110, 130),
  ref_risk_lower = c(100, 70, 110),
  hb_risk_type = c("low_hb", "high_hb", "low_hb"),
  ln_rr = c(1.0, -1.5, 2.0)
)

testthat::test_that("ln_rr values are correctly inverted", {
  expected_ln_rr <- c(-1.0, 1.5, -2.0)
  result <- invert_and_swap_risks(dat)
  testthat::expect_equal(result$ln_rr, expected_ln_rr)
})

testthat::test_that("Risk bounds are correctly swapped", {
  expected_alt_upper <- c(130, 110, 130)
  expected_alt_lower <- c(100, 70, 110)
  
  result <- invert_and_swap_risks(dat)
  testthat::expect_equal(result$alt_risk_upper, expected_alt_upper)
  testthat::expect_equal(result$alt_risk_lower, expected_alt_lower)
})

testthat::test_that("hb_risk_type values are correctly updated", {
  expected_hb_risk_type <- c("high_hb", "low_hb", "high_hb")
  
  result <- invert_and_swap_risks(dat)
  testthat::expect_equal(result$hb_risk_type, expected_hb_risk_type)
})
