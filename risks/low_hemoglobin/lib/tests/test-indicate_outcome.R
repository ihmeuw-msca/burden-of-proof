source("risk_hemoglobin/lib/indicate_outcome.R")

dat <- data.table::data.table(
  bio_type_1 = c(rep("Hemoglobin", 5), rep("hemoglobin", 5), rep("Hematocrit", 5), rep("Packed_Cell_Volume", 4), NA, NA, NA),
  bio_type_2 = c(rep(NA, 15), rep("%", 4), NA, NA, NA),
  out_type = c("Stillbirth", "Neonatal mortality", "Perinatal mortality", "Maternal mortality", "Pre-eclampsia", 
                        "Postpartum Hemorrhage - Primary", "Stillbirth", "Neonatal mortality", "Perinatal mortality", "Maternal mortality",
                        "Pre-eclampsia", "Postpartum Hemorrhage - Primary", "Stillbirth", "Neonatal mortality", "Perinatal mortality",
                        "Maternal mortality", "Pre-eclampsia", "Antenatal Hemorrhage", "Postpartum Hemorrhage - Primary", "Postpartum Hemorrhage - Unspecified", "Maternal sepsis", "Neonatal mortality"),
  out_timepoint = c("at delivery", "0-6days of life", "7-27 days of life", "0-27 days of life", "unspecified", NA, 
                    "at delivery", "0-6days of life", "7-27 days of life", "0-27 days of life", "unspecified", NA,
                    "at delivery", "0-6days of life", "7-27 days of life", "0-27 days of life", "unspecified", NA,
                    "at delivery", "0-6days of life", "7-27 days of life", "0-27 days of life"),
  out_unit = c("percentile", rep("grams", 7), rep("weeks", 7), rep("percentile", 7)),
  out_lower = c(90, 100, 36.9, 4500, 10, 100, 20, 32, 31.9, 1000, 28, 32, 27.9, 28, 2500, 32, NA, 100, 27, 32, 1500, 50),
  out_upper = c(1500, 1000, 2500, 7, 28, 37, 32, 1500, 36.9, 2500, 10, 28, 28, 32, 1500, 5, 10, 32, 28, 37, 100, 28),
  prenataldep = c(rep(0, 10), rep(1, 10), NA, NA),
  postnataldep = c(rep(1, 10), rep(0, 10), NA, NA),
  analysis_exclcrx = c(0, rep(1, 20), NA)
)


#' Ensure the expected handling of missing outcome_severity for LBW
testthat::test_that("Missing outcome_severity for LBW returns overall LBW", {
  result <- indicate_outcome(dat, outcome = "lbw", outcome_severity = "")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(result$out_upper >= 2499 &
                          result$out_upper <= 2500)
})


#' Ensure returned output is a data table
testthat::test_that("Parent function indicate_outcome() returns a data table", {
  result <- indicate_outcome(dat, outcome = "lbw")
  testthat::expect_s3_class(result, "data.table")
})


# Test LBW outcome indication ----

# Test to ensure that macrosomia observations are removed
testthat::test_that("Macrosomia observations are removed for LBW outcomes", {
  result <- indicate_lbw(dat, outcome = "lbw")
  testthat::expect_true(all(result$macrosomia_ind == 0))
  testthat::expect_false(any(result$out_lower >= 4000))
})

# Test to ensure LBW observations are correctly identified
testthat::test_that("LBW observations are correctly identified", {
  result <- indicate_lbw(dat, outcome = "lbw")
  testthat::expect_true(all(result$lbw_ind == 1))
  testthat::expect_false(any(result$out_upper > 2500, na.rm = TRUE))
})

# Test to ensure VLBW observations are correctly identified
testthat::test_that("VLBW observations are correctly identified when specified", {
  result <- indicate_lbw(dat, outcome = "lbw", outcome_severity = "vlbw")
  testthat::expect_true(all(result$vlbw_ind == 1))
  testthat::expect_false(any(result$out_upper > 1500, na.rm = TRUE))
})

# Test to ensure ELBW observations are correctly identified
testthat::test_that("ELBW observations are correctly identified when specified", {
  result <- indicate_lbw(dat, outcome = "lbw", outcome_severity = "elbw")
  testthat::expect_true(all(result$elbw_ind == 1))
  testthat::expect_false(any(result$out_upper > 1000, na.rm = TRUE))
})

# Test SGA/LGA outcome indication ----

# Ensure SGA observations are output as expected
testthat::test_that("SGA outcomes are output as expected", {
  result <- indicate_sga_lga(dat, outcome = "sga")
  testthat::expect_true(all(result$sga_ind == 1 &
                              result$out_upper >= 9.9 &
                              result$out_upper <= 10))
})


# Ensure that LGA observations are output as expected
testthat::test_that("LGA outcomes are output as expected", {
  result <- indicate_sga_lga(dat, outcome = "lga")
  testthat::expect_true(all(result$lga_ind == 1 &
                              result$out_lower >= 90 &
                              result$out_lower <= 90.1))
})


# Test maternal hemorrhage outcome indication ----

#' Ensure that antenatal hemorrhage observations are output as expected
testthat::test_that("Antenatal hemorrhage outcomes are output as expected", {
  result <- indicate_maternal_outcome(dat, outcome = "mat_hem", outcome_severity = "antenatal_hem")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type == "Antenatal Hemorrhage"))
})

#' Ensure that postpartum hemorrhage observations are output as expected
testthat::test_that("Postpartum hemorrhage outcomes are output as expected", {
  result <- indicate_maternal_outcome(dat, outcome = "mat_hem")  # no severity
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type %in% c(
    "Postpartum Hemorrhage - Primary",
    "Postpartum Hemorrhage - Unspecified")))
})


# Test maternal mortality outcome indication ----

#' Ensure that maternal mortality observations are output as expected
testthat::test_that("Maternal mortality outcomes are output as expected", {
  result <- indicate_maternal_outcome(dat, outcome = "mat_mort")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type == "Maternal mortality"))
})

# Test maternal sepsis outcome indication ----

# Ensure that maternal sepsis observations are output as expected
testthat::test_that("Maternal sepsis outcomes are output as expected", {
  result <- indicate_maternal_outcome(dat, outcome = "mat_sepsis")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type == "Maternal sepsis"))
})


# Test neonatal mortality outcome indication ----

#' Ensure that neonatal mortality obervations are output as expected
testthat::test_that("Neonatal mortality outcomes are output as expected", {
  result <- indicate_neo_mortality(dat, outcome = "neo_mort")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type %in% c("total neonatal",
                                                   "early neonatal",
                                                   "late neonatal")))
})

# Test neonatal sepsis outcome indication ----

# Ensure neonatal sepsis observations are output as expected
testthat::test_that("Neonatal sepsis outcomes are output as expected", {
  result <- indicate_neo_sepsis(dat, outcome = "neo_sepsis")
  testthat::expect_true(data.table::is.data.table(result))
})


# Test PPD outcome indication ----

# Ensure PPD observations are output as expected
testthat::test_that("PPD outcomes are output as expected", {
  # Test for prenatal depression
  result <- indicate_maternal_outcome(dat, outcome = "ppd", outcome_severity = "prenatal_dep")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$prenataldep == 1))
  
  # Test for postnatal depression
  result <- indicate_maternal_outcome(dat, outcome = "ppd", outcome_severity = "postnatal_dep")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$postnataldep == 1))
})


# Test preeclampsia outcome indication ----

# Ensure preeclampsia observations are output as expected
testthat::test_that("Preeclampsia outcomes are output as expected", {
  dat_preecl <- dat[, outsub_htn := out_type]
  result <- indicate_maternal_outcome(dat_preecl, outcome = "preecl")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$outsub_htn == "Pre-eclampsia"))
})


# Test PTB outcome indication ----

# Ensure PTB observations are output based on severity
testthat::test_that("Moderately PTB outcomes are output as expected", {
  result <- indicate_outcome(dat, outcome = "ptb", outcome_severity = "mptb")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_unit == "weeks"))
  testthat::expect_true(all(result$out_upper >= 36.9 &
                              result$out_upper <= 37))
  testthat::expect_true(all(result$out_lower >= 31.9 &
                              result$out_lower <= 32))
})

testthat::test_that("Very PTB outcomes are output as expected", {
  result <- indicate_outcome(dat, outcome = "ptb", outcome_severity = "vptb")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_unit == "weeks"))
  testthat::expect_true(all(result$out_upper >= 31.9 &
                              result$out_upper <= 32))
  testthat::expect_true(all(result$out_lower >= 27.9 &
                              result$out_lower <= 28))
})

testthat::test_that("Extremely PTB outcomes are output as expected", {
  result <- indicate_outcome(dat, outcome = "ptb", outcome_severity = "eptb")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_unit == "weeks"))
  testthat::expect_true(all(result$out_upper >= 27.9 &
                              result$out_upper <= 28))
})

testthat::test_that("Overall PTB outcomes are output as expected", { 
  result <- indicate_outcome(dat, outcome = "ptb")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_unit == "weeks"))
  testthat::expect_true(all(result$out_upper >= 36.9 &
                              result$out_upper <= 37))
})


# Test stillbirth outcome indication ----

#' Ensure that stillbirth observations are output as expected
testthat::test_that("indicate_neonatal_outcome as expected subsets data", {
  result <- indicate_stillbirth(dat, outcome = "stillbirth")
  testthat::expect_true(data.table::is.data.table(result))
  testthat::expect_true(all(result$out_type_original == "Stillbirth" |
                              (result$out_type_original %in% c("Neonatal mortality",
                                                               "Perinatal mortality") &
                                 result$out_timepoint == "at delivery") |
                              result$out_lower == 20))
})
