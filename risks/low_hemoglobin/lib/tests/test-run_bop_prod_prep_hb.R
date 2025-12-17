library(testthat)

source("risk_hemoglobin/lib/run_bop_prod_prep_hb.R")

#' Test add_indicator_cols()
test_that("add_indicator_cols adds correct columns and values", {
  temp_data <- data.table::data.table(id = 1:3, value = c(10, 20, 30))
  temp_path <- "/data_elbw_trim1_cov_rev_01.csv"
  
  result <- add_indicator_cols(temp_data, temp_path)
  
  # Check that new columns are added
  expect_true("out_sub" %in% colnames(result))
  expect_true("trim1" %in% colnames(result))
  expect_true("trim2" %in% colnames(result))
  expect_true("trim3" %in% colnames(result))
  expect_true("cov_rev_agnostic" %in% colnames(result))
  expect_true("cov_con_socio_0" %in% colnames(result))
  expect_true("incl_undef_ptb" %in% colnames(result))
  expect_true("incl_undef_lbw" %in% colnames(result))
  expect_true("incl_undef_sga" %in% colnames(result))
  expect_true("overall" %in% colnames(result))
  
  # Check relevant column values based on data path
  expect_equal(result$out_sub[1], "elbw")
  expect_equal(result$trim1[1], 1)
  expect_equal(result$trim2[1], 0)
  expect_equal(result$trim3[1], 0)
  expect_equal(result$cov_rev_agnostic[1], 1)
  expect_equal(result$cov_con_socio_0[1], 0)
  expect_equal(result$incl_undef_ptb[1], 0)
  expect_equal(result$incl_undef_lbw[1], 0)
  expect_equal(result$incl_undef_sga[1], 0)
  expect_equal(result$overall[1], 0)
})


#' Test run_bop_prod_prep_hb()
test_that("run_bop_prod_prep_hb processes valid paths as expected", {
  temp_dir <- tempdir()
  
  # Create test data file
  temp_data <- data.table::data.table(id = 1:3, value = c(10, 20, 30))
  temp_path <- file.path(temp_dir, "PATH TO DATASET")
  if (!dir.exists(dirname(temp_path))) {
    dir.create(dirname(temp_path), recursive = TRUE, mode = "775")
  }
  data.table::fwrite(temp_data, temp_path)
  
  # Create test map of valid data file path
  map_data <- data.frame(filepath = temp_path)
  map_path <- file.path(temp_dir, "Maternal_summarytable.csv")
  data.table::fwrite(map_data, map_path)
  
  # Check that valid paths are processed
  expected_path <- "PATH TO DATASET"
  run_bop_prod_prep_hb(paths = map_path)
  
  expect_true(file.exists(expected_path))
  
  # Create test map of invalid data file path
  map_data <- data.frame(filepath = "invalid/path.csv")
  map_path <- file.path(temp_dir, "Maternal_summarytable.csv")
  data.table::fwrite(map_data, map_path)
  
  # Check that invalid paths return a warning
  expect_warning(run_bop_prod_prep_hb(paths = map_path))

  file.remove(expected_path)
})
