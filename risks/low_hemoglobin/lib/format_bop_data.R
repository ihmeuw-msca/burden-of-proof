#' Low Hb function to format data for BoP
#' 
#' @param dat Data object output by outcome-specific cleaning function.
#' @param outcome String containing outcome of interest. Must be one of:
#'   "lbw", "lga", "mat_hem", "mat_mort", "mat_sepsis", "preecl", neo_mort",
#'   "neo_sepsis", "ppd", "ptb", "sga", "stillbirth".
#' @param clean_data_path Path to cleaned and formatted extractions file.
#' @param outcome_severity String indicating outcome severity. Must align with
#'   `outcome` values and be one of: "lbw", "vlbw", "elbw", "ptb", "mptb",
#'   "vptb", "eptb", "early_neonatal", "late_neonatal", "total_neonatal",
#'   "antenatal_hem", "prenatal_dep", "postnatal_dep".
#' @param trimester Number indicating trimester of pregnancy. Default `NULL`.
#' @cov_con_socio Desired value by which to subset `cov_con_socio`. Must be 0,
#'   1 or 0:1.
#' @cov_rev Desired value by which to subset `cov_rev`. Must be 0, 1 or 0:1.
#' 
#' @note The second time the file is written purposely overwrites an existing
#'   file to preserve file name, which facilitates running sensitivity analyses.

run_format_bop_data <- function(dat,
                                outcome,
                                clean_data_path,
                                outcome_severity = NULL,
                                trimester = NULL,
                                cov_con_socio = NULL,
                                cov_rev = NULL) {
  
  box::use(
    risk_hemoglobin/lib/params,
    risk_hemoglobin/lib/indicate_outcome[indicate_outcome],
    risk_hemoglobin/lib/format_bop_bundle[format_bop_bundle],
    risk_hemoglobin/lib/subset_by_trimester[subset_by_trimester],
    risk_hemoglobin/lib/invert_and_swap_risks[invert_and_swap_risks],
    risk_hemoglobin/lib/create_sensitivity_specific_string[
      create_sensitivity_specific_string
    ],
    data.table
  ) |> 
    suppressWarnings()
  
  cli::cli_progress_step("Beginning to format cleaned data for BoP...")
  
  # Subset to specified outcome or, if requested, outcome sensitivity subgroup
  dat <- indicate_outcome(dat = dat,
                          outcome = outcome,
                          outcome_severity = outcome_severity)
  
  # For PPD only, 
  if (outcome == "ppd") {
    dat <- format_bop_bundle(input_df = dat, release_id = 16)
  }
  
  # Impute missing hemoglobin reference and alternate values
  checkmate::assert_character(dat$bio_unit_1, min.chars = 1)
  checkmate::assert_true("conf_altitude_adjusted" %in% names(dat))
  cli::cli_progress_step("Imputing missing hemoglobin reference and alternate values...")
  dat_list <- bophf::impute_missing_hb_main(
    input_df = dat,
    biomarker_column_name = "bio_type_1",
    biomarker_unit_column_name = "bio_unit_1",
    ref_risk_upper_column_name = "ref_upper_1",
    ref_risk_lower_column_name = "ref_lower_1",
    alt_risk_upper_column_name = "alt_upper_1",
    alt_risk_lower_column_name = "alt_lower_1",
    apply_elevation_adjustment = TRUE,
    conf_elevation_adjust_column_name = "conf_altitude_adjusted"
  )
  
  # no_cutoff_df <- data.table::copy(df_list$no_cutoff_df)
  # high_hb_df <- data.table::copy(df_list$high_hb)
  # bad_unit_df <- data.table::copy(df_list$bad_conversion_df)
  # bad_bounds_df <- data.table::copy(df_list$bad_bounds_df)
  dat <- data.table::copy(dat_list$cleaned_df)
  
  # Cap alt_upper at 160 and drop observations where alt_lower > alt_upper
  dat$alt_risk_upper <- pmin(dat$alt_risk_upper, 160)
  dat <- dat[dat$alt_risk_lower <= dat$alt_risk_upper, ]
  
  # Transform mean and standard error columns
  cli::cli_progress_step("Transforming mean and standard error columns...")
  dat_list <- bophf::transform_mean_se_bop(
    input_df = dat,
    mean_column_name = "effect_estimate",
    upper_ui_column_name = "effect_upper",
    lower_ui_column_name = "effect_lower",
    se_column_name = "effect_se"
  )
  
  # Remove invalid rows
  #bad_transform_dat <- data.table::copy(dat_list$na_df)
  dat <- data.table::copy(dat_list$df)
  
  # Invert ln_rr, swap exposure bounds, and relabel risk type where applicable
  dat <- invert_and_swap_risks(dat)
  
  # Subset to given trimester, if specified
  if (!is.null(trimester) && !is.na(trimester)) {
    dat <- subset_by_trimester(dat = dat, trimester = trimester)
  }
  
  # Drop problematic bias covs
  cli::cli_progress_step("Dropping bias covariate columns with linear dependence...")
  if (nrow(dat) >= 4) {
    dat <- bophf::check_for_linear_dependent_covs(dat, drop = TRUE)
  } else if (nrow(dat) < 4) {
    dat <- dat |>
      dplyr::select(-starts_with("cov_"))
  }
  
  # Subset to overall or mutually-exclusive estimates based on exp_timepoint
  assertthat::assert_that(nrow(dat)> 0)
  if ((is.null(trimester) || is.na(trimester)) && outcome != "ppd") {
    cli::cli_progress_step("Identifying overall estimates...")
    dat <- bophf::identify_overall_estimates(dat)
  }
  
  # Inflate SE based on multiple exposure assessments
  dat <- dat |>
    dplyr::group_by(nid,
                    age_start,
                    age_end) |>
    dplyr::mutate(
      inflate_se = dplyr::case_when(
        dplyr::n_distinct(exp_timepoint) > 1 ~ dplyr::n_distinct(exp_timepoint),
        TRUE ~ 1
      ),
      ln_rr_se = ln_rr_se * sqrt(inflate_se)
    ) |>
    dplyr::ungroup() |>
    data.table::setDT()
  
  # Get final BoP data table
  cli::cli_progress_step("Preparing and writing out data file to use for BoP...",
                         msg_failed = "File writing unsuccessful",
                         msg_done = "Preparing and writing complete.")

  bop_dat <- bophf::get_final_bop_df(
    input_df = dat,
    pipeline_type = params$pipeline_type,
    other_cols_to_add = c("exp_timepoint", "hb_risk_type", "out_unit", "out_lower",
                          "out_upper", "out_type", "out_timepoint")
  ) |>
    data.table::setDT()
  
  # Write clean data snapshot file with informative suffix (for documentation)
  data.table::fwrite(
    bop_dat,
    clean_data_path
  )
  
  # Dynamically create string for outcome-sensitivity specific data directory
  dir_string <- create_sensitivity_specific_string(
    outcome = outcome,
    outcome_severity = outcome_severity,
    trimester = trimester,
    cov_con_socio = cov_con_socio,
    cov_rev = cov_rev
  )
  
  dir <- glue::glue("{params$risk_dir}/hb-{dir_string}/bop_model_sandbox/data")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    Sys.chmod(dir, "775", use_umask = FALSE)
  }
  
  # If needed, copy settings.yaml for a sensitivity subdir
  settings_path <- file.path(dir, "settings.yaml")
  if (!file.exists(settings_path)) {
    file.copy(glue::glue(
      "{params$risk_dir}/hb-{outcome}/bop_model_sandbox/data/settings.yaml"
    ),
    settings_path)
  }
  
  # Write same clean data file to appropiate BoP model sandbox
  rei_name <- nch::name_for("rei", params$rei_id)
  ro_pair <- paste(tolower(gsub(" ", "_", rei_name)), dat$acause[1], sep = "-")
  
  data.table::fwrite(
    bop_dat,
    glue::glue("{dir}/{ro_pair}.csv")
  )
  
  return(bop_dat)
}
