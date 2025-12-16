#' Low Hb BoP data cleaning
#' 
#' @param outcome String containing outcome of interest. Must be one of:
#'   "lbw", "lga", "mat_hem", "mat_mort", "mat_sepsis", "preecl", neo_mort",
#'   "neo_sepsis", "ppd", "ptb", "sga", "stillbirth".
#' @param extracted_data_path Path to extractions data file.
#' @param outcome_severity String indicating outcome severity. Must align with
#'   `outcome` values and be one of: "lbw", "vlbw", "elbw", "ptb", "mptb",
#'   "vptb", "eptb", "early_neonatal", "late_neonatal", "antenatal_hem",
#'   "prenatal_dep", "postnatal_dep".
#' @param covs_to_drop Character string of columns to drop. Default `NULL`.

run_clean_ro_specific_data <- function(outcome,
                                       extracted_data_path,
                                       outcome_severity = NULL,
                                       covs_to_drop = NULL,
                                       cov_con_socio = NULL,
                                       cov_rev = NULL) {
  box::use(
    risk_hemoglobin/lib/params,
    risk_hemoglobin/lib/add_outcome_specific_cols[add_outcome_specific_cols],
    risk_hemoglobin/lib/drop_cols[drop_cols],
    risk_hemoglobin/lib/outcome_specific_cleaning[outcome_specific_cleaning],
    risk_hemoglobin/lib/subset_by_given_cov_value[subset_by_given_cov_value],
    data.table,
    transform_data/impute_age[impute_missing_ages]
  ) |>
    suppressWarnings()
  
  # Read in extraction data to be cleaned ---------------------------------
  
  cli::cli_progress_step("Reading in extracted data to be cleaned...")
  dat <- data.table::fread(extracted_data_path,
                           na.strings = c(NA_character_, ""))
  
  # REDCap cleaning, if needed --------------------------------------------
  
  # Coalesce on exposure units
  if (outcome %in% c("mat_mort", "mat_sepsis", "preecl", "ghtn", "hdop")) {
    unit_cols_to_coalesce <- dplyr::select(
      dat,
      c("bio_unit_1_hgb", "bio_unit_1_hct", "bio_unit_1_pcv")
    )
    unit_cols_to_coalesce <- bundleprep::coalesce_bundle_cols(
      unit_cols_to_coalesce,
      prefix = "bio_unit_1"
    )
    dat$bio_unit_1 <- unit_cols_to_coalesce$bio_unit_1
  }
  
  # Rename and dynamically create columns as needed to run BoP ------------
  if ("record_nid" %in% names(dat)) {
    dat[, nid := record_nid]
  }
  
  columns_to_rename <- grep("altitude|elevation", names(dat), value = TRUE)
  
  for (col in columns_to_rename) {
    suppressWarnings(data.table::setnames(dat,
                                          old = col,
                                          new = "conf_altitude_adjusted"))
  }
  
  old_new_names <- c("outcome_type" = "out_type",
                     "Standard.deviation" = "effect_se",
                     "lower" = "effect_lower",
                     "upper" = "effect_upper",
                     "mean" = "effect_estimate",
                     "Parameter.Measure" = "effect_measure",
                     "biomarker_type_1" = "bio_type_1",
                     "unit_of_measure_1" = "bio_unit_1",
                     "ref_upper_limit_b1" = "ref_upper_1",
                     "ref_lower_limit_b1" = "ref_lower_1",
                     "alt_upper_limit_1" = "alt_upper_1",
                     "alt_lower_limit_1" = "alt_lower_1",
                     "alt_risk_lower" = "alt_lower_1",
                     "alt_risk_upper" = "alt_upper_1",
                     "ref_risk_lower" = "ref_lower_1",
                     "ref_risk_upper" = "ref_upper_1",
                     "Outcome.timepoint" = "out_timepoint")
  
  lapply(names(old_new_names), function(old_name) {
    if (old_name %in% names(dat)) {
      data.table::setnames(
        dat,
        old = old_name,
        new = old_new_names[old_name]
      )
    }
  })
  
  req_cols <- c("out_type",
                "effect_se",
                "effect_lower",
                "effect_upper",
                "effect_estimate",
                "effect_measure",
                "bio_type_1",
                "bio_unit_1",
                "ref_upper_1",
                "ref_lower_1",
                "alt_upper_1",
                "alt_lower_1",
                "out_timepoint")
  
  checkmate::assert_names(
    names(dat),
    must.include = req_cols
  )
  
  if ("out_upper" %in% names(dat)) {
    dat$out_upper <- as.numeric(dat$out_upper)
  }
  
  # Perform outcome-specific one-off cleaning
  dat <- dat[!(is.na(bio_unit_1))]
  dat <- outcome_specific_cleaning(dat, outcome)
  
  # Drop columns, if requested
  if (!is.null(covs_to_drop) && !is.na(covs_to_drop)) {
    dat <- drop_cols(dat, covs_to_drop)
  }
  
  # Perform subsetting ----------------------------------------------------
  
  cli::cli_progress_step("Cleaning and subsetting data...")
  
  # Drop observations w/ missing effect_lower or effect_upper
  dat <- dat[!is.na(effect_lower) | !is.na(effect_upper), ]
  
  # Subset by both_anemia_overall_and_levels and, if present, included_in_model
  if ("both_anemia_overall_and_levels" %in% colnames(dat)) {
    dat <- dat[dat$both_anemia_overall_and_levels == 0 |
                 is.na(dat$both_anemia_overall_and_levels), ]
  }
  
  if ("included_in_model" %in% names(dat)) {
    dat <- dat[dat$included_in_model == 1 | is.na(dat$included_in_model), ]
  }
  
  # Sum number of "conf_" cols by group and subset to the most adjusted models
  if ("conf_altitude_adjusted" %in% names(dat)) {
    suppressWarnings(
      dat[is.na(conf_altitude_adjusted), conf_altitude_adjusted := 0]
    )
  } else {
    suppressWarnings(
      dat[, conf_altitude_adjusted := 0]
    )
  }
  
  if (outcome %in% c("lbw",
                     "lga",
                     "neo_mort",
                     "neo_sepsis",
                     "ptb",
                     "sga",
                     "stillbirth",
                     "mat_hem")) {
    conf_vars <- names(dat)[grep("^conf_", names(dat))]
    dat[, conf_count := rowSums(.SD == 1, na.rm = TRUE), .SDcols = conf_vars]
    
    dat[, retain := ifelse(conf_count == max(conf_count), 1, 0), 
        by = .(nid, age_start, age_end, rei, exp_timepoint, bio_type_1,
               alt_upper_1, alt_lower_1, ref_upper_1, ref_lower_1, out_type,
               out_lower, out_upper, effect_measure)]
    
    dat <- dat[retain == 1]  # only most adjusted models from each analysis
  }
  
  # Drop inapplicable effect_measures
  dat <- dat[!effect_measure == 'Chi Square Test']
  dat <- dat[!effect_measure == 'Mean difference']
  
  # Rid of inapplicable bio_units
  dat <- dat[!(bio_unit_1 %in% c('Z-score', "z-score"))]
  
  # Rid of inapplicable bio_forms
  if ("bio_form_1" %in% colnames(dat)) {
    dat <- dat[!(bio_form_1 %in% c("Continuous", "Continous"))]
  }
  
  # Subset by given value(s) of select covariate column(s), if specified
  dat <- subset_by_given_cov_value(dat, cov_con_socio, cov_rev)
  
  # Additional subsetting
  if (outcome %in% c("lbw",
                     "lga",
                     "neo_mort",
                     "neo_sepsis",
                     "ptb",
                     "sga",
                     "stillbirth",
                     "mat_hem")) {
    dat <- dat[exp_neonatal==0 | is.na(exp_neonatal)]    
  }
  
  # Perform additional cleaning -------------------------------------------
  
  # Prepare data for hb imputation
  if (outcome %in% c("mat_mort", "mat_sepsis", "preecl", "ghtn", "hdop")) {
    dat$rei <- ifelse(
      dat$risk == "Low hemoglobin (anemia or iron deficiency anemia)",
      376,
      95)
    dat$bio_type_1 <- ifelse(
      dat$bio_type_1 == "Packed cell volume",
      "packed_cell_volume",
      dat$bio_type_1) 
  }
  
  # Rename and dynamically create columns as needed to run BoP
  data.table::setnames(dat,
                       old = c("rei", "study_design", "out_measure"),
                       new = c("rei_id", "design", "measure"),
                       skip_absent = TRUE)
  
  dat[, rei := data.table::fifelse(
    rei_id == params$rei_id,
    "low_hgb",
    data.table::fifelse(rei_id == 95, "nutrition_iron", NA_character_)
  )]
  
  dat <- add_outcome_specific_cols(dat, outcome)
  
  # Subset data for the risk-outcome pair being modeled
  dat <- dat[cause_id == unique(dat$cause_id) &
               rei_id == params$rei_id &
               tolower(effect_measure) %in% tolower(params$effect_measure)]
  
  # Create and reformat additional columns as needed to run BoP
  dat$is_outlier <- 0
  dat$cv_pregnant <- 1  # add column to indicate pregnancies
  dat$measure <- tolower(dat$measure)
  dat$bio_unit_1 <- tolower(dat$bio_unit_1)
  dat$seq <- 1:nrow(dat)
  
  # Impute missing age values ---------------------------------------------
  
  cli::cli_progress_step("Imputing missing age values...",
                         msg_done = "Risk-outcome specific cleaning complete.")
  
  dat <- impute_missing_ages(dat)
  
  checkmate::assert_names(
    names(dat),
    must.include = req_cols
  )
  
  return(dat)
}
