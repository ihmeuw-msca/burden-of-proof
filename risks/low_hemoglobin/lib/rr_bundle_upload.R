#' Save RR bundles, bundle versions, and crosswalk versions
#'
#' @param bun_id A bundle ID to process.
#' @param data_paths A vector of data file paths to be processed and uploaded.
#' @param xwalk_path A path to which output files are saved.
#'
#' @note The code below can be used for troubleshooting when current xwalk data
#'   needs to be inspected:
#'   crosswalk_version_id <- 44564
#'   crosswalk_version_df <- ihme::get_crosswalk_version(crosswalk_version_id)
#'   crosswalk_version_df_tsh <- crosswalk_version_df |>
#'     as.data.frame() |>
#'     dplyr::filter(
#'     alt_risk_upper <= alt_risk_lower | ref_risk_upper <= ref_risk_lower
#'     )
#'
#' @example
#'   source("PATH TO SCRIPT")
#'   rr_bundle_upload(
#'     bun_id = 10507,
#'     data_paths = c(),
#'     xwalk_path = "PATH TO PARENT FOLDER"
#'   )

rr_bundle_upload <- function(bun_id, analysis_vec, data_paths, xwalk_path) {
  cli::cli_progress_step(
    msg = paste("Starting RR bundle upload step for bundle ID", bun_id),
    msg_done = paste("Completed RR bundle upload step for bundle ID", bun_id),
    msg_failed = paste("RR bundle upload step failed for bundle ID", bun_id),
  )

  bundle_version_id <- upload_bundle_and_save_bv(
    bun_id = bun_id,
    data_paths = data_paths,
    xwalk_path = xwalk_path
  )$bundle_version_id

    prep_and_save_xwalk_versions(
    bv_id = bundle_version_id,
    analysis_vec = analysis_vec,
    xwalk_path = xwalk_path
  )
}


#' Helper function to prep hemoglobin bundle data
prep_bundle_data <- function(data_paths, bun_id) {
  # Initialize an empty data frame
  bundle_data <- data.frame()

  # Read and combine data from multiple paths
  bundle_data <- lapply(data_paths, data.table::fread) |>
    data.table::rbindlist(use.names = TRUE, fill = TRUE)

  # Print unique values for manual review
  print(unique(bundle_data$acause))
  print(nrow(bundle_data))

  # Filter and print unique values for manual review
  bundle_data <- bundle_data |>
    dplyr::filter(cause_id != 1)
  print(unique(bundle_data$acause))
  print(unique(bundle_data$out_sub))
  print(unique(bundle_data$out_type))
  print(unique(bundle_data$cause_id))

  # Rename and transform columns
  bundle_data <- bundle_data |>
    dplyr::mutate(
      bundle_id = bun_id,
      seq = NA
    )

  # Add missing columns and set them to NA
  add_cols <- c(
    "underlying_nid",
    "upper",
    "lower",
    "input_type",
    "source_type_id"
  )
  bundle_data[, add_cols] <- NA

  # Rename columns and perform other prep
  if (!("cov_representativeness" %in% names(bundle_data))) {
    if ("cov_representative" %in% names(bundle_data)) {
      bundle_data <- bundle_data |>
        dplyr::rename(cov_representativeness = cov_representative)
    } else {
      bundle_data$cov_representativeness <- 1
    }
  }

  bundle_data <- bundle_data |>
    dplyr::mutate(
      design = tolower(design),
      is_outlier = 0,
      measure = "relrisk",
      standard_error = ln_rr_se,
      mean = ln_rr,
      effect_size_measure = "relative risk",
      effect_size_unit = "log",
      source_type = "Case notifications - infectious",
      sex = "Female"
    ) |>
    # Handle missing values in columns containing 'cov_'
    dplyr::mutate(dplyr::across(
      dplyr::contains('cov_'),
      tidyr::replace_na,
      0
    )) |>
    # Drop unnecessary columns
    dplyr::select(-dplyr::any_of(c("cov_representative", "bundle_version_id")))

  names(bundle_data)[names(bundle_data) == "study_id"] <- "nid"
  
  cli::cli_progress_done()

  return(bundle_data)
}


#' Helper function to upload data to bundle and save new bundle version
upload_bundle_and_save_bv <- function(bun_id, data_paths, xwalk_path) {
  # Save a backup of current bundle data
  current_bundle <- ihme::get_bundle_data(bundle_id = bun_id)
  path <- glue::glue(
    "{xwalk_path}/backup_current_bundle_{bun_id}_{Sys.Date()}.csv"
  )
  data.table::fwrite(current_bundle, path)

  # Empty the bundle if necessary, skip if there is no need to empty the bundle
  # TODO: recreate error and revise clear_bundle() to handle empty bundles
  if (nrow(current_bundle) > 0) {
    bundleprep::clear_bundle(bundle_id = bun_id)
  }

  ## Call prep_bundle_data() to combine R-O data ----
  bundle_data <- prep_bundle_data(data_paths = data_paths, bun_id = bun_id)

  # Write to XLSX file with sheet named "extraction"
  temp_directory <- withr::local_tempdir()
  temp_path <- file.path(temp_directory, "bundle.xlsx")
  writexl::write_xlsx(list(extraction = bundle_data), path = temp_path)

  # Upload data to bundle
  ihme::upload_bundle_data(
    bundle_id = bun_id,
    filepath = temp_path
  )
  current_bundle <- ihme::get_bundle_data(bundle_id = bun_id)
  checkmate::assertTRUE(nrow(current_bundle) == nrow(bundle_data))

  # Save bundle version ----
  bundle_version_id <- ihme::save_bundle_version(bundle_id = bun_id)

  return(bundle_version_id)
}


#' Helper function to prep xwalk versions and write XLSX files given a named list
prep_and_save_xwalk_versions <- function(bv_id, analysis_vec, xwalk_path) {
  # Prepare data from new bundle version for crosswalk upload
  bundle_version <- ihme::get_bundle_version(bundle_version_id = bv_id)
  bundle_version$crosswalk_parent_seq <- bundle_version$seq
  bundle_version$cov_reverse_causation <- NULL
  bundle_version$cov_outcome_quality <- NULL
  bundle_version$cov_exposure_quality <- NULL
  bundle_version$cov_confounder_quality <- NULL
  bundle_version$cov_selection_bias <- NULL

  xwalk_version_ids <- c() # store all version IDs
  
  # Save crosswalk versions by outcome and analysis as XLSX files for upload
  for (i in analysis_vec) {
    path <- file.path(xwalk_path, glue::glue("xwalk_{i}.xlsx"))
    bundle_version_filtered <- bundle_version |>
      dplyr::filter(.data[[i]] == 1) # filtered version only for writing
    
    if (nrow(bundle_version_filtered) == 0) {
      warning(glue::glue("No data found for analysis '{i}' in bundle version {bv_id}. Skipping."))
      next
    }
    
    openxlsx::write.xlsx(
      bundle_version_filtered,
      file = path,
      sheetName = "extraction"
    )
    
    description <- basename(path)
    xwalk_version <- ihme::save_crosswalk_version(
      bundle_version_id = bv_id,
      data_filepath = path,
      description = description
    )
    xwalk_version_id <- xwalk_version$crosswalk_version_id
    message(glue::glue(
      "Saved crosswalk version {xwalk_version_id} for bundle version ID {bv_id} analysis '{i}'."
    ))
    xwalk_version_ids <- c(xwalk_version_ids, xwalk_version_id)
  }
  return(xwalk_version_ids)
}
