#' Function to upload bundle and save crosswalk versions for BoP-prod

save_bundles_and_xwalks <- TRUE
path <- "PATH TO PARAMETER MAP"
paths_map <- data.table::fread(path)

# Upload crosswalk versions ----

if (save_bundles_and_xwalks) {
  box::use(
    risk_hemoglobin / lib / run_bop_prod_prep_hb,
    risk_hemoglobin / lib / rr_bundle_upload
  )
  
  cli::cli_progress_step(
    msg = "Preparing files for bop-prod",
    msg_done = "Completed file preparation",
    msg_failed = "Failed file preparation"
  )
  
  ## Copy files from bop-dev runs and add indicator columns to input data ----
  for (i in seq_len(nrow(paths_map))) {
    map_row <- paths_map[i, ]
    dat_path <- run_bop_prod_prep_hb$run_bop_prod_prep_hb(
      category = map_row$category,
      outcome = map_row$outcome,
      analysis = map_row$analysis,
      bun_id = map_row$bundle_id,
      trial_dir = map_row$trial_dir
    )
    
    paths_map[i, prod_dat_path := dat_path]
  }
  
  xwalk_version_ids <- c()
  
  ## Upload bundles and save crosswalk versions ----
  for (i in seq_len(nrow(paths_map))) {
    paths <- paths_map[i]
    
    cli::cli_progress_step(
      msg = paste("Saving xwalk version for bundle ", unique(paths$bundle_id)),
      msg_done = "Successfully saved crosswalk version",
      msg_failed = "Failed to save crosswalk version"
    )
    
    xwalk_version_id <- rr_bundle_upload$rr_bundle_upload(
      bun_id = paths$bundle_id,
      data_paths = paths$prod_dat_path,
      analysis_vec = paths$analysis,
      xwalk_path = "PATH TO PARENT FOLDER"
    )
    
    xwalk_version_ids <- c(xwalk_version_ids, xwalk_version_id)
  }
  paths_map$xwalk_version_id <- xwalk_version_ids
  
  data.table::fwrite(paths_map, path)
}