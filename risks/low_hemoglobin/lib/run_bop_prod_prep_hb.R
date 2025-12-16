#' Prepare crosswalk versions for hemoglobin BoP
#'
#' @example
#' path_map <- "PATH TO PARAMETER MAP"
#' map <- data.table::fread(path_map)[1, ]
#' run_bop_prod_prep_hb <- function(
#'     category = map$category,
#'     outcome = map$outcome,
#'     analysis = map$analysis,
#'     bun_id = map$bun_id,
#'     trial_dir = map$trial_dir
#' )

run_bop_prod_prep_hb <- function(
  category,
  outcome,
  analysis,
  bun_id,
  trial_dir
) {
  # Prepare the bop-prod dir and copy over data and settings files from bop-dev
  dat_path <- copy_bop_dev_files(
    category = category,
    outcome = outcome,
    analysis = analysis,
    bun_id = bun_id,
    trial_dir = trial_dir
  )

  return(dat_path)
}


#' Add columns needed for crosswalking based on output directory name
add_indicator_cols <- function(dat, outcome, analysis) {
  dat <- dat |>
    dplyr::mutate(
      # add col indicating outcome subtype
      out_sub = outcome,
      # add cols indicating any trimester sensitivity
      trim1 = as.integer(grepl("trim1", analysis)),
      trim2 = as.integer(grepl("trim2", analysis)),
      trim3 = as.integer(grepl("trim3", analysis)),
      # add cols indicating any bias cov sensitivity (reverse causality)
      reverse_causal_agnostic_incl_severe_preecl_ecl = as.integer(
        grepl("reverse_causal_agnostic_incl_severe_preecl_ecl", analysis)
      ),
      reverse_causal_agnostic_incl_severe_preecl = as.integer(
        grepl("reverse_causal_agnostic_incl_severe_preecl", analysis) &
          !grepl("reverse_causal_agnostic_incl_severe_preecl_ecl", analysis)
      ),
      reverse_causal_agnostic = as.integer(
        grepl("reverse_causal_agnostic", analysis) &
          !grepl("reverse_causal_agnostic_incl_severe_preecl", analysis) &
          !grepl("reverse_causal_agnostic_incl_severe_preecl_ecl", analysis)
      ),
      # add cols indicating any bias cov sensitivity (socioeconomic confounding)
      socio_conf_optimal = as.integer(grepl("socio_conf_optimal", analysis)),
      # add cols indicating any undefined outcome sensitivity
      incl_undef = as.integer(
        (grepl("incl_undef", analysis) | grepl("timing_20", analysis)) &
          (grepl("ptb", outcome) |
            grepl("lbw", outcome) |
            grepl("sga", outcome) |
            grepl("lga", outcome) |
            grepl("stillbirth", outcome))
      ),

      overall = as.integer(
        trim1 == 0 &
          trim2 == 0 &
          trim3 == 0 &
          reverse_causal_agnostic == 0 &
          reverse_causal_agnostic_incl_severe_preecl == 0 &
          reverse_causal_agnostic_incl_severe_preecl_ecl == 0 &
          socio_conf_optimal == 0 &
          incl_undef == 0
        # 1 if none of the other indicator cols have 1s
      )
    )

  return(dat)
}


#' Prepare the bop-prod dir and copy over data and settings files from bop-dev
copy_bop_dev_files <- function(category, outcome, analysis, bun_id, trial_dir) {
  # Check that required input files exist at given bop-dev path
  checkmate::assertDirectoryExists(trial_dir)
  files <- list.files(
    path = trial_dir,
    pattern = "\\.csv$",
    full.names = TRUE
  )
  dev_dat_path <- files[grepl("raw", files)]
  checkmate::assertFileExists(dev_dat_path, .var.name = "dev_dat_path")
  checkmate::assertFileExists(file.path(trial_dir, "settings.yaml"))

  # Create the bop-prod dir, if needed
  bop_prod_dir <- file.path(
    "PATH TO PARENT FOLDER",
    basename(dirname(trial_dir))
  )
  if (!dir.exists(bop_prod_dir)) {
    dir.create(bop_prod_dir, recursive = TRUE)
    Sys.chmod(bop_prod_dir, "775", use_umask = FALSE)
  }
  checkmate::assertDirectoryExists(bop_prod_dir)

  # Clear any existing files in the bop-prod dir
  bop_prod_data_dir <- file.path(bop_prod_dir, "data")
  unlink(bop_prod_data_dir, recursive = TRUE)
  bop_prod_results_dir <- file.path(bop_prod_dir, "results")
  unlink(file.path(bop_prod_dir, "results"), recursive = TRUE)

  # Copy over settings.yaml file
  dir.create(bop_prod_data_dir, recursive = TRUE)
  Sys.chmod(bop_prod_data_dir, "775", use_umask = FALSE)
  checkmate::assertDirectoryExists(bop_prod_data_dir)

  file.copy(
    file.path(trial_dir, "settings.yaml"),
    bop_prod_data_dir,
    recursive = TRUE,
    overwrite = TRUE
  )
  settings_path <- file.path(bop_prod_data_dir, "settings.yaml")
  update_settings_yaml(settings_path, trial_dir)
  checkmate::assertFileExists(settings_path)

  # Read input data from bop-dev location, revise if needed, then write to bop-prod location
  dat <- data.table::fread(dev_dat_path)
  dat[,
    alt_risk_lower := data.table::fifelse(
      alt_risk_lower == alt_risk_upper,
      alt_risk_lower - 0.0000001, # small change to pass validation
      alt_risk_lower
    )
  ]

  dat_path <- file.path(bop_prod_data_dir, basename(dev_dat_path))
  data.table::fwrite(dat, dat_path)
  checkmate::assertFileExists(dat_path)

  # Create a results dir
  dir.create(bop_prod_results_dir, recursive = TRUE)
  Sys.chmod(bop_prod_results_dir, "775", use_umask = FALSE)
  checkmate::assertDirectoryExists(bop_prod_results_dir)

  # If needed, add indicator columns and overwrite file at input data path
  dat <- data.table::fread(dat_path)
  if (!all(colnames(dat) %in% c("out_sub", "trim1", "trim2", "trim3"))) {
    dat <- add_indicator_cols(
      dat = dat,
      outcome = outcome,
      analysis = analysis
    )

    data.table::fwrite(dat, dat_path)
  }

  return(dat_path)
}

# Edit settings.yaml format for BoP-prod run

update_settings_yaml <- function(settings_path, trial_dir) {
  ## Create new sections
  settings_yaml <- readLines(settings_path)

  ro_pair <- strsplit(trial_dir, "/")[[1]][7]
  ro_pair <- sub("hemoglobin", "hgb", ro_pair)
  divide <- which(stringr::str_detect(settings_yaml, "^metadata:"))

  top_section <- settings_yaml[1:(divide - 1)]
  bottom_section <- settings_yaml[divide:length(settings_yaml)]
  
  top_section <- gsub("^([ ]*num_samples:)[ ]*[0-9]+", "\\1 50", top_section)
  
  top_section_indented <- paste0("  ", top_section)
  
  bottom_section_indented <- paste0("  ", bottom_section)

  ## Combine into final output & overwrite
  combined_file <- c(
    paste0(ro_pair, ":"),
    top_section_indented,
    "",
    "default:",
    bottom_section_indented
  )

  writeLines(combined_file, settings_path)
}
