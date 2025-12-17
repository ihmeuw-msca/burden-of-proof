#' Function to be parallel-run for BoP-prod and generating Omnibus outputs

library(data.table)

path <- "PATH TO PARAMETER MAP"
paths_map <- fread(path)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
map_row <- paths_map[task_id, ]

process_path_row <- function(map_row) {
  box::use(risk_hemoglobin/lib/run_bop_prod_post_hb)
  
  shell_path <- "PATH TO SHELL SCRIPT"
  
  bop_command <- glue::glue(
    "'continuous_pipeline -c {map_row$xwalk_version_id} -b {map_row$bundle_id} -i ./data -o ./results'"
  )
  
  bop_prod_dir <- file.path(
    "PATH TO PARENT FOLDER",
    basename(dirname(map_row$trial_dir))
  )
  
  system(paste(shell_path, bop_prod_dir, bop_command))
  
  ro_pair <- strsplit(map_row$trial_dir, "/")[[1]][7]
  ro_pair <- sub("hemoglobin", "hgb", ro_pair)
  
  results_files <- list.files(file.path(bop_prod_dir, "results", ro_pair), full.names = TRUE)
  results_updated_dir <- file.path(bop_prod_dir, "results_updated", ro_pair, "")
  
  unlink(results_updated_dir, recursive = TRUE)
  dir.create(results_updated_dir, recursive = TRUE)
  Sys.chmod(results_updated_dir, "775", use_umask = FALSE)
  file.copy(results_files, results_updated_dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
  
  run_bop_prod_post_hb$run_bop_prod_post_hb(
    dir = results_updated_dir,
    ro_pair = ro_pair,
    analysis = map_row$analysis,
    update_output = TRUE
  )
  
  outcome_severity <- ifelse(grepl("elbw|vlbw|mptb|eptb|vptb", map_row$outcome), map_row$outcome, 999)
  trimester <- ifelse(grepl("trim", map_row$analysis), stringr::str_extract(map_row$analysis, "\\d$"), 999)
  cov_con_socio <- ifelse(grepl("socio_conf_optimal", map_row$analysis), 0, 1)
  cov_rev <- ifelse(grepl("reverse_causal_agnostic", map_row$analysis), 999, 0)
  incl_undef <- ifelse(grepl("incl_undef", map_row$analysis), 1, 0)
  
  results_data_path <- list.files(
    path = results_updated_dir,
    pattern = "\\.csv$",
    full.names = TRUE
  ) |> 
    (
      \(x) x[!grepl("raw|draws|quantiles", x)]
    )()
  
  dat <- fread(results_data_path)
  df_list <- bophf:::flag_high_low_hb(dat)
  dat <- rbindlist(list(df_list$low_hb_dat, df_list$high_hb_dat), use.names = TRUE, fill = TRUE)
  plot_path <- paste0(results_updated_dir,"/dat_for_plotting.csv")
  fwrite(dat, plot_path)
  
  outcome <- map_row$outcome
  analysis <- map_row$analysis
  trial <- sub(".*trial_", "", map_row$trial_dir)
  
  nch::submit_job(
    script = "risk_hemoglobin/lib/plot_diagnostics_before_bop.R",
    script_args = c(
      plot_path, 
      outcome, 
      analysis,
      trial, 
      outcome_severity,
      trimester, 
      cov_con_socio, 
      cov_rev, 
      incl_undef,
      outputs_dir = results_updated_dir
    ),
    job_name = "plot_diagnostics_before_bop.R",
    memory = 5, 
    ncpus = 1, 
    time = 5, 
    partition = "long.q"
  )
}

process_path_row(map_row)
