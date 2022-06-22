#
# run_pipeline_parallel.R
#

library(dplyr)
library(parallel)

#####
# user params
#
WORK_DIR <- "/ihme/homes/jiaweihe/msca/mrbrt/evidence_score_pipeline"
source(paste0(WORK_DIR, "/config.R"))
source(paste0(WORK_DIR, "/src/utils/prep_diet_data_function.R"))
source(paste0(WORK_DIR, "/src/utils/qsub_function.R"))


#####
# create directories
#

if (!dir.exists(OUT_DIR)) {
  if (!dir.exists(dirname(OUT_DIR))) {
    dir.create(dirname(OUT_DIR))
  }
  dir.create(OUT_DIR)
} else {
  warning("Directory '", OUT_DIR, "' already exists")
}

for (dir in SUB_DIRS) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  } else {
    warning("Directory '", dir, "' already exists")
  }
}


submit_jobs <- function(pair, WORK_DIR) {
  # stage 1, create signal
  submit_sub_job(pair, "01_create_template.R", "_01_template", WORK_DIR)
  qwait("01_create_template_models", pair)

  # # stage 2, loglinear model
  submit_sub_job(pair, "02_loglinear_models.R", "_02_loglinear", WORK_DIR)
  qwait("02_loglinear_models", pair)

  # stage 3, covariate selection
  submit_sub_job(pair, "03_covariate_selection.R", "_03_cov_selection", WORK_DIR)
  qwait("03_covariate_selection_models", pair)

  # stage 4, final model
  submit_sub_job(pair, "04_mixed_effects_models.R", "_04_mixed_effects", WORK_DIR)
}

#####
# data prep for diet risks
#

stage0_results <- lapply(RO_PAIRS, function(ro_pair) {
  x <- try({
    prep_diet_data(
      ro_pair = ro_pair,
      obs_var = OBS_VAR,
      obs_se_var = OBS_SE_VAR,
      ref_vars = REF_EXPOSURE_COLS,
      alt_vars = ALT_EXPOSURE_COLS,
      allow_ref_gt_alt = FALSE,
      diet_dir = INPUT_DATA_DIR,
      study_id_var = "nid",
      verbose = TRUE
    )
  })
  
  saveRDS(x, paste0(OUT_DIR, "00_prepped_data/", ro_pair, ".RDS"))
  return(x)
})

names(stage0_results) <- RO_PAIRS
saveRDS(stage0_results, paste0(OUT_DIR, "stage0_results.RDS"))

# Submit stage jobs for each pair
mclapply(RO_PAIRS, function(pair) {
  submit_jobs(pair, WORK_DIR)
}, mc.cores = length(RO_PAIRS))

