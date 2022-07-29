rm(list = ls())
source("/ihme/homes/xdai88/gbd_tobacco/gbd2020_smoking/evidence_score_pipeline/src/upload_dichotomous.R")

ARCHIVE <- "[directory to the archive folder]"
out_dir <- "[directory to the outputs folder]"

pair_info <- list(
  smoking_hip_fracture = list(
    rei_id = "99",
    cause_id = "878",
    model_path = paste0(out_dir,"fracture_model.pkl")
  ),
  smoking_non_hip_fracture = list(
    rei_id = "99",
    cause_id = "923",
    model_path = paste0(out_dir,"fracture_model.pkl")
  )
)

for (pair in names(pair_info)) {
  print(paste0("upload pair=", pair))
  results_folder <- file.path(ARCHIVE, pair)
  if (!dir.exists(results_folder)) {
    dir.create(results_folder)
  }
  do.call(upload_results, c(pair_info[[pair]], list(results_folder = results_folder)))
}
