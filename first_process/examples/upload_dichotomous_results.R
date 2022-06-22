rm(list = ls())
source("src/upload_dichotomous.R")

results_folder <- "/mnt/team/msca/pub/archive/evidence-score-test/gbd2020-results"

pair_info <- list(
  opioid_suicide = list(
    rei_id = "unknown",
    cause_id = "unkown",
    model_path = file.path(results_folder, "opioid_suicide", "model.pkl")
  ),
  idu_hepB = list(
    rei_id = "unknown",
    cause_id = "unkown",
    model_path = file.path(results_folder, "idu_hepB", "model.pkl")
  ),
  idu_hepC = list(
    rei_id = "unknown",
    cause_id = "unkown",
    model_path = file.path(results_folder, "idu_hepC", "model.pkl")
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
