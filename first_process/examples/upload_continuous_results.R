rm(list = ls())
source("src/upload_continuous.R")

results_folder <- "/mnt/team/msca/pub/archive/evidence-score-test/gbd2020-results"

pair_info <- list(
  dairy_diabetes = list(
    rei_id = "unknown",
    cause_id = "unkown",
    risk_unit = "unknown",
    signal_model_path = file.path(results_folder, "dairy_diabetes", "signal_model.pkl"),
    linear_model_path = file.path(results_folder, "dairy_diabetes", "linear_model.pkl")
  ),
  air_pmhap_lri = list(
    rei_id = "unknown",
    cause_id = "unkown",
    risk_unit = "unknown",
    signal_model_path = file.path(results_folder, "air_pmhap_lri", "signal_model.pkl"),
    linear_model_path = file.path(results_folder, "air_pmhap_lri", "linear_model.pkl")
  ),
  fpg_neo_liver = list(
    rei_id = "unknown",
    cause_id = "unkown",
    risk_unit = "unknown",
    signal_model_path = file.path(results_folder, "fpg_neo_liver", "signal_model.pkl"),
    linear_model_path = file.path(results_folder, "fpg_neo_liver", "linear_model.pkl")
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
