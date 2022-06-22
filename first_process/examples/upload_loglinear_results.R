rm(list = ls())
source("src/upload_loglinear.R")


pair_info <- list(
  air_no2_resp_asthma = list(
    rei_id = 404,
    cause_id = 515,
    risk_unit = "ppb",
    model_path = "/ihme/erf/GBD2020/air_no2/rr/models/20/model.pkl"
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
