#' Launch script for parallel runs of BoP-prod 

path <- "PATH TO PARAMETER MAP"
paths_map <- data.table::fread(path)

nch::submit_job(
  script = "risk_hemoglobin/02_run_bop_prod_and_generate_outputs.R",
  array = seq_len(nrow(paths_map)),
  # memory = 10,
  # ncpus = 1,
  time = 60*3,
  partition = "long.q"
)
