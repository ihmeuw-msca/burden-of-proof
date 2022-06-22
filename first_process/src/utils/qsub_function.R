
run_script <- function(script, img, args = c("")) {
  cmd <- paste(
    "/ihme/singularity-images/rstudio/shells/execRscript.sh",  
    "-i", img,  "-s", script, paste(args, collapse = " ")
  )
  system(cmd)
}

submit_qsub <- function(script, job_name, img, proj = PROJ,
                        queue = "long.q", hours = 6, threads = 1, 
                        error_logs = paste0("/share/temp/sgeoutput/", USER, "/errors"),
                        output_logs = paste0("/share/temp/sgeoutput/", USER, "/output"),
                        memory = "8G", args = "", verbose = TRUE) {
  cmd <- paste(
    "qsub -terse -N", job_name, 
    "-q", queue,
    paste0("-l fthread=", threads),
    paste0("-l m_mem_free=", memory),
    paste0("-l h_rt=", hours, ":00:00"),
    paste0("-l archive=TRUE"),
    "-P", proj,
    "-e", error_logs,
    "-o", output_logs,
    "/ihme/singularity-images/rstudio/shells/execRscript.sh ", 
    "-i", img,
    "-s", script, 
    paste(args, collapse = " ")
  )
  
  if (verbose) cat(cmd, "\n")
  system(cmd)
}

# Wait for upstream job to finish for each pair
qwait <- function(sub_dir, pair){
  outfile <- paste0(OUT_DIR, sub_dir, "/", pair, ".RDS")
  while (!file.exists(outfile)) {
    Sys.sleep(1)
  }
}

# Submit job for each stage for each pair
submit_sub_job <- function(pair, script, job_name_suffix, script_dir) {
  submit_qsub(
    script = paste0(CODE_PATH, script),
    job_name = paste0(pair, job_name_suffix), 
    img = SINGULARITY_IMG,
    args = c(pair, OUT_DIR, script_dir)
  )
}

# Submit job for plotting risk function and derivative fit.
submit_plot_job <- function(pair) {
  cmd <- paste0(
    paste0("sh ", CODE_PATH, "submit_qsub_python.sh "),
    paste0(CODE_PATH, " "),
    paste0(OUT_DIR, "04_monospline_pkl_files/ "),
    paste0(OUT_DIR, "05_monospline_pdf/ "),
    paste0(pair, " "),
    paste0(USER, " "),
    PROJ
  )
  system(cmd)
}

# Submit job for generating evidence score.
submit_score_job <- function(pair) {
  cmd <- paste0(
    paste0("sh ", CODE_PATH, "submit_score_qsub_python.sh "),
    paste0(CODE_PATH, " "),
    paste0(OUT_DIR, "05_monospline_pdf/ "),
    paste0(OUT_DIR, "06_evidence_score/ "),
    paste0(pair, " "),
    paste0(USER, " "),
    PROJ
  )
  system(cmd)
}