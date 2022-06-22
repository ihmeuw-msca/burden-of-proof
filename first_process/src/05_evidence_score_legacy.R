#
# 05_evidence_score.R
#
#
library(dplyr)
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")

args <- commandArgs(trailingOnly = TRUE)

ro_pair <- args[1]
out_dir <- args[2]
WORK_DIR <- args[3]
setwd(WORK_DIR)
source("./config.R")


# Load signal_model and final_model
signal_model <- py_load_object(filename=paste0(out_dir, "01_template_pkl_files/", ro_pair, ".pkl"), 
  pickle = "dill")

final_model <- py_load_object(filename=paste0(out_dir, "04_mixed_effects_pkl_files/", ro_pair, ".pkl"), 
  pickle = "dill")

# using the scorelator

# need to run 'repl_python()' to open an interactive Python interpreter,
# then immediately type 'exit' to get back to the R interpreter
# -- this helps to load a required Python package
repl_python()
# -- type 'exit' or hit escape

evidence_score <- import("mrtool.evidence_score.scorelator")
scorelator <- evidence_score$ContinuousScorelator(signal_model = signal_model, final_model = final_model, 
                                                  alt_cov_names= as.list(ALT_EXPOSURE_COLS), 
                                                  ref_cov_names = as.list(REF_EXPOSURE_COLS),
                                                  name=ro_pair)
scorelator$plot_model(folder = paste0(out_dir, "05_evidence_score/"))
score <- scorelator$get_score()
low_score <- scorelator$get_score(use_gamma_ub=TRUE)
