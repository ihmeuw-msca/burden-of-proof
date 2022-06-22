#
# 01_create_template.R
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

# xiaochen's example
# model <- py_load_object(
# 	filename="/ihme/homes/xdai88/gbd_tobacco/gbd2019_alcohol/evidence_score/testing/test_run1_2020_09_05/04_monospline_pkl_files/lung_cancer_0.9_ensemble.pkl",
# 	pickle = "dill")
# data = model$data
# df <- data$to_df()

# diet example
# data <- readRDS(paste0(out_dir, "00_prepped_data/", ro_pair, ".RDS"))
# df <- data$df

library(readxl)
df_meta <- read_excel("/ihme/homes/jiaweihe/msca/mrbrt/evidence_score_pipeline/evidence_score.xlsx")
infile <- df_meta[df_meta$ro_pair==ro_pair, 'data']
df <- read.csv(infile[[1]])

# Deal with inconsistency of naming
OBS_VAR <- "ln_effect"
OBS_SE_VAR <- "ln_se"
STUDY_ID_VAR <- "nid"

if (!(STUDY_ID_VAR %in% names(df))){
	STUDY_ID_VAR <- "study_id"
}

if (!(OBS_VAR %in% names(df))){
	OBS_VAR <- "obs"
	if (!(OBS_VAR %in% names(df))){
		OBS_VAR <- "log_rr"
	}
}

if (!(OBS_SE_VAR %in% names(df))){
	OBS_SE_VAR <- "obs_se"
	if (!(OBS_SE_VAR %in% names(df))){
		OBS_SE_VAR <- "log_se"
	}
}

# Specify all the columns you need for your application
cov_names <- c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)

mrdata <- MRData()

mrdata$load_df(
  data = df, 
  col_obs = OBS_VAR,
  col_obs_se = OBS_SE_VAR, 
  col_study_id = STUDY_ID_VAR, 
  col_covs = as.list(cov_names)
)

monotonicity <- DIRECTION[ro_pair][[1]]
# if (is.na(monotonicity)){
# 	monotonicity <- NULL
# }

N_I_KNOTS <- 3
PRIOR_VAR_RSLOPE = 1e-6
PRIOR_VAR_MAXDER <- 1e-4

ensemble_cov_model <- LogCovModel(
	alt_cov = ALT_EXPOSURE_COLS,
  	ref_cov = REF_EXPOSURE_COLS,
	use_spline = TRUE,
	use_re = FALSE,
	spline_degree = 3L,
    spline_knots_type = 'domain',
    spline_r_linear = TRUE,
    prior_spline_funval_uniform = array(c(-1 + 1e-6, 19)),
    prior_spline_num_constraint_points = 150L,
    spline_knots = array(seq(0, 1, length.out = N_I_KNOTS + 2)),
    prior_spline_maxder_gaussian = cbind(rbind(rep(0, N_I_KNOTS), 
        rep(sqrt(PRIOR_VAR_MAXDER), N_I_KNOTS)), c(0, sqrt(PRIOR_VAR_RSLOPE))),
    prior_spline_der2val_gaussian = NULL,
    prior_spline_der2val_gaussian_domain = array(c(0.0, 1.0)),
    prior_spline_monotonicity = monotonicity
)

# Create knot samples
knots <- import("mrtool.core.model")
knots_samples <- knots$create_knots_samples(
  data = mrdata, l_zero = TRUE, num_splines = 50L, 
  num_knots = 5L, width_pct = 0.2,
  alt_cov_names = ALT_EXPOSURE_COLS,
  ref_cov_names = REF_EXPOSURE_COLS
)

# Ensemble model with exposure only 
signal_model <- MRBeRT(mrdata,
                      ensemble_cov_model=ensemble_cov_model,
                      ensemble_knots=knots_samples,
                      inlier_pct=0.9)

signal_model$fit_model(inner_print_level=5L, inner_max_iter=200L, 
	outer_step_size=200L, outer_max_iter=100L)

# create "new covariates" for later use
signal <- signal_model$predict(mrdata, predict_for_study=FALSE)

# Extract weights of data point
w <- t(do.call(rbind, 
          lapply(1:length(signal_model$sub_models), 
                 function(i){signal_model$sub_models[[i]]$w_soln}))
       ) %*% signal_model$weights

df_data <-  mrdata$to_df()
# Assign signal to data for use in later stage
df_data$signal <-  signal
# Drop data trimmed
df_data <-  df_data[w >= 0.1,]

# Save data and model
py_save_object(object = signal_model, 
	filename = paste0(out_dir, "01_template_pkl_files/", ro_pair, ".pkl"), 
	pickle = "dill")

out <- append(data, list(df_data=df_data))
saveRDS(out, paste0(out_dir, "01_template_models/", ro_pair, ".RDS"))
