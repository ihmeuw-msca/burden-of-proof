#
# 04_mixed_effects_models.R
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

# Extract selected covariates
data <- readRDS(paste0(out_dir, "03_covariate_selection_models/", ro_pair, ".RDS"))
df_data <- data$df_data
df_tmp <- data$df
# Only keep rows that are not trimmed
df_tmp <- df_tmp[as.numeric(rownames(df_data)),]

cov_names <- data$selected_covs
bias_covs <- cov_names[!cov_names == "exposure_linear"]

# Add interaction
for (cov in bias_covs) df_data[, cov] <- df_data$signal * df_tmp[, cov]

# Selected bias covariates plus signal
covs <- c("signal", bias_covs)

mrdata <- MRData()
mrdata$load_df(
  df_data,
  col_obs = c('obs'),
  col_obs_se = c('obs_se'), 
  col_study_id = c('study_id'),
  col_covs=as.list(covs)
)


loglinear_model <- readRDS(paste0(out_dir, "02_loglinear_models/", ro_pair, ".RDS"))

# Beta prior from first loglinear model results.
beta_gprior_std <- loglinear_model$beta_std

# Combine cov models
cov_models <- list()
for (cov in bias_covs) cov_models <- append(cov_models, 
  list(
    do.call(
      LinearCovModel, 
      list(
        alt_cov=cov,
        beta_gprior_std=BETA_PRIOR_MULTIPLIER * beta_gprior_std
      )
    )
  )
)

# Mixed effects model
cov_models <- append(cov_models, LinearCovModel('signal', use_re=TRUE, 
  prior_beta_uniform=array(c(1.0, 1.0))))

model <- MRBRT(
  data=mrdata,
  cov_models = cov_models,
  inlier_pct = 1.0
)

model$fit_model(inner_print_level=5L, inner_max_iter=200L, 
  outer_step_size=200L, outer_max_iter=100L)

# Load signal model and data in Stage 1
signal_model <- py_load_object(filename=paste0(out_dir, "01_template_pkl_files/", ro_pair, ".pkl"), 
  pickle = "dill")
orig_data <- readRDS(paste0(out_dir, "01_template_models/", ro_pair, ".RDS"))
df <- orig_data$df

# This should be provided by the user
NUM_POINTS <- 100L
exposure_lower <- min(df[,c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)])
exposure_upper <- max(df[,c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)])
exposure <- seq(exposure_lower, exposure_upper, length.out=NUM_POINTS)
min_cov <- rep(exposure_lower, NUM_POINTS)

# Deal with Sarah's data
if ('a_0' %in% REF_EXPOSURE_COLS){
  df_signal_pred <- data.frame(a_0=min_cov, a_1=min_cov, b_0=exposure, b_1=exposure)
} else {
  df_signal_pred <- data.frame(a_0=min_cov, b_0=exposure)
  names(df_signal_pred) <- c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)
}

# Predict using signal model and gridded exposure
data_signal_pred <- MRData()
data_signal_pred$load_df(
  df_signal_pred,
  col_covs = as.list(c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS))
  )
signal_pred <- signal_model$predict(data_signal_pred)

# TODO: data of selected covariates to be added
df_final_pred <- data.frame(signal=signal_pred)
data_final_pred <- MRData()
data_final_pred$load_df(
  df_final_pred,
  col_covs = as.list(c("signal"))
)

# create draws and prediction
sampling <- import("mrtool.core.other_sampling")
num_samples <- 1000L
beta_samples <- sampling$sample_simple_lme_beta(num_samples, model)
gamma_samples <- rep(model$gamma_soln, num_samples) * matrix(1, num_samples)

curve <- model$predict(data_final_pred)
draws <- model$create_draws(
  data_final_pred,
  beta_samples=beta_samples,
  gamma_samples=gamma_samples
)

# Save model
py_save_object(object = model, 
  filename = paste0(out_dir, "04_mixed_effects_pkl_files/", ro_pair, ".pkl"), 
  pickle = "dill")

# OBS_VAR <- "ln_effect"
# OBS_SE_VAR <- "ln_se"

# if (!(OBS_VAR %in% names(df))){
#   OBS_VAR <- "obs"
# }

# if (!(OBS_SE_VAR %in% names(df))){
#   OBS_SE_VAR <- "obs_se"
# }

# Sanity check
# pdf(paste0(out_dir, "04_mixed_effects_models/", ro_pair, ".pdf"))

# if (length(ALT_EXPOSURE_COLS) == 1){
#   plot(df[,ALT_EXPOSURE_COLS] - df[,REF_EXPOSURE_COLS], 
#      df[, OBS_VAR], cex=1/(7*df[, OBS_SE_VAR]), xlab="exposure", ylab="ln_effect",
#      main=ro_pair, col=c('blue'))
# } else {
#   plot(apply(df[,ALT_EXPOSURE_COLS], 1, mean) - apply(df[,REF_EXPOSURE_COLS], 1, mean), 
#      df[, OBS_VAR], cex=1/(7*df[, OBS_SE_VAR]), xlab="exposure", ylab="ln_effect",
#      main=ro_pair, col=c('blue'))
# }

# lines(exposure, curve)
# dev.off()