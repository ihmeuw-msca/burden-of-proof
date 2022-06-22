#
# 06_publication_bias.R
#
#
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")
args <- commandArgs(trailingOnly = TRUE)


### Running settings
# ro_pair <- args[1]
# out_dir <- args[2]
# WORK_DIR <- args[3]
ro_pair <- c("lpa_ihd")
out_dir <- ""
work_dir <- "/ihme/homes/zhengp/Repositories/evidence_score_pipeline"

setwd(work_dir)
source("./config.R")
source("./src/utils/loglinear_functions.R")

model_path <- paste0("/home/j/temp/zhengp/escore/sim_loglinear.pkl")
ref_covs <- NULL
alt_covs <- c("exp")


### Load model objects
model <- py_load_object(filename = model_path, pickle = "dill")

data_info <- extract_data_info(model,
                               ref_covs = ref_covs,
                               alt_covs = alt_covs)
data_info$ro_pair <- ro_pair
df <- data_info$df

### Detect publication bias
df_no_outlier <- df[!df$outlier,]
egger_model_all <- egger_regression(df$residual, df$residual_se)
egger_model <- egger_regression(df_no_outlier$residual, df_no_outlier$residual_se)
has_pub_bias <- egger_model$pval < 0.05

### Adjust for publication bias
if (has_pub_bias) {
  df_fill <- get_df_fill(df[!df$outlier,])
  num_fill <- nrow(df_fill)
} else {
  num_fill <- 0
}

# fill the data if needed and refit the model
if (num_fill > 0) {
  df <- rbind(df, df_fill)
  data_info$df <- df
  
  # refit the model
  data = MRData()
  data$load_df(
    data=df[!df$outlier,],
    col_obs='obs',
    col_obs_se='obs_se',
    col_covs=as.list(model$cov_names),
    col_study_id='study_id'
  )
  model_fill <- MRBRT(data, cov_models=model$cov_models)
  model_fill$fit_model()
} else {
  model_fill <- NULL
}

### Extract scores
uncertainty_info <- get_uncertainty_info(data_info, model)
if (is.null(model_fill)) {
  uncertainty_info_fill <- NULL
} else {
  uncertainty_info_fill <- get_uncertainty_info(data_info, model_fill)
}


### Output diagnostics
# figures
title <- paste0(ro_pair, ": egger_mean=", round(egger_model$mean, 3),
                ", egger_sd=", round(egger_model$sd,3), ", egger_pval=", 
                round(egger_model$pval, 3))
plot_residual(df, title)

plot_model(data_info,
           uncertainty_info,
           model,
           uncertainty_info_fill,
           model_fill)

# summary
summary <- summarize_model(data_info,
                           uncertainty_info,
                           model,
                           egger_model,
                           egger_model_all,
                           uncertainty_info_fill,
                           model_fill)
summary

draws <- get_draws(data_info, model)
