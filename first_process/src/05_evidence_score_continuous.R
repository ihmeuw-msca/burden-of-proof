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
source("./src/utils/continuous_functions.R")

linear_model_path <- paste0("/home/j/temp/zhengp/escore/", ro_pair, "_linear.pkl")
signal_model_path <- paste0("/home/j/temp/zhengp/escore/", ro_pair, "_signal.pkl")
ref_covs <- c("a_0", "a_1")
alt_covs <- c("b_0", "b_1")


### Load model objects
linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")
signal_model <- py_load_object(filename = signal_model_path, pickle = "dill")

data_info <- extract_data_info(signal_model,
                               linear_model,
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
    col_covs=as.list(linear_model$cov_names),
    col_study_id='study_id'
  )
  linear_model_fill <- MRBRT(data, cov_models=linear_model$cov_models)
  linear_model_fill$fit_model()
} else {
  linear_model_fill <- NULL
}

### Extract scores
uncertainty_info <- get_uncertainty_info(data_info, linear_model)
if (is.null(linear_model_fill)) {
  uncertainty_info_fill <- NULL
} else {
  uncertainty_info_fill <- get_uncertainty_info(data_info, linear_model_fill)
}


### Output diagnostics
# figures
title <- paste0(ro_pair, ": egger_mean=", round(egger_model$mean, 3),
                ", egger_sd=", round(egger_model$sd,3), ", egger_pval=", 
                round(egger_model$pval, 3))
plot_residual(df, title)

plot_model(data_info,
           uncertainty_info,
           linear_model,
           signal_model,
           uncertainty_info_fill,
           linear_model_fill)

# summary
summary <- summarize_model(data_info,
                           uncertainty_info,
                           linear_model,
                           signal_model,
                           egger_model,
                           egger_model_all,
                           uncertainty_info_fill,
                           linear_model_fill)

draws <- get_draws(data_info, linear_model)
