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
ro_pair <- c("sim")
out_dir <- ""
work_dir <- "/ihme/homes/zhengp/Repositories/evidence_score_pipeline"

setwd(work_dir)
source("./config.R")
source("./src/utils/dichotomous_functions.R")

model_path <- "/home/j/temp/zhengp/escore/sim_dicho.pkl"


### Load model objects
model <- py_load_object(filename = model_path, pickle = "dill")


### Extract data
df <- extract_data_info(model)


### Detect publication bias
egger_model_all <- egger_regression(df$residual, df$residual_se)
egger_model <- egger_regression(df[!df$outlier,]$residual, df[!df$outlier,]$residual_se)
has_pub_bias <- egger_model$pval < 0.05


### Adjust for publication bias
if (has_pub_bias) {
  df_fill <- get_df_fill(df)
  num_fill <- nrow(df_fill)
} else {
  num_fill <- 0
}

# fill the data if needed and refit the model
if (num_fill > 0) {
  df <- rbind(df, df_fill)
  
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
uncertainty_info <- get_uncertainty_info(model)
if (is.null(model_fill)) {
  uncertainty_info_fill <- NULL
} else {
  uncertainty_info_fill <- get_uncertainty_info(model_fill)
}


### Output diagnostics
plot_model(df, uncertainty_info, model, uncertainty_info_fill, model_fill, ro_pair)
summary <- summarize_model(ro_pair, model, model_fill, egger_model, egger_model_all, uncertainty_info)
draws <- get_draws(model)
