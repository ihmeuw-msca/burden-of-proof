#------------------------------------------------------------------
# Purpose: Estimate binary RR of smoking-fracture and create draws
# Author: Xiaochen Dai
# Date: 07/25/2022
#-----------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)
library(dplyr)
library(msm)
library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
source("/ihme/cc_resources/libraries/current/r/get_age_metadata.R")
args <- commandArgs(trailingOnly = TRUE)

## NEED TO CHANGE THE RO_PAIR HERE

ro_pair <- "fractures" # only works for fractures 
cov_setting <- "cov_finder_no_sex" # option: ['cov_finder', 'cov_finder_no_sex', 'no_cov','percent_male_only','self_selected'(no percent_male)]
trm <- 0.9
out_dir <- "[change to output directory]"

# load the config file
source("config.R")

cause <- ro_pair

# read in cleaned data
dt_sub <- fread(file.path(out_dir, paste0(cause, '_all_extraction.csv')))

# remove cv_exposure_selfreport
dt_sub[, cv_exposure_selfreport:= NULL]

# remove column with all 0s and all 1s
dt_sub <- dt_sub[, colSums(dt_sub != 0) > 0, with = FALSE]
dt_sub <- dt_sub[, colSums(dt_sub != 1) > 0, with = FALSE]

# 1. use cov finder to select significant covariates --------------------------------------------------------------------------------------------------------------------------------------------------------------

# covariate selection
if(cov_setting=="cov_finder_no_sex"){
  candidate_covs <-c(names(dt_sub)[grepl('cv_', names(dt_sub))])
} else {
  candidate_covs <-c(names(dt_sub)[grepl('cv_', names(dt_sub))],"percent_male")
}

candidate_covs <- candidate_covs[!grepl('confounding', candidate_covs)]

data <- MRData()
data$load_df(
  data=dt_sub,
  col_obs='ln_effect', #log space since bounded by 0
  col_obs_se='ln_se', 
  col_covs=as.list(candidate_covs),
  col_study_id='nid' 
)


if(cov_setting=="no_cov"){
  new_covs <- c()
  
  pred_data <- as.data.table(expand.grid('intercept'=c(1)))
  
} else if(cov_setting %in% c("cov_finder", 'cov_finder_no_sex')) {
  covfinder <- CovFinder(
    data = data,
    covs = as.list(candidate_covs),
    pre_selected_covs = list('intercept'),
    normalized_covs = FALSE,
    num_samples = 1000L,
    power_range = list(-4, 4),
    power_step_size = 0.05,
    inlier_pct = trm,
    laplace_threshold = 1e-5
  )
  
  covfinder$select_covs(verbose = FALSE)
  
  new_covs <- covfinder$selected_covs
  new_covs <- new_covs[!grepl('intercept', new_covs)]
  print(new_covs)
  
  if('percent_male' %in% new_covs){
    pred_data_m <- pred_data_f <- as.data.table(expand.grid("intercept"=c(1)))
    
    for(i in 1:length(new_covs)){
      pred_data_m <- cbind(pred_data_m,0)
      pred_data_f <- cbind(pred_data_f,0)
    }
    names(pred_data_m) <- c("intercept", new_covs)
    names(pred_data_f) <- c("intercept", new_covs)
    pred_data_m[, percent_male := 1]
    
  } else {
    pred_data <- as.data.table(expand.grid("intercept"=c(1)))
    for(i in 1:length(new_covs)){
      pred_data <- cbind(pred_data,0)
    }
    names(pred_data) <- c("intercept", new_covs)
  }
  
  
} else if(cov_setting=="percent_male_only"){
  new_covs <- c('percent_male')
  
  pred_data_m <- as.data.table(expand.grid("intercept"=c(1), "percent_male"=1))
  pred_data_f <- as.data.table(expand.grid("intercept"=c(1), "percent_male"=0))
  
} else {
  new_covs <- c('cv_adj_L1')
  
  pred_data <- as.data.table(expand.grid("intercept"=c(1), "cv_adj_L1"=0)) 
}

# 2. Prep data with covariates selected and run standard mixed effects model ------------------------------------------------------------------------------------------------------------------------------

if(cov_setting=="no_cov"){
  
  data1 <- MRData()
  data1$load_df(
    data=dt_sub,
    col_obs='ln_effect', #log space since bounded by 0
    col_obs_se='ln_se', 
    col_study_id='nid'
    )
  
} else {
  
  data1 <- MRData()
  data1$load_df(
    data=dt_sub,
    col_obs='ln_effect', #log space since bounded by 0
    col_obs_se='ln_se', 
    col_covs=as.list(new_covs),
    col_study_id='nid'
    )
  
} 

#Using use_re with covariates

cov_models <- list(
  LinearCovModel("intercept", use_re = T)
)

for (cov in new_covs) cov_models <- append(cov_models,
                                           list(do.call(
                                             LinearCovModel,
                                             c(list(alt_cov=cov, use_re = F)
                                             ))))

model <- MRBRT(
  data=data1, 
  cov_models=cov_models,
  inlier_pct=trm
)

model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

# 3. Make predictions ----------------------------------------------------------------------------------------------------------------------------------------------

if(cov_setting=="no_cov"){
  data_pred1 <- MRData()
  data_pred1$load_df(
    data = pred_data,
    col_covs = as.list("intercept")
    )
} else if(cov_setting %in% c("cov_finder", "percent_male_only")) {
  data_pred_m <- MRData()
  data_pred_f <- MRData()
  
  # predict male
  data_pred_m$load_df(
    data = pred_data_m,
    col_covs = as.list(new_covs)
  )
  
  # predict female
  data_pred_f$load_df(
    data = pred_data_f,
    col_covs = as.list(new_covs)
  )
} else {
  data_pred1 <- MRData()
  
  # predict 
  data_pred1$load_df(
    data = pred_data,
    col_covs = as.list(new_covs)
  )
}

# creating draws
n_samples <- 1000L
samples <- model$sample_soln(sample_size=n_samples)

if(cov_setting %in% c("no_cov", "self_selected", "cov_finder_no_sex")){
  draws <- model$create_draws(
    data = data_pred1, 
    beta_samples = samples[[1]],
    gamma_samples = samples[[2]],
    random_study=TRUE)
  
  pred_data$Y_mean <- model$predict(
    data = data_pred1, 
    predict_for_study = T, 
    sort_by_data_id = T) 
  
  pred_data$Y_mean_lo_re <- apply((draws), 1, function(x) quantile(x, 0.025)) 
  pred_data$Y_mean_hi_re <- apply((draws), 1, function(x) quantile(x, 0.975)) 
  
  #summary
  exp(pred_data$Y_mean) %>% print
  exp(pred_data$Y_mean_lo_re) %>% print
  exp(pred_data$Y_mean_hi_re) %>% print
  
} else if(cov_setting %in% c("cov_finder", "percent_male_only")){
  # create draws for male
  draws_m <- model$create_draws(
    data = data_pred_m, 
    beta_samples = samples[[1]],
    gamma_samples = samples[[2]],
    random_study=TRUE)
  
  pred_data_m$Y_mean <- model$predict(
    data = data_pred_m, 
    predict_for_study = T, 
    sort_by_data_id = T) 
  
  pred_data_m$Y_mean_lo_re <- apply((draws_m), 1, function(x) quantile(x, 0.025)) 
  pred_data_m$Y_mean_hi_re <- apply((draws_m), 1, function(x) quantile(x, 0.975))
  
  draws_fe_m <- model$create_draws(
    data = data_pred_m, 
    beta_samples = samples[[1]],
    gamma_samples = samples[[2]],
    random_study=FALSE)
  
  pred_data_m$Y_mean_lo_fe <- apply((draws_fe_m), 1, function(x) quantile(x, 0.025)) 
  pred_data_m$Y_mean_hi_fe <- apply((draws_fe_m), 1, function(x) quantile(x, 0.975)) 
  
  #summary
  print(exp(pred_data_m$Y_mean))
  print(exp(pred_data_m$Y_mean_lo_re))
  print(exp(pred_data_m$Y_mean_hi_re))
  
  # create draws for female
  draws_f <- model$create_draws(
    data = data_pred_f, 
    beta_samples = samples[[1]],
    gamma_samples = samples[[2]],
    random_study=TRUE)
  
  pred_data_f$Y_mean <- model$predict(
    data = data_pred_f, 
    predict_for_study = T, 
    sort_by_data_id = T) 
  
  pred_data_f$Y_mean_lo_re <- apply((draws_f), 1, function(x) quantile(x, 0.025)) 
  pred_data_f$Y_mean_hi_re <- apply((draws_f), 1, function(x) quantile(x, 0.975))
  
  draws_fe_f <- model$create_draws(
    data = data_pred_f, 
    beta_samples = samples[[1]],
    gamma_samples = samples[[2]],
    random_study=FALSE)
  
  pred_data_f$Y_mean_lo_fe <- apply((draws_fe_f), 1, function(x) quantile(x, 0.025)) 
  pred_data_f$Y_mean_hi_fe <- apply((draws_fe_f), 1, function(x) quantile(x, 0.975)) 

  
  #summary
  print(exp(pred_data_f$Y_mean))
  print(exp(pred_data_f$Y_mean_lo_re))
  print(exp(pred_data_f$Y_mean_hi_re))
  
}

# 4. Format and save  observation data for plotting ----------------------------------------------------------------------------------------------------------------------------------------
if (cov_setting %in% c('cov_finder', 'percent_male_only')){
  obs_data <- cbind("val" = data1$obs , "se"= data1$obs_se, "study" = data1$study_id, "included" = model$w_soln, 'sample_sex' = data1$covs$percent_male) %>% data.table
} else {
  obs_data <- cbind("val" = data1$obs , "se"= data1$obs_se, "study" = data1$study_id, "included" = model$w_soln) %>% data.table
}

obs_data[, lower:=val-1.96*se]
obs_data[, upper:=val+1.96*se]
obs_data <- obs_data[order(val)]
obs_data[, data:= 1]
obs_data[included > 0 & included < 1, included:=0.5]

if(cov_setting %in% c('cov_finder', 'percent_male_only')){
  # add results for male
  results_m <- data.table("val" = exp(pred_data_m$Y_mean), "study" = c("Result w/ gamma"), 
                          "lower" = exp(pred_data_m$Y_mean_lo_re), 
                          "upper" = exp(pred_data_m$Y_mean_hi_re))
  results_m[,data:= 2]
  results_m[, included := 1]
  results_m[, sample_sex := 1]
  
  # add results for female
  results_f <- data.table("val" = exp(pred_data_f$Y_mean), "study" = c("Result w/ gamma"), 
                          "lower" = exp(pred_data_f$Y_mean_lo_re), 
                          "upper" = exp(pred_data_f$Y_mean_hi_re))
  results_f[,data:= 2]
  results_f[, included := 1]
  results_f[, sample_sex := 0]
  
  obs_data <- rbindlist(list(results_m, results_f, obs_data), fill = T)
  obs_data[, row := 1:nrow(obs_data)]
  
} else {
  # add results for male
  results <- data.table("val" = exp(pred_data$Y_mean), "study" = c("Result w/ gamma"), 
                          "lower" = exp(pred_data$Y_mean_lo_re), 
                          "upper" = exp(pred_data$Y_mean_hi_re))
  results[,data:= 2]
  results[, included := 1]
  
  obs_data <- rbindlist(list(results, obs_data), fill = T)
  obs_data[, row := 1:nrow(obs_data)]
  
}

# save model objects and obs_data dt

save_dir <- "/mnt/team/team/pub/sub_risks/tobacco/code/xdai88/gbd2020_smoking/relative_risk_curves/binary_risk/fracture_binary/"
write.csv(obs_data, paste0(save_dir, cause,"_", cov_setting, '_', trm, '.csv'), row.names=F)
py_save_object(object = model, filename = paste0(save_dir, cause, '_', cov_setting, '_', trm, '.pkl'), pickle = "dill")

# 5. Save Draws --------------------------------------------------------------------------------------------------------------------------------------------------

# for GBD 2020, use draws that incorporate between-study heterogeneity (i.e. with gamma)

if(cov_setting %in% c("no_cov", "self_selected", "cov_finder_no_sex")){
  save_draws <- copy(draws) %>% data.table
  write.csv(save_draws, paste0(save_dir,cause, '_', cov_setting, '_', trm,'_','draws_raw.csv'), row.names = F)
} else {
  save_draws_m <- copy(draws_m) %>% data.table
  save_draws_f <- copy(draws_f) %>% data.table
  
  save_draws_m[, sex_id := 1]
  save_draws_f[, sex_id := 2]
  
  save_draws <- rbindlist(list(save_draws_m, save_draws_f), fill = T)
  
  write.csv(save_draws, paste0(save_dir,cause, '_', cov_setting, '_', trm,'_', 'draws_raw.csv'), row.names = F)
}
