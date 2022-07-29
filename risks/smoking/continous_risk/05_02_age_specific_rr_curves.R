#---------------------------------------------------
# Purpose: create final age_specific relative risk curves for CVD outcomes
# Author: Xiaochen Dai, adapted from Haley's codes
# Date: 07/25/2022
#---------------------------------------------------

rm(list = ls())

# System info
os <- Sys.info()[1]
user <- Sys.info()[7]

# Drives
j <- if (os == "Linux") "/home/j/" else if (os == "Windows") "J:/"
h <- if (os == "Linux") paste0("/homes/", user, "/") else if (os == "Windows") "H:/"

code_dir <- '[path to age_rr_utils.R]'
save_dir <- "[path to results]"
age_rr_dt_dir <- '[path to cleaned age-stratified data and AF draws]'
plot_dir <- "[path to plots]"

library(dplyr)
library(ggplot2)
library(data.table)
library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
source("get_ids.R")
source("get_age_metadata.R")
source(paste0(code_dir, "age_rr_utils.R"))
source("helper_functions.R")
np <- import("numpy")
np$random$seed(as.integer(123))

# Set up arguments
if(interactive()){
  ro_pair <- "[enter CVD risk-outcome pair of interest]"
  mean_factor <- F
  log_af <- F # using AF calculated based on log_rr
  level_100 <- F
} else {
  args <- commandArgs(trailingOnly = TRUE)
  ro_pair <- args[1]
  mean_factor <- as.logical(args[2])
  log_af <- as.logical(args[3])
  level_100 <- as.logical(args[4])
}

if(level_100){
  output_dir <- "[path to raw draws for 100 exposure levels]"
  new_dir <- "[path to final formated draws for 100 exposure levels]"
} else {
  output_dir <- "[path to raw draws for 1000 exposure levels]"
  new_dir <- "[path to final formated draws for 1000 exposure levels]"
}

ages <- get_age_metadata(19)
setnames(ages, c("age_group_years_start", "age_group_years_end"), c("age_start", "age_end"))
ages <- ages[,.(age_start, age_end, age_group_id)]

# get the reference age group
data <- readRDS(paste0(save_dir, "01_template_models/", ro_pair, ".RDS"))
df_data <- data$df_data
age_ref <- df_data$age_ref %>% mean
age_ref_group <- ages[age_start <= age_ref & age_end >= age_ref, age_group_id]

# load rr draws and age pattern draws
rr_draws <- fread(paste0(output_dir, "smoking_", ro_pair, ".csv"))
setnames(rr_draws, "risk", "exposure")
age_pattern_draws <- fread(file.path(age_rr_dt_dir, paste0("attenuation_pct_draws_", ro_pair, "_",age_ref_group,".csv")))

col_names <- c("exposure", paste0("draw_", 0:999))

# change variable name
names(rr_draws) <- col_names

# load mean of age pattern
age_pattern_mean <- fread(file.path(age_rr_dt_dir, paste0("attenuation_pct_summary_", ro_pair, "_",age_ref_group,".csv")))

# apply attenuation factors
if(mean_factor){
  # apply mean of attenuation factors only
  plot_path <- paste0(plot_dir, "age_spec_smoking_", ro_pair,"_af_mean.pdf")
  age_spec_risk_curve <- apply_age_pattern_mean_af(ro_pair =ro_pair, 
                                                   risk_curve_draws_df = rr_draws, 
                                                   age_pattern_mean_df = age_pattern_mean, 
                                                   age_pattern_mean_log_df=age_pattern_mean_log, 
                                                   log_af = log_af, 
                                                   plot = T, 
                                                   plot_path = plot_path)
  
} else {
  plot_path <- paste0(plot_dir, "age_spec_smoking_", ro_pair,"_af_draws_no_gamma.pdf")
  age_spec_risk_curve <- apply_age_pattern(ro_pair =ro_pair, 
                                           risk_curve_draws_df = rr_draws, 
                                           age_pattern_draws_df = age_pattern_draws, 
                                           #draws_in_log = T, 
                                           #return_draws_log = F, 
                                           plot = T, 
                                           plot_path = plot_path)
  
}

# re-shape the dataset
age_spec_rr <- melt(age_spec_risk_curve, id.vars = c("exposure", "age_group_id"), variable.name = "draw", value.name = "rr")
age_spec_rr[, draw:= as.numeric(draw)-1]

age_spec_rr_full <- copy(age_spec_rr)
for(age_id in c(6:8)){
  temp <- age_spec_rr[age_group_id==9]
  temp[, age_group_id := age_id]
  age_spec_rr_full <- rbindlist(list(temp, age_spec_rr_full), use.names = T)
}
age_spec_rr_full[, age_group_id] %>% unique

# add sex_id
age_spec_rr_full_m <- copy(age_spec_rr_full)
age_spec_rr_full_f <- copy(age_spec_rr_full)
age_spec_rr_full_m[, sex_id := 1]
age_spec_rr_full_f[, sex_id := 2]

age_spec_rr_full <- rbindlist(list(age_spec_rr_full_m, age_spec_rr_full_f), use.names = T)

setorder(age_spec_rr_full, "exposure","draw", "sex_id", "age_group_id")

# save the draws
message("saving draws...")
write.csv(age_spec_rr_full, paste0(new_dir, ro_pair, ".csv"), row.names = F)
