#-----------------------------------------------------------------
# Purpose: cleaning the newly extracted tobacco data for mr-brt es
# Author: Xiaochen Dai
# Date: 07/25/2022
#-----------------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(data.table)
library(readxl)
library(stringr)

source('prep_data_function.R')
source("get_age_metadata.R")

# new data folder
new_rr_dir <- "[path to newly extracted data]"

# list files
new_rr_files <- list.files(new_rr_dir, full.names = T)

# remove temp files
new_rr_files <- new_rr_files[!grepl("~$",new_rr_files, fixed = T) & !grepl("template|map",new_rr_files)]

# read in the excel file
new_rr_data <- lapply(new_rr_files,read_excel, sheet = "extraction") %>% rbindlist(use.names = T, fill = T) %>% data.table

new_rr_data <- new_rr_data[!is.na(nid) & !grepl('uwnet',extractor, )]
new_rr_data[grepl('uwnet',extractor)]

# subset the risk-outcomes of interest
outcomes <- new_rr_data[, outcome] %>% unique %>% sort
write.csv(outcomes, paste0(new_rr_dir,'/outcome_map_v2.csv'), row.names = F)

# read in the outcome map
outcome_map <- fread(paste0(new_rr_dir,'/outcome_map.csv'))

# merge the ro_pair info
new_rr_data_ro <- merge(new_rr_data, outcome_map, by="outcome", all.x = T)

cause_map <- fread('map_causes.csv')
pairs <- cause_map[, cause_choice] %>% unique

pairs <- pairs[!grepl("fracture", pairs)]

# ad hoc fix, NID 343580 is pack-year not cig/day for prostate cancer
new_rr_data_ro <- new_rr_data_ro[!(ro_map=="prostate_cancer" & nid==343580),]

# fix a wrong nid for esophageal cancer
new_rr_data_ro <- new_rr_data_ro[nid==502488, nid := 502489]

for(ro_pair in pairs){
  print(ro_pair)
  # subset the data
  new_rr_data_sub <- new_rr_data_ro[ro_map==ro_pair,]
  if (nrow(new_rr_data_sub)==0) message(ro_pair, " has no new extraction")

  # ad hoc fix for bladder cancer
  if(ro_pair=="bladder_cancer"){
    # fix the lower effect size=0, the crude lower is 0.065. The extracted adjusted lower is 0.0. 
    # So, the adjusted lower should be around 0.03-0.05
    new_rr_data_sub[lower==0, lower := 0.04]
  }
  
  #-----------------
  ## clean the data
  #-----------------
  
  confounders <- names(new_rr_data_sub)[grepl('confounders_', names(new_rr_data_sub)) & !grepl('_other', names(new_rr_data_sub))]
  
  # change the data to numeric
  int_vars <- c("nid", "location_id", "year_start_study", "year_end_study", "age_start", "age_end", "rep_population", "sex_issue", "age_issue",
                confounders, "CI_uncertainty_type_value", "measure_issue", "effect_size_calculated", "subgroup_analysis", "uncertainty_issue",
                "effect_size_multi_location", "pooled_cohort", "confirm_unexp_reference", "cohort_number_events_exp", "cohort_number_events_unexp",
                "cohort_number_events_total", "cohort_sample_size_exp", "cohort_sample_size_unexp", "cohort_sample_size_total", "cc_community",
                "cc_cases_exp", "cc_control_exp", "cc_cases", "cc_sample_size")
  
  num_vars <- c("age_mean", "age_sd", "percent_male", "value_of_duration_fup", "effect_size", "lower", "upper", "exp_level_lower", "exp_level_upper",
                "exp_level_mid", "unexp_level_lower", "unexp_level_upper", "unexp_level_mid", "cohort_person_years_exp", "cohort_person_years_unexp", "cohort_person_years_total"
  )

  new_rr_data_sub[, (int_vars) := lapply(.SD, str_trim), .SDcols=int_vars]
  new_rr_data_sub[, (num_vars) := lapply(.SD, str_trim), .SDcols=num_vars]
  
  new_rr_data_sub[, (int_vars) := lapply(.SD, as.integer), .SDcols=int_vars]
  new_rr_data_sub[, (num_vars) := lapply(.SD, as.numeric), .SDcols=num_vars]
  
  # remove former smoking exposure
  new_rr_data_sub <- new_rr_data_sub[!grepl("former", exp_temporality, ignore.case = T)]
  
  #---------------------------
  # log-transform the effects
  #---------------------------
  
  new_rr_data_sub[, effect_size_log := log(effect_size)]
  new_rr_data_sub[, se_effect_log := abs((log(upper)-log(lower))/3.92)] 
  new_rr_data_sub[, se_effect_log] %>% summary
  
  # remove missing SE
  new_rr_data_sub <- new_rr_data_sub[!is.na(effect_size_log)]
  new_rr_data_sub <- new_rr_data_sub[!is.na(se_effect_log)]
  
  # creating bias covariates
  new_rr_data_sub[,cv_subpopulation := ifelse(rep_population==1, 0, 1)]
  new_rr_data_sub[,cv_exposure_selfreport := ifelse(exp_method_1=="Self-report (human/environment)",1,0)]
  new_rr_data_sub[,cv_outcome_selfreport := ifelse(outcome_assess_1=="Self-report",1,0)]
  new_rr_data_sub[,cv_exposure_study := ifelse(exp_assess_period=="only at baseline",1,0)]
  
  if(all(is.na(new_rr_data_sub[, age_start]))){
    new_rr_data_sub[!is.na(age_mean),cv_older := ifelse(age_mean>50,1,0)]
  } else {
    new_rr_data_sub[!is.na(age_start),cv_older := ifelse(age_start>50,1,0)]
  }
  
  new_rr_data_sub[,cv_non_smoker := ifelse(unexp_def=="never smokers",0,1)]
  new_rr_data_sub[,cv_risk_measure := ifelse(effect_size_measure=="Relative risk (RR)" | effect_size_measure=="Hazard ratio (HR)",0,1)]
  
  cvs <- names(new_rr_data_sub)[grepl('cv_', names(new_rr_data_sub))]
  
  # for confounders, replace NA with 0
  for (j in c(confounders)){
    set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,0)
  }
  
  # for cv covariates, replacing NA to 1 if necessary (very rare)
  for (j in c(cvs)){
    set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,1)
  }
  
  # replace missing age_start and age_end with the median value (this is fine for non-cvd outcomes since age is not relavent. May be problematic for CVDs though)
  set(new_rr_data_sub, which(is.na(new_rr_data_sub[["age_start"]])),"age_start", median(new_rr_data_sub[["age_start"]], na.rm = T))
  set(new_rr_data_sub, which(is.na(new_rr_data_sub[["age_end"]])),"age_end", median(new_rr_data_sub[["age_end"]], na.rm = T))
  new_rr_data_sub[is.na(percent_male), percent_male := 0.5]
  
  # check the exposure upper missing values
  new_rr_data_sub[, .(exp_level_lower, exp_level_upper, unexp_level_lower, unexp_level_upper, exp_level_mid, unexp_level_mid)]
  new_rr_data_sub[is.na(exp_level_lower),]
  
  new_rr_data_sub[is.na(unexp_level_lower) & is.na(unexp_level_upper) & !is.na(unexp_level_mid), c("unexp_level_lower", "unexp_level_upper") := unexp_level_mid]
  
  # change NA to 0 for unexposed group
  for (j in c("unexp_level_lower", "unexp_level_upper")){
    set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,0)
  }
  
  # change NA in exp_upper to 1.5*lower if missing
  new_rr_data_sub[is.na(exp_level_upper), exp_level_upper := 1.5*exp_level_lower]
  
  # count numbers of confounders
  new_rr_data_sub[, adjustment := apply(.SD, 1, sum), .SDcols=confounders]
  
  # count confounders in confounde_other
  new_rr_data_sub[, confounders_other := gsub(",", ";", confounders_other)]

  if(nrow(new_rr_data_sub[!is.na(confounders_other)])!=0){
    new_rr_data_sub[!is.na(confounders_other), num := unlist(lapply(str_split(confounders_other, ";"), length))]
  }
  
  new_rr_data_sub[is.na(confounders_other), num := 0]
  new_rr_data_sub[, adjustment:= adjustment + num]
  
  # creating cascading dummies
  new_rr_data_sub[ ,cv_adj := as.numeric()]
  
  new_rr_data_sub[,adj_age_sex := 0]
  new_rr_data_sub[confounders_age==1 & confounders_sex==1, adj_age_sex := 1]
  
  # cv_adj=3 if age or sex is not adjusted
  new_rr_data_sub[adj_age_sex==0, cv_adj := 3]
  
  # cv_adj=2 if only age and sex are adjusted
  new_rr_data_sub[adj_age_sex==1 & adjustment==2, cv_adj := 2]
  
  # cv_adj=1 if age+sex + 3 more covariates are adjusted
  new_rr_data_sub[adj_age_sex==1 & adjustment>2 & adjustment<=5, cv_adj := 1]
  
  # cv_adj=0 if age+sex + more than 3 covariates are adjusted
  new_rr_data_sub[adj_age_sex==1 & adjustment>5, cv_adj := 0]
  
  # check whether there is missing in cv_adj
  message("there is ", nrow(new_rr_data_sub[is.na(cv_adj)]), " missing values in cv_adj")
  
  # add cascading dummies
  new_rr_data_sub[, cv_adj_L0 := 1]
  new_rr_data_sub[, cv_adj_L1 := 1]
  new_rr_data_sub[, cv_adj_L2 := 1]
  
  # if cv_adj==0, change all dummies to be 0
  new_rr_data_sub[cv_adj==0, c("cv_adj_L0", "cv_adj_L1", "cv_adj_L2") := 0]
  # if cv_adj==1, change cv_adj_L1 and cv_adj_L2 to be 0
  new_rr_data_sub[cv_adj==1, c("cv_adj_L1", "cv_adj_L2") := 0]
  # if cv_adj==2, change cv_adj_L2 to be 0
  new_rr_data_sub[cv_adj==2, c("cv_adj_L2") := 0]
  
  # remove cv_adj
  new_rr_data_sub[, cv_adj := NULL]
  
  # change the exposed and unexposed name
  new_rr_data_sub <- setnames(new_rr_data_sub, old = c("exp_level_lower", "exp_level_upper", "unexp_level_lower", "unexp_level_upper"), 
                              new = c("b_0", "b_1", "a_0", "a_1"))
  
  new_rr_data_sub <- setnames(new_rr_data_sub, old = c("effect_size_log", "se_effect_log"), 
                              new = c("ln_effect", "ln_se"))
  
  # add mean exposure
  new_rr_data_sub[, mean_exp := (b_0+b_1)/2]
  
  
  # check data before saving
  bias_covs <- names(new_rr_data_sub)[grepl('cv_', names(new_rr_data_sub))]
  
  cvd_ro <- c("ihd", "stroke", "afib_and_flutter", "peripheral_artery_disease", "aortic_aneurism")
  
  if(ro_pair %in% cvd_ro){
    incl_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'age_start', 'age_end', "age_mean", "age_sd", "duration_fup_units",
                   "value_of_duration_fup", "duration_fup_measure", 'percent_male', bias_covs)
    nesc_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'percent_male', 'age_start', 'age_end', bias_covs)
  } else {
    incl_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'percent_male', bias_covs)
    nesc_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'percent_male', bias_covs)
  }
  
  
  if(nrow(new_rr_data_sub[ln_se<0,])>0) message(paste0("ln_se < 0 for ", nrow(new_rr_data_sub[ln_se<0,]), " row" ))
  if(nrow(new_rr_data_sub[lower>upper,])>0) message(paste0("lower and upper swapped for ", nrow(new_rr_data_sub[lower>upper,]), " rows" ))
  
  # ad hoc fixes
  if(ro_pair=="rheumatoid_arthritis"){
    new_rr_data_sub <- new_rr_data_sub[!is.na(b_1),]
  }
  
  # remove one study of cataracts since the exposure may be problematic
  if(ro_pair=="cataracts"){
    new_rr_data_sub <- new_rr_data_sub[nid!=501887,]
  }
  
  # remove one study of cancer since the exposure may be problematic
  if(ro_pair=="lung_cancer"){
    new_rr_data_sub <- new_rr_data_sub[nid!=502206,]
  }
  
  if(new_rr_data_sub[, nesc_vars, with=F] %>% is.na %>% any){
    message("there are missing values in the required variables... please check")
  } else {
    message("there are no missing values! Good to go!")
  }
  
  # keep usable variables
  new_data <- new_rr_data_sub[,incl_vars, with=F]
  
  # write the cleaned data
  write.csv(new_data, 
            paste0("[output path]"),
            row.names = F)
  
}