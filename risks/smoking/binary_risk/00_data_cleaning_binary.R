#------------------------------------------------------------------
# Purpose: Cleaning binary RR data of smoking-fracture for MRBRT
# Author: Xiaochen Dai
# Date: 07/25/2022
#-----------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)
library(dplyr)
library(msm)
library(readxl)

library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
source("/ihme/cc_resources/libraries/current/r/get_age_metadata.R")
library(gargle, lib.loc = '/homes/xdai88/rlibs')
library("googledrive", lib.loc = '/homes/xdai88/rlibs')

cause <- "fractures"
out_dir <- "/mnt/team/team/pub/sub_risks/tobacco/code/xdai88/gbd2020_smoking/relative_risk_curves/cleaned_new_data/"

#--------------------
# load extracted data
#--------------------
jason_path <- "[directory to Jason's data]"
gabi_path <- "[directory to Gabi's data]"
chuks_path <- "[directory to Chukwuma's data]"
xdai_path <- "[directory to Xiaochen's data]"

dt_jason <- read_xlsx(jason_path, sheet = "extraction") %>% data.table
dt_gabi <- read_xlsx(gabi_path, sheet = "extraction") %>% data.table
dt_chuks <- read_xlsx(chuks_path, sheet = "extraction") %>% data.table
dt_xdai <- read_xlsx(xdai_path, sheet = "extraction") %>% data.table

# remove the first row and missing rows
dt_jason <- dt_jason[-1,]; dt_gabi <- dt_gabi[-1,]; dt_chuks <- dt_chuks[-1,]; dt_xdai <- dt_xdai[-1,]
dt_jason <- dt_jason[!is.na(nid),]; dt_gabi <- dt_gabi[!is.na(nid),]; dt_chuks <- dt_chuks[!is.na(nid),]; dt_xdai <- dt_xdai[!is.na(nid),]

#-------------------------------------
# load in new extraction for fracture
#-------------------------------------
new_rr_dir <- "[directory to newly extracted data]"

# list files
new_rr_files <- list.files(new_rr_dir, full.names = T)

# remove temp files
new_rr_files <- new_rr_files[!grepl("~$",new_rr_files, fixed = T) & !grepl("template|map",new_rr_files)]

# read in the excel file
new_rr_data <- lapply(new_rr_files,read_excel, sheet = "extraction") %>% rbindlist(use.names = T, fill = T) %>% data.table
new_rr_data <- new_rr_data[!is.na(nid) & !grepl('uwnet',extractor, )]

# read in the outcome map (please manually check this map)
outcome_map <- fread(paste0(new_rr_dir,'/outcome_map.csv'))
# merge the ro_pair info
new_rr_data_ro <- merge(new_rr_data, outcome_map, by="outcome", all.x = T)

# subset the data
ro_pair <- "fractures"
new_rr_fracture <- new_rr_data_ro[ro_map=="fractures",]

same_cols <- names(dt_xdai)[names(dt_xdai) %in% names(new_rr_fracture)]
diff_cols_old <- names(dt_xdai)[!names(dt_xdai) %in% names(new_rr_fracture)]
diff_cols_new <- names(new_rr_fracture)[!names(new_rr_fracture) %in% names(dt_xdai)]

# cols in new data need to change name
new_names <- c("rep_population", "year_start_study", "year_end_study", "sex_categorical", "effect_size", 
               "CI_uncertainty_type_value", "exposed_def","unexp_def")
change_to_old_names <- c("representativeness", "year_start", "year_end", "sex", "mean", 
                     "uncertainty_type_value","cohort_exposed_def", "cohort_unexp_def")

setnames(new_rr_fracture, new_names, change_to_old_names)

# add an indicator of new data
new_rr_fracture[, new_extraction := 1]

# combine data together
dt_old <- rbindlist(list(dt_jason, dt_gabi, dt_chuks, dt_xdai), use.names = T)
dt_old[, new_extraction := 0]
dt_orig <- rbindlist(list(dt_old, new_rr_fracture), use.names = T, fill = T)

confounders <- names(dt_orig)[grepl('confounders_', names(dt_orig)) & !grepl('_other', names(dt_orig))]

names(dt_orig)
str(dt_orig)

# transform str to numbers
var_int <- c("nid", "underlying_nid", "location_id", "smaller_site_unit", "representativeness", "proxy_responses", "year_start", "year_end", "year_issue", 
             "age_start", "age_end", "age_issue", confounders, "exp_assess_num", "subgroup_analysis", "effect_size_multi_location",
             "most_adj_model", "least_adj_model", "uncertainty_issue", "effect_size_calculated", "included_excluded", "reverse_causation",
             "cohort_person_years_exp",	"cohort_person_years_unexp",	"cohort_person_years_total",	"cohort_number_events_exp",	"cohort_number_events_unexp",	"cohort_number_events_total",	
             "cohort_sample_size_exp",	"cohort_sample_size_unexp",	"cohort_sample_size_total",
             "cc_community",	"cc_cases",	"cc_control",	"cc_response_rate",	"cc_response_assess")

var_num <- c("age_mean", "age_sd", "percent_male", "exp_recall_period_value", "value_of_duration_fup", "mean",	"lower",	"upper",	"uncertainty_type_value", "standard_error",
             "cc_response_rate", "cohort_retention_rate")

dt_orig[, (var_int) := lapply(.SD, as.integer), .SDcols=var_int]
dt_orig[, (var_num) := lapply(.SD, as.numeric), .SDcols=var_num]

dt_orig[,nid] %>% unique

# for 394006, remove the duplicated data
dt_394006_used <- dt_orig[nid==394006 & !is.na(included_excluded)]
dt_orig_excl_394006 <- dt_orig[nid!= 394006]
dt_orig <- rbindlist(list(dt_orig_excl_394006, dt_394006_used), use.names = T)

#---------------------------- 
#create risk of bias covs
#----------------------------
# representativeness
if(dt_orig[is.na(exp_assess_level)] %>% nrow){
  message("there is missing value in cv_population")
} else {
  message("there is no missing value in cv_population. Good to go!")
}

dt_orig[!is.na(representativeness), cv_subpopulation := ifelse(representativeness==0, 1, 0)]

# Exposure population vs individual
if(dt_orig[is.na(exp_assess_level)] %>% nrow){
  message("there is missing value in cv_exposure_population")
} else {
  message("there is no missing value in cv_exposure_population. Good to go!")
}

dt_orig[!is.na(exp_assess_level), cv_exposure_population := ifelse(exp_assess_level=="At the individual", 0, 1)]

# Exposure self-report vs. not
if(dt_orig[is.na(exp_method_1) & is.na(exp_method_2) & is.na(exp_method_3),] %>% nrow){
  message("there is missing value in cv_exposure_selfreport")
} else {
  message("there is no missing value in cv_exposure_selfreport Good to go!")
}

dt_orig[!is.na(exp_method_1) & is.na(exp_method_2) & is.na(exp_method_3), cv_exposure_selfreport := ifelse(exp_method_1=="Self-report (human/environment)" |
                                                                                                           exp_method_2=="Self-report (human/environment)" |
                                                                                                           exp_method_3=="Self-report (human/environment)", 1, 0)]

# exposure measured multiple times or not
if(dt_orig[is.na(exp_assess_period),] %>% nrow){
  message("there is missing value in cv_exposure_study")
} else {
  message("there is no missing value in cv_exposure_study. Good to go!")
}

dt_orig[!is.na(exp_assess_period), cv_exposure_study := ifelse(exp_assess_period=="only at baseline", 1, 0)]

#	Outcome self-report vs. not
if(dt_orig[is.na(outcome_assess_1) & is.na(outcome_assess_2) & is.na(outcome_assess_3),] %>% nrow){
  message("there is missing value in cv_outcome_selfreport")
} else {
  message("there is no missing value in cv_outcome_selfreport. Good to go!")
}

dt_orig[!is.na(exp_method_1) & is.na(exp_method_2) & is.na(exp_method_3), cv_outcome_selfreport := ifelse(outcome_assess_1=="Self-report" |
                                                                                                            outcome_assess_2=="Self-report" |
                                                                                                            outcome_assess_3=="Self-report", 1, 0)]
# Reverse causation or not
if(dt_orig[is.na(reverse_causation),] %>% nrow){
  message("there is missing value in cv_reverse_causation")
} else {
  message("there is no missing value in cv_reverse_causation Good to go!")
}

dt_orig[!is.na(reverse_causation), cv_reverse_causation := ifelse(reverse_causation==1, 1, 0)]

# selection_bias or not
if(dt_orig[is.na(cc_response_rate) & is.na(cohort_retention_rate),] %>% nrow){
  message(paste0("there is missing value in cv_selection_bias. ",nrow(dt_orig[is.na(cc_response_rate) & is.na(cohort_retention_rate),]), " rows missing"))
} else {
  message("there is no missing value in cv_selection_bias Good to go!")
}
# missing values are considered to have cv_selection_bias=1
# using cascading dummies. 
dt_orig[, cv_selection_bias_cat1 := ifelse((cc_response_rate>=0.95 | cohort_retention_rate>=0.95), 0, 1)]
dt_orig[, cv_selection_bias_cat2 := ifelse((cc_response_rate>=0.85 | cohort_retention_rate>=0.85), 0, 1)]
dt_orig[is.na(cv_selection_bias_cat1), cv_selection_bias_cat1:=1]
dt_orig[is.na(cv_selection_bias_cat2), cv_selection_bias_cat2:=1]

# adding cv_older to indicate estiamtes only among older people. 
dt_orig[ , cv_older := ifelse(age_start>50, 1, 0)]
dt_orig[age_start>50, cv_older := 1]

# add cv_risk_measure to indicate RR or not, RR is 0
dt_orig[, cv_risk_measure := ifelse(effect_size_measure=="relative risk", 0, 1)]

# add non-smoker indicator
dt_orig[, cv_non_smoker := ifelse((grepl("non", cc_unexposed_def, ignore.case=T) | grepl("non", cohort_unexp_def, ignore.case=T)), 1, 0)]

# add cv_bmd to indicate whether estimates adjusted for BMD
dt_orig[, cv_bmd := ifelse(confounders_bmd==0, 0, 1)]


# numbers of confounders controlled
# for confounders, replace NA with 0
for (j in c(confounders)){
  set(dt_orig,which(is.na(dt_orig[[j]])),j,0)
}

# count numbers of confounders
dt_orig[, adjustment := apply(.SD, 1, sum), .SDcols=confounders]

# count confounders in confounde_other
dt_orig[, confounders_other := gsub(",", ";", confounders_other)]
# dt_orig[, confounders_other] %>% str_split(., ";", simplify=F) %>% lapply(., length)
dt_orig[!is.na(confounders_other), num := unlist(lapply(str_split(confounders_other, ";"), length))]
dt_orig[is.na(confounders_other), num := 0]

dt_orig[, adjustment:= adjustment + num]


# adding adjustment level covariate
# if age+sex+num2 adjusted, cv_adj=0, 
# if age+sex+num1 adjusted, cv_adj=1, 
# if age+sex adjusted, cv_adj=2,
# if age or sex not adjusted, cv_adj=3

dt_orig[ , cv_adj := as.numeric()]

dt_orig[,adj_age_sex := 0]
dt_orig[confounders_age==1 & confounders_sex==1, adj_age_sex := 1]

# cv_adj=3 if age or sex is not adjusted
dt_orig[adj_age_sex==0, cv_adj := 3]

# cv_adj=2 if only age and sex are adjusted
dt_orig[adj_age_sex==1 & adjustment==2, cv_adj := 2]

# cv_adj=1 if age+sex + 3 more covariates are adjusted
dt_orig[adj_age_sex==1 & adjustment>2 & adjustment<=5, cv_adj := 1]

# cv_adj=0 if age+sex + more than 3 covariates are adjusted
dt_orig[adj_age_sex==1 & adjustment>5, cv_adj := 0]

# check whether there is missing in cv_adj
message("there is ", nrow(dt_orig[is.na(cv_adj)]), " missing values in cv_adj")

# add cascading dummies
dt_orig[, cv_adj_L0 := 1]
dt_orig[, cv_adj_L1 := 1]
dt_orig[, cv_adj_L2 := 1]

# if cv_adj==0, change all dummies to be 0
dt_orig[cv_adj==0, c("cv_adj_L0", "cv_adj_L1", "cv_adj_L2") := 0]
# if cv_adj==1, change cv_adj_L1 and cv_adj_L2 to be 0
dt_orig[cv_adj==1, c("cv_adj_L1", "cv_adj_L2") := 0]
# if cv_adj==2, change cv_adj_L2 to be 0
dt_orig[cv_adj==2, c("cv_adj_L2") := 0]

# remove cv_adj
dt_orig[, cv_adj := NULL]

#-----------------------------------------------------------------
# make the RR consistent, having never smoker/non-smokers as ref
#-----------------------------------------------------------------
dt_orig[,cohort_exposed_def] %>% unique
dt_orig[,cohort_unexp_def] %>% unique
dt_orig[,cc_exposed_def] %>% unique
dt_orig[,cc_unexposed_def] %>% unique
dt_orig[, risk_def] %>% unique

# NOTE: do not have never/non-smoker as the exposed group

#---------------------------
# log-transform the effects
#---------------------------

dt_orig[, effect_size_log := log(mean)]
dt_orig[, se_effect_log := abs((log(upper)-log(lower))/3.92)]
dt_orig[, se_effect_log] %>% summary

# remove missing SE and former smokers
dt_orig <- dt_orig[!is.na(effect_size_log)]
dt_orig <- dt_orig[!is.na(se_effect_log)]

dt_orig <- dt_orig[!(grepl("former", cc_exposed_def, ignore.case=T) | grepl("former", cohort_exposed_def, ignore.case=T)),]

# keep variables that are necessary

cvs <- names(dt_orig)[grepl('cv_', names(dt_orig))]

var_names <- c('extractor','nid', 'effect_size_log', 'se_effect_log', 
               "age_start", "age_end", "adjustment", "percent_male", 
               "most_adj_model", "least_adj_model", confounders, cvs, 
               "cc_exposed_def", "cohort_exposed_def", "new_extraction")

dt_sub <- dt_orig[, var_names, with=F]

# get the max of adjustment by NID and sex. 
dt_sub[, adj_max := max(adjustment), by=c("nid", "cc_exposed_def", "cohort_exposed_def", "age_end", "age_start", "percent_male")]
dt_sub[ ,fully_adj_ind := ifelse(adjustment==adj_max, 1, 0)]

# use the most_adj_model
dt_sub <- dt_sub[most_adj_model==1]

# check the missingness of the data again
summary(dt_sub)

# for confounders, replace NA with 0
for (j in c(confounders, "adjustment")){
  set(dt_sub,which(is.na(dt_sub[[j]])),j,0)
}

# for cv covariates, 1 is bad, so replacing NA to 1 if necessary (very rare)
for (j in c(cvs)){
  set(dt_sub,which(is.na(dt_sub[[j]])),j,1)
}

set(dt_sub, which(is.na(dt_sub[["age_start"]])),"age_start", median(dt_sub[["age_start"]], na.rm = T))
set(dt_sub, which(is.na(dt_sub[["age_end"]])),"age_end", median(dt_sub[["age_end"]], na.rm = T))
dt_sub[is.na(percent_male), percent_male := 0.5]

# check the missingness of the data again
summary(dt_sub)

# change the exposed and unexposed name
dt_sub <- setnames(dt_sub, old = c("effect_size_log", "se_effect_log"), 
                   new = c("ln_effect", "ln_se"))

# create ref and alt
dt_sub[, ref := 0]; dt_sub[, alt := 1]
dt_sub[new_extraction==1, nid] %>% unique %>% length

# save the data
write.csv(dt_sub, file.path(out_dir, paste0(cause, '_all_extraction.csv')), row.names = F)
