## HEADER #################################################################
# Purpose: Formatting for MR-BRT
#          

## SET-UP #################################################################
# Clear memory
rm(list=ls())

# Runtime configuration ========================================================================
date <- format(Sys.Date(), "%d_%b_%y")

# Load packages, and install if missing ========================================================================
packages <- c("data.table","tidyverse","ggplot2", "readxl", "openxlsx", "msm")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Load files and source code =========================================================================
bundle <- 9089
rei <- 332
cause_ids <- c("neo_esophageal" = 411, "neo_mouth" = 444, "neo_larynx" = 423, "neo_nasopharynx" = 447, "neo_otherpharynx" = 450, 
               "cvd_stroke" = 494, "cvd_ihd" = 493)

data_fp <- paste0("filepath")
cleaned_data_file <- c("ihd_filename.csv", "stroke_filename.csv", "hnc_filename.csv")

## SCRIPT ##################################################################
data <- data.table()
for(i in cleaned_data_file){
  data_temp <- fread(paste0(data_fp, i))
  data <- rbindlist(list(data, data_temp), use.names = T, fill = T)
  rm(data_temp)
  rm(i)
}

########### SET UP BIAS COVARIATE ##################################################
# bc_loc_representative (0 = geographically representative, 1 = not geographically representative)
data[, bc_loc_representative := ifelse(rep_geography == "No", 1, 0)]

# bc_representative (0 = representative of study location, 1 = not representative of study location)
data[, bc_representative := ifelse(rep_participants == "No" | rep_prev_disease == "Yes", 1, 0)]

# bc_exposure_study (0 = multiple exposure measures, 1 = only baseline)
data[, bc_exposure_study := ifelse(`Exposure assessment period` == "Only at baseline", 1, 0)]

# bc_outcome_selfreport (0 = some form of non-self-report outcome ascertainment was used, 1 = self-reported outcome)
data[, bc_outcome_selfreport := ifelse(`Outcome assessment method` %like% "Self-report", 1, 0)]

# bc_controlling
data[, cv_adj_alt := ifelse(cv_age == 0 | cv_sex == 0, 3, 
                                ifelse(cv_age == 1 & cv_sex == 1 & cv_smoking == 0, 2, 
                                       ifelse(cv_age == 1 & cv_sex == 1 & cv_smoking == 1 & adjustment_count == 3, 1, 
                                              ifelse(cv_age == 1 & cv_sex == 1 & cv_smoking == 1 & adjustment_count > 3, 0, NA))))]
## Create cascading dummies: 
data[, bc_adj_L0 := ifelse(cv_adj_alt == 0, 0, 1)]
data[, bc_adj_L1 := ifelse(cv_adj_alt %in% c(0, 1), 0, 1)]
data[, bc_adj_L2 := ifelse(cv_adj_alt %in% c(0, 1, 2), 0, 1)]

### Other covariates
# bc_subpopulation (0 = not a sample subgroup, 1 = subgroup)
data[, bc_subpopulation := ifelse(is.na(subgroup) | subgroup == "", 0, 1)]

# bc_exp_temporality (0 = current exposure, 1 = ever exposure)
data[`Exposure temporality` != "Former", bc_exp_temporality := ifelse(`Exposure temporality` == "Ever", 1, 0)]

# bc_exp_definition (0 = exact, 1 = subtype)
data[, bc_exp_definition := ifelse(`Effect risk mapping` != "exact", 1, 0)]

# bc_aggregate_outcome (0 = exact or subtype outcome, 1 = aggregate outcome)
data[, bc_aggregate_outcome := ifelse(`Multi-outcome` == "Yes", 1, 0)]

############# SELECT MODELING ROWS #######################################################
# Subsetting to current chewing or ever chewing data
total_rows <- nrow(data)
total_nids <- length(unique(data$NID))
data <- data[`Exposure temporality` != "Former" & `Temporality of unexposed group` != "Former"]
print(paste0("Dropping the former chewing effect sizes dropped ", total_rows - nrow(data), " rows of data and ", total_nids - length(unique(data$NID)), " NIDs."))

# check that all rows have been reviewed
if(nrow(data[is.na(to_use_alt) | to_use_alt == ""]) != 0){
  stop("There is a problem where some rows of data have not been reviewed and selected for use!")
}
new_rows <- nrow(data)
new_nids <- length(unique(data$NID))
modelset <- copy(data) %>% .[to_use_alt == 1]
print(paste0("Dropping rows of data that were not selected dropped ", new_rows - nrow(modelset), " rows of data and ", new_nids - length(unique(modelset$NID)), " NIDs."))

############# SET UP FOR MR-BRT ##########################################################
version_id <- 0

modeling_data <- copy(modelset)
modeling_data[, seq := 1:.N, by = "acause"]
modeling_data[, rei_id := rei]
modeling_data[, bundle_id := bundle]
modeling_data[, bundle_version_id := version_id]
modeling_data[, risk_type := "dichotomous"]

for(i in unique(modeling_data$acause)){
  cause <- cause_ids[[i]]
  modeling_data[acause == i, cause_id := cause]
}

setnames(modeling_data, c("NID"), 
         c("study_id"))

# Transform effect sizes and standard error
modeling_data[, ln_rr := log(`Effect size`)]
modeling_data[se_calc == "uses confidence interval", ln_rr_se := (log(`Upper CI`) - log(`Lower CI`))/3.92]
imputed_se <- quantile(modeling_data$ln_rr_se, .98, na.rm = T)
modeling_data[se_calc == "uses 98th percentile of observed standard errors", ln_rr_se := imputed_se]

# Applying downweighting
modeling_data[!is.na(overlapping), observation_count := .N, by = c("acause", "study_id", "overlapping")]
modeling_data[!is.na(overlapping), ln_rr_se := sqrt(observation_count)*ln_rr_se]

# Clean up sample size
modeling_data[`Study design` %like% "ase-", `:=` (`Total cases` = ifelse(is.na(`Total cases`), `Exposed cases`+`Unexposed cases`, `Total cases`), 
                                                  `Total controls` = ifelse(is.na(`Total controls`), `Exposed controls`+`Unexposed controls`, `Total controls`))]
modeling_data[`Study design` %like% "ase-", sample_size := `Total cases`+`Total controls`]
modeling_data[!(`Study design` %like% "ase-"), sample_size := ifelse(is.na(`Total sample size`), `Exposed sample size`+`Unexposed sample size`, `Total sample size`)]

# Renaming columns used for vetting
setnames(modeling_data, c("Study design", "Age start", "Age end", "Age mean", "Age SD", "Effect risk mapping", "Smoking status", "Percentage male", "Outcome definition", "Exposure temporality", "Temporality of unexposed group", "Study ID", "Exposed group definition"), 
         c("study_design", "age_start", "age_end", "age_mean", "age_sd", "risk_def", "smoking_status", "percent_male", "outcome_def", "exp_temp", "unexp_temp", "study_label", "exp_definition_free"))

if(F){ # Testing location-specific models
  modeling_data <- modeling_data %>%
    dplyr::select(`Location 1 IHME Location ID`, seq, se_calc, rei_id, bundle_id, bundle_version_id, risk_type, cause_id, study_id, ln_rr, ln_rr_se, bc_loc_representative, bc_representative, bc_exposure_study, bc_outcome_selfreport, bc_adj_L0, bc_adj_L1, bc_adj_L2, 
                  bc_subpopulation, bc_exp_temporality, bc_exp_definition, bc_aggregate_outcome, `Total cases`, `Total controls`,`Exposed cases`, `Unexposed cases`, `Exposed controls`, `Unexposed controls`,
                  acause, study_design, age_start, age_end, age_mean, age_sd, risk_def, smoking_status, percent_male, outcome_def, exp_temp, unexp_temp, study_label, sample_size, exp_definition_free,
                  cv_age, cv_sex, cv_income, cv_smoking, cv_alcohol, cv_religion, cv_geography, cv_other, cv_other_smokeless_tob, cv_occupation, cv_oral_hygiene, cv_race, cv_time, cv_urban, cv_nontob_chewing, cv_diet, cv_genotype, cv_ses, cv_medical_history, cv_language, cv_shs)
  setnames(modeling_data, "Location 1 IHME Location ID", "ihme_loc_id")
  modeling_data[, country_id := str_sub(ihme_loc_id, start = 1L, end = 3L)]
  
  source("filepath")
  locs <- get_location_metadata(35, release_id = 9)
  modeling_data <- merge(modeling_data, locs[, .(ihme_loc_id, region_name)], by.x = "country_id", by.y = "ihme_loc_id")
  
  # Limiting to only Asian countries
  modeling_data <- modeling_data[region_name %like% "Asia"]
} else {
  modeling_data <- modeling_data %>%
    dplyr::select(seq, se_calc, rei_id, bundle_id, bundle_version_id, risk_type, cause_id, study_id, ln_rr, ln_rr_se, bc_loc_representative, bc_representative, bc_exposure_study, bc_outcome_selfreport, bc_adj_L0, bc_adj_L1, bc_adj_L2, 
                  bc_subpopulation, bc_exp_temporality, bc_exp_definition, bc_aggregate_outcome, `Total cases`, `Total controls`,`Exposed cases`, `Unexposed cases`, `Exposed controls`, `Unexposed controls`,
                  acause, study_design, age_start, age_end, age_mean, age_sd, risk_def, smoking_status, percent_male, outcome_def, exp_temp, unexp_temp, study_label, sample_size, exp_definition_free,
                  cv_age, cv_sex, cv_income, cv_smoking, cv_alcohol, cv_religion, cv_geography, cv_other, cv_other_smokeless_tob, cv_occupation, cv_oral_hygiene, cv_race, cv_time, cv_urban, cv_nontob_chewing, cv_diet, cv_genotype, cv_ses, cv_medical_history, cv_language, cv_shs, reasoning_alt, overlapping)
}

#Sensitivity tests
if(F){ # No aggregate outcomes
  modeling_data <- modeling_data[bc_aggregate_outcome == 0]
  modeling_data$bc_aggregate_outcome <- NULL
}

if(F){ # Threshold for number of exposed individuals
  length(unique(modeling_data$study_id))
  unknown <- modeling_data[is.na(`Exposed cases`) | is.na(`Exposed controls`) | is.na(`Unexposed cases`) | is.na(`Unexposed controls`) & !(study_design %like% "Prospective cohort")]
  issues <- modeling_data[!((`Exposed cases` >= 10 & `Exposed controls` >= 10 & `Unexposed cases` >= 10 & `Unexposed controls` >= 10) | study_design %like% "Prospective cohort")]
  temp <- modeling_data[(`Exposed cases` >= 10 & `Exposed controls` >= 10 & `Unexposed cases` >= 10 & `Unexposed controls` >= 10) | study_design %like% "Prospective cohort"]
  
  modeling_data <- rbindlist(list(unknown, temp), use.names = T)
  modeling_data <- unique(modeling_data)
}

if(F){ # Dropping ever users
  modeling_data <- modeling_data[bc_exp_temporality == 0]
  modeling_data$bc_exp_temporality <- NULL
}

if(F){ # No covariates
  modeling_data <- modeling_data %>%
    dplyr::select(seq, se_calc, rei_id, bundle_id, bundle_version_id, risk_type, cause_id, study_id, ln_rr, ln_rr_se, `Total cases`, `Total controls`,`Exposed cases`, `Unexposed cases`, `Exposed controls`, `Unexposed controls`,
                  acause, study_design, age_start, age_end, age_mean, age_sd, risk_def, smoking_status, percent_male, outcome_def, exp_temp, unexp_temp, study_label, sample_size, exp_definition_free,
                  cv_age, cv_sex, cv_income, cv_smoking, cv_alcohol, cv_religion, cv_geography, cv_other, cv_other_smokeless_tob, cv_occupation, cv_oral_hygiene, cv_race, cv_time, cv_urban, cv_nontob_chewing, cv_diet, cv_genotype, cv_ses, cv_medical_history, cv_language, cv_shs)
}

if(F){ # Only non-smokers
  modeling_data <- modeling_data[smoking_status == "Never smokers" | smoking_status == "Non-smokers"]
}

# Renaming bias covariates
col_names <- colnames(modeling_data)[colnames(modeling_data) %like% "bc_"]
new_col_names <- str_replace(col_names, "^bc_", "cov_")
setnames(modeling_data, col_names, new_col_names)

paste("Saving file now!")
for(i in unique(modelset$acause)){
  cause <- cause_ids[[i]]
  subset <- modeling_data[cause_id == cause]

  all_cvs <- names(subset)[names(subset) %like% 'cov_']
  relevant_cols <- names(colSums(subset[, ..all_cvs])[which(colSums(subset[, ..all_cvs]) <= 1 | colSums(subset[, ..all_cvs]) >= nrow(subset) - 1)])
  subset <- subset %>% select(!all_of(relevant_cols))

  fwrite(subset, paste0(data_fp, "filename-", i, ".csv"))
}

if(F){ # Sex-specific models
  modelset <- copy(data) %>% .[to_use == 1]
  print(paste0("Dropping rows of data that were not selected dropped ", new_rows - nrow(modelset), " rows of data and ", new_nids - length(unique(modelset$NID)), " NIDs."))
  
  version_id <- 0
  
  modeling_data <- copy(modelset)
  modelset[, seq := 1:.N, by = "acause"]
  modelset[, rei_id := rei]
  modelset[, bundle_id := bundle]
  modelset[, bundle_version_id := version_id]
  modelset[, risk_type := "dichotomous"]
  
  for(i in unique(modelset$acause)){
    cause <- cause_ids[[i]]
    modelset[acause == i, cause_id := cause]
  }
  
  setnames(modelset, c("NID"), 
           c("study_id"))
  
  # Transform effect sizes and standard error
  modelset[, ln_rr := log(`Effect size`)]
  modelset[se_calc == "uses confidence interval", ln_rr_se := (log(`Upper CI`) - log(`Lower CI`))/3.92]
  imputed_se <- quantile(modeling_data$ln_rr_se, .98, na.rm = T)
  modelset[se_calc == "uses 98th percentile of observed standard errors", ln_rr_se := imputed_se]
  
  # Applying downweighting
  modelset[!is.na(overlapping), observation_count := .N, by = c("acause", "study_id", "overlapping")]
  modelset[!is.na(overlapping), ln_rr_se := sqrt(observation_count)*ln_rr_se]
  
  # Clean up sample size
  modelset[`Study design` %like% "ase-", `:=` (`Total cases` = ifelse(is.na(`Total cases`), `Exposed cases`+`Unexposed cases`, `Total cases`), 
                                               `Total controls` = ifelse(is.na(`Total controls`), `Exposed controls`+`Unexposed controls`, `Total controls`))]
  modelset[`Study design` %like% "ase-", sample_size := `Total cases`+`Total controls`]
  modelset[!(`Study design` %like% "ase-"), sample_size := ifelse(is.na(`Total sample size`), `Exposed sample size`+`Unexposed sample size`, `Total sample size`)]
  
  # Renaming columns used for vetting
  setnames(modelset, c("Study design", "Age start", "Age end", "Age mean", "Age SD", "Effect risk mapping", "Smoking status", "Percentage male", "Outcome definition", "Exposure temporality", "Temporality of unexposed group", "Study ID", "Exposed group definition"), 
           c("study_design", "age_start", "age_end", "age_mean", "age_sd", "risk_def", "smoking_status", "percent_male", "outcome_def", "exp_temp", "unexp_temp", "study_label", "exp_definition_free"))
  
  # Renaming bias covariates
  col_names <- colnames(modelset)[colnames(modelset) %like% "bc_"]
  new_col_names <- str_replace(col_names, "^bc_", "cov_")
  setnames(modelset, col_names, new_col_names)
  
  male_modelset <- modelset[percent_male == 1]
  female_modelset <- modelset[percent_male == 0]
  
  female_modelset <- female_modelset %>%
    dplyr::select(seq, se_calc, rei_id, bundle_id, bundle_version_id, risk_type, cause_id, study_id, ln_rr, ln_rr_se, cov_loc_representative, cov_representative, cov_exposure_study, cov_outcome_selfreport, cov_adj_L0, cov_adj_L1, cov_adj_L2, 
                  cov_subpopulation, cov_exp_temporality, cov_exp_definition, cov_aggregate_outcome, `Total cases`, `Total controls`,`Exposed cases`, `Unexposed cases`, `Exposed controls`, `Unexposed controls`,
                  acause, study_design, age_start, age_end, age_mean, age_sd, risk_def, smoking_status, percent_male, outcome_def, exp_temp, unexp_temp, study_label, sample_size, exp_definition_free,
                  cv_age, cv_sex, cv_income, cv_smoking, cv_alcohol, cv_religion, cv_geography, cv_other, cv_other_smokeless_tob, cv_occupation, cv_oral_hygiene, cv_race, cv_time, cv_urban, cv_nontob_chewing, cv_diet, cv_genotype, cv_ses, cv_medical_history, cv_language, cv_shs)
  male_modelset <- male_modelset %>%
    dplyr::select(seq, se_calc, rei_id, bundle_id, bundle_version_id, risk_type, cause_id, study_id, ln_rr, ln_rr_se, cov_loc_representative, cov_representative, cov_exposure_study, cov_outcome_selfreport, cov_adj_L0, cov_adj_L1, cov_adj_L2, 
                  cov_subpopulation, cov_exp_temporality, cov_exp_definition, cov_aggregate_outcome, `Total cases`, `Total controls`,`Exposed cases`, `Unexposed cases`, `Exposed controls`, `Unexposed controls`,
                  acause, study_design, age_start, age_end, age_mean, age_sd, risk_def, smoking_status, percent_male, outcome_def, exp_temp, unexp_temp, study_label, sample_size, exp_definition_free,
                  cv_age, cv_sex, cv_income, cv_smoking, cv_alcohol, cv_religion, cv_geography, cv_other, cv_other_smokeless_tob, cv_occupation, cv_oral_hygiene, cv_race, cv_time, cv_urban, cv_nontob_chewing, cv_diet, cv_genotype, cv_ses, cv_medical_history, cv_language, cv_shs)
}

paste("Saving file now!")
for(i in unique(male_modelset$acause)){
  cause <- cause_ids[[i]]
  subset <- male_modelset[cause_id == cause]
  
  all_cvs <- names(subset)[names(subset) %like% 'cov_']
  relevant_cols <- names(colSums(subset[, ..all_cvs])[which(colSums(subset[, ..all_cvs]) <= 1 | colSums(subset[, ..all_cvs]) >= nrow(subset) - 1)])
  subset <- subset %>% select(!all_of(relevant_cols))
  
  fwrite(subset, paste0(data_fp, "filename-", i, ".csv"))
}
paste("Saving file now!")
for(i in unique(female_modelset$acause)){
  cause <- cause_ids[[i]]
  subset <- female_modelset[cause_id == cause]
  
  all_cvs <- names(subset)[names(subset) %like% 'cov_']
  relevant_cols <- names(colSums(subset[, ..all_cvs])[which(colSums(subset[, ..all_cvs]) <= 1 | colSums(subset[, ..all_cvs]) >= nrow(subset) - 1)])
  subset <- subset %>% select(!all_of(relevant_cols))
  
  fwrite(subset, paste0(data_fp, "filename-", i, ".csv"))
}


