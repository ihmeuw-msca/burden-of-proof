## SET-UP #################################################################
# Clear memory
rm(list=ls())

# Load packages, and install if missing ========================================================================
packages <- c("data.table","tidyverse","ggplot2", "googlesheets4")

for(p in packages){

    if(p %in% rownames(installed.packages())==FALSE){
      install.packages(p)
    }
    library(p, character.only = T) 

}

# Set parameters ===================================================================
age_of_exposure <- "childhood" # Options: adulthood, childhood, all
exposure_breadth <- "exact" # Options: exact, inclusive

make_figures <- F # Do you want to make and save the heatmaps of RO pairs to explore? Yes = T
include_riskfactors <- F # Do you want to include risk factors in the models you can run? Yes = T
include_impairments <- F # Do you want to include GBD impairments in the models you can run? Yes = T
include_aggregateoutcomes <- F # Do you want to include aggregate outcome definitions in the data? Yes = T
include_pregnancyrecall <- T # Do you want to include violence during pregnancy recall types? Yes = T
gbd_only <- T # Do you want to run models only for outcomes included in the GBD? Yes = T

dropping_missing_se <- T # Do you want to include observations with missing SEs? Yes = F

description <- "main_models"
readme_description <- "These models are the primary models reported in the study with updated data."

if(age_of_exposure == "childhood"){
  exposures <- c("physical", "psychological", "neglect", "sexual") 
} else {
  exposures <- c("sexual", "physical", "psychological", "ipv")
}

run_sensitivity <- T
sensitivity_analyses <- c("perp_specific", "nonperp_specific", "anyperp_specific", "nopregnancyrecall", "onlypregnancy", "maleonly", "femaleonly", "noadjustment")
cleaned_dataset_fp <- "filepath/data_cleaned.csv"
reference_fp <- "filepath"
previous_fp <- "filepath"
output_fp <- "filepath"

not_mapping_outcomes <- T
date <- format(Sys.Date(), "%m-%d-%y")

# Load files and source code =========================================================================
source("filepath/99.0_helpful_functions.R")
cleaned_dataset <- fread(cleaned_dataset_fp)
exclusions <- read_sheet('url', sheet = "Sheet1")

mapped_outcomes <- fread(paste0(reference_fp, "combined_map.csv"))

## SCRIPT ##################################################################
if(dir.exists(paste0(output_fp, description, "/"))){
  if(!dir.exists(paste0(output_fp, description, "/results/"))){
    dir.create(paste0(output_fp, description, "/results/"))
  }
  if(!dir.exists(paste0(output_fp, description, "/data/"))){
    dir.create(paste0(output_fp, description, "/data/"))
  }
} else {
  dir.create(paste0(output_fp, description, "/"))
  dir.create(paste0(output_fp, description, "/results/"))
  dir.create(paste0(output_fp, description, "/data/"))
}

setDT(exclusions)
exclusions[value_of_interest == "NULL", value_of_interest := NA_character_]
exclusions <- exclusions[rowSums(is.na(exclusions)) != ncol(exclusions),]

# Checking to see if exclusions have changed from last time
old_exclusions <- fread(paste0(reference_fp, "current_exclusions.csv")) %>% .[!is.na(V1)]
if(nrow(exclusions) != nrow(old_exclusions)){
  print("Saving new exclusions")
  temporary <- copy(exclusions)
  temporary[, `:=` (`added by` = NULL, added_by_molly = NULL, changed_by_molly = NULL, notes_on_change = NULL, check = NULL, check_notes = NULL)]
  temporary <- apply(temporary,2,as.character)
  write.csv(old_exclusions, paste0(previous_fp, "exclusions_up_to_", date, ".csv"))
  write.csv(temporary, paste0(reference_fp, "current_exclusions.csv"))
}

cleaned_dataset[cleaned_dataset == ""] <- NA

cleaned_dataset <- cleaned_dataset[!is.na(most_adjusted_mean) | !is.na(raw_mean)]

# selecting the effect size to use
cleaned_dataset[, `:=` (most_adjusted_lower = ifelse(most_adjusted_lower == 0, 0.01, most_adjusted_lower), 
                        raw_lower = ifelse(raw_lower == 0, 0.01, raw_lower))]
cleaned_dataset[!is.na(most_adjusted_mean), `:=` (mean = as.numeric(most_adjusted_mean), lower = as.numeric(most_adjusted_lower), upper = as.numeric(most_adjusted_upper), 
                                                  most_adjusted_other_uncertainty_value = as.numeric(most_adjusted_other_uncertainty_value),
                                                  effect_size = most_adjusted_mean, effect_size_type = "adjusted", ln_rr = log(most_adjusted_mean))]
cleaned_dataset[is.na(most_adjusted_mean), `:=` (mean = as.numeric(raw_mean), lower = as.numeric(raw_lower), upper = as.numeric(raw_upper), 
                                                 raw_other_uncertainty_value = as.numeric(raw_other_uncertainty_value),
                                                 effect_size = raw_mean, effect_size_type = "raw", ln_rr = log(raw_mean))]

# cleaning standard errors
cleaned_dataset[Uncertainty_type == "confidence interval" & (!is.na(lower) & lower > 0) & !is.na(upper), ln_rr_se := ifelse(Confidence_interval_level == "95", (log(upper)-log(lower))/3.92, 
                                                                              ifelse(Confidence_interval_level == "90", (log(upper)-log(lower))/3.29,
                                                                                     ifelse(Confidence_interval_level == "99", (log(upper)-log(lower))/5.15, NA)))]

cleaned_dataset[Uncertainty_type %like% "standard error" & is.na(ln_rr_se), ln_rr_se := ifelse(!is.na(most_adjusted_mean), as.numeric(most_adjusted_other_uncertainty_value)/mean, as.numeric(raw_other_uncertainty_value)/mean)]

# Filling in or dropping missing uncertainty
if(dropping_missing_se == T){
  cleaned_dataset <- cleaned_dataset[!is.na(ln_rr_se)]
} else {
  high_quantile_ln_se <- quantile(cleaned_dataset$ln_rr_se, na.rm = T, probs = 0.98)
  cleaned_dataset[is.na(ln_rr_se), `:=` (se_filled = "missing uncertainty", ln_rr_se = high_quantile_ln_se) ]
}

# cleaning some other variables
cleaned_dataset[, violence_type_combination := ifelse(violence_type_combination == "either emotional abuse or neglect or both, but not physical or sexual abuse. ", "and/or", 
                                                      ifelse(violence_type_combination %like% "sole sexual maltreatment type", "and/or", 
                                                             ifelse(violence_type_combination %like% "with or without emotional abuse or neglect", "not only", 
                                                                    ifelse(violence_type_combination == "" | is.na(violence_type_combination) | violence_type_combination == "other:", "unknown", 
                                                                           ifelse(violence_type_combination == "or (but not both)", "or", violence_type_combination)))))]
cleaned_dataset[is.na(violence_type_combination), violence_type_combination := "unknown"]

# Fixes to perpetrator type for sensitivity analyses
cleaned_dataset[exposure_timing == "Child maltreatment", perpetrator_type := ifelse(((Study_ID %in% c("hovens 2015") & exposure_definition %like% "parental attention") | 
                                                                                       (Study_ID %in% c("brown 1999") & exposure_definition %like% "adult caretaker") | 
                                                                                       (Study_ID %in% c("lemasters 2021") & exposure_definition %like% "caregiver") | 
                                                                                       (Study_ID %in% c("mall 2020") & exposure_definition %like% "someone in my family") |
                                                                                       (Study_ID %in% c("abajobir 2017") & exposure_definition %like% "those who were taking care of a child") | 
                                                                                       (Study_ID %in% c("enns 2006") & exposure_definition %like% "people at home") |
                                                                                       (Study_ID %in% c("rajapakse 2020") & exposure_definition %like% "parents/guardians")) 
                                                                                    & violence_type == "neglect", "caregiver", 
                                                                                    ifelse(((Study_ID %in% c("zhang 2023", "galloeag 2017", "dong 2004") & exposure_definition %like% "domestic violence") | 
                                                                                              (Study_ID %in% c("chapman 2004") & exposure_definition %like% "battered mother") | 
                                                                                              (Study_ID %in% c("rajapakse 2020") & exposure_definition %like% "violence against household members")) 
                                                                                           & (violence_type %like% "witnessing" | violence_type %like% "witnessed"), "household member", perpetrator_type))]
cleaned_dataset[Study_ID == "leung 2002" & perpetrator_type == "partner; former partner; family member or caregiver", perpetrator_type := "partner; former partner"]

# cleaning up perpetrator types into different groups we might want to look at
cleaned_dataset[, mapped_perps := ifelse(perpetrator_type %like% "anyone or not-specified" | perpetrator_type %like% "anyone or not specified" | perpetrator_type %in% c("stranger; known person", "any adult or one 5 years older than respondent", "any adult or person older than yourself", "adult or someone 5 years or more older", "adult or older person", "adult or someone five years or more older", "family member or caregiver; stranger; acquaintance", "adult", "parent or adult"), "anyone/unspecified", 
                                         ifelse(perpetrator_type == "stranger", "stranger", 
                                                ifelse(perpetrator_type %in% c("partner", "partner; former partner", "former partner"), "partner", 
                                                       ifelse(perpetrator_type %in% c("non-partner", "non-partner; any adult or person older than yourself"), "non-partner", 
                                                              ifelse(perpetrator_type %in% c("acquaintance", "coworker", "acquaintance; any adult or person older than yourself", "acquaintance; non-family juvenile", "acquaintance; grown-up, another child or kid"), "acquaintance", 
                                                                     ifelse((perpetrator_type %like% "partner; former partner" & perpetrator_type %like% "family" & !(perpetrator_type %like% "stranger")) | perpetrator_type == "known person", "broad non-stranger", 
                                                                            ifelse((perpetrator_type %like% "partner; former partner" & (perpetrator_type %like% "family" | perpetrator_type %like% "non-partner")) | perpetrator_type %in% c("partner; someone important", "partner; any adult or person older than yourself", "partner; another male family member"), "extremely broad", 
                                                                                   ifelse(perpetrator_type %like% "extrafamilial" | perpetrator_type %in% c("someone other than a parent or caretaker", "non-family member"), "non-family/household member", 
                                                                                          ifelse(perpetrator_type %in% c("mother", "in-laws", "guardian", "father", "relative") | perpetrator_type %like% "household" | perpetrator_type %like% "home" | perpetrator_type %like% "caregiver" | perpetrator_type %like% "parent", "family/household member", NA)))))))))]

if(nrow(cleaned_dataset[is.na(mapped_perps)]) > 0){
  message("Perpetrators left to map!")
}
 
cleaned_dataset[, perp_group_strangers := ifelse(mapped_perps == "stranger", "stranger", 
                                                 ifelse(mapped_perps %in% c("anyone/unspecified", "extremely broad", "non-family/household member", "non-partner"), "anyone/unspecific", 
                                                        ifelse(mapped_perps %in% c("family/household member", "partner", "broad non-stranger", "acquaintance"), "non-stranger", NA)))]
cleaned_dataset[, perp_group_partners := ifelse(mapped_perps == "partner", "partner", 
                                                 ifelse(mapped_perps %in% c("anyone/unspecified", "extremely broad", "non-family/household member", "broad non-stranger"), "anyone/unspecific", 
                                                        ifelse(mapped_perps %in% c("non-partner", "acquaintance", "stranger", "family/household member"), "non-partner", NA)))]
cleaned_dataset[, perp_group_family := ifelse(mapped_perps == "family/household member", "family/household member", 
                                                ifelse(mapped_perps %in% c("anyone/unspecified", "extremely broad", "broad non-stranger", "non-partner", "partner"), "anyone/unspecific", 
                                                       ifelse(mapped_perps %in% c("non-family/household member", "acquaintance", "stranger"), "non-family/household member", NA)))]

# adding and cleaning columns needed for modeling
names(cleaned_dataset) <- gsub("bc_", "cov_", names(cleaned_dataset))
setnames(cleaned_dataset, "Covidence_#", "study_id")
cleaned_dataset[, `:=` (rei_id = 0, bundle_id = 0, bundle_version_id = 0, risk_type = "dichotomous", risk_unit = "exposed", ref_risk_lower = 0, ref_risk_upper = 0, alt_risk_lower = 0, alt_risk_upper = 0)]

# First we want to merge onto mapped outcomes and check that all of the outcomes have been mapped properly
if(not_mapping_outcomes == F){
  cleaned_dataset <- merge(cleaned_dataset, mapped_outcomes, by = c("cause", "outcome_name"), all.x = T)
  if(nrow(cleaned_dataset[is.na(done) | done == ""]) > 0){
    pending_mapping <- cleaned_dataset[is.na(done) | done == ""]
    fwrite(unique(pending_mapping[, .(cause, outcome_name)]), paste0(output_fp, "outcomes_to_map.csv"))
    stop("Go and map the pending causes before continuing!")
  }
}
# Deal with aggregate outcomes
if(include_aggregateoutcomes){
  cleaned_dataset[, cov_includes_other_outcomes := ifelse(outcome %like% ";" | outcome_grouping_level4 %like% ";" | outcome_grouping_level3 %like% ";" | outcome_grouping_level2 %like% ";", 1, 0)]
  cleaned_dataset <- separate_rows(cleaned_dataset, outcome, sep = "; ")
  cleaned_dataset <- separate_rows(cleaned_dataset, outcome_grouping_level4, outcome_grouping_level3, outcome_grouping_level2, sep = "; ")
} else {
  cleaned_dataset[, `:=` (outcome = ifelse(outcome %like% ";", NA, outcome), 
                          outcome_grouping_level4 = ifelse(outcome_grouping_level4 %like% ";", NA, outcome_grouping_level4), 
                          outcome_grouping_level3 = ifelse(outcome_grouping_level3 %like% ";", NA, outcome_grouping_level3), 
                          outcome_grouping_level2 = ifelse(outcome_grouping_level2 %like% ";", NA, outcome_grouping_level2))]
  cleaned_dataset[, cov_includes_other_outcomes := 0]
}

setDT(cleaned_dataset)
cleaned_dataset[(is.na(outcome_grouping_level4) | outcome_grouping_level4 == "") & type %like% "non_gbd", outcome_grouping_level4 := outcome]
cleaned_dataset[(is.na(outcome_grouping_level3) | outcome_grouping_level3 == "") & type %like% "non_gbd", outcome_grouping_level3 := outcome_grouping_level4]
cleaned_dataset[(is.na(outcome_grouping_level2) | outcome_grouping_level2 == "") & type %like% "non_gbd", outcome_grouping_level2 := outcome_grouping_level3]
print(nrow(cleaned_dataset))
cleaned_dataset <- cleaned_dataset[!is.na(outcome_grouping_level2) & outcome_grouping_level2 != ""]
print(nrow(cleaned_dataset))


# Adding a temporary fix for a study we know is incorrectly labeled
cleaned_dataset <- cleaned_dataset[!(Study_design %like% "cross" & exposure_timing == "Adulthood violence exposure" & Study_ID != "fonck 2005")] # We don't include cross-sectional studies in our adulthood models
cleaned_dataset <- cleaned_dataset[!(Study_ID == "thirumalai 2018" & study_id == 134206)] # Dropped because it is a dissertation and, as a result, not peer-reviewed
cleaned_dataset <- cleaned_dataset[!(Study_ID == "sweet 2013" & violence_type == "physical")]

# Removing pregnancy-specific recalls
if(include_pregnancyrecall == F){
  cleaned_dataset <- cleaned_dataset[!(exposure_recall_type %like% "during pregnancy")]
}

# subset by exposure lifetime
if(age_of_exposure == "adulthood"){
  life_data <- cleaned_dataset[exposure_timing == "Adulthood violence exposure"]
  life_data[, cov_perp := ifelse(perp_group_partners %like% "partner", 1, 0)]
  life_data[, cov_exp_inclyouth := ifelse(exposure_age_lower <= 18, 1, 0)]
  life_data[, cov_pregnancy_recall := ifelse(exposure_recall_type %like% "during pregnancy", 1, 0)]
  
} else if(age_of_exposure == "childhood"){
  life_data <- cleaned_dataset[exposure_timing == "Child maltreatment"]
  life_data[, cov_perp := ifelse(perp_group_family %like% "family/household member", 1, 0)]
  if("sexual" %in% exposures){
    life_data[, cov_exp_before_ageabove15:=ifelse(exposure_age_upper > 15, 1, 0)]
    life_data[, cov_exp_before_agebelow15:=ifelse(exposure_age_upper < 15, 1, 0)]    
  } else {
    life_data[, cov_exp_before_agebelow18 := ifelse(exposure_age_upper < 18, 1, 0)]
    
  }
  life_data[, cov_teen_exp := ifelse(exposure_age_lower >= 13, 1, 0)]
  
} else if(age_of_exposure == "all"){
  life_data <- copy(cleaned_dataset)
  life_data[, cov_perp := ifelse(perp_group_strangers %like% "stranger", 1, 0)]
  life_data[, cov_childhood := ifelse(exposure_timing == "Child maltreatment", 1, 0)]
  
} else {
  stop("Don't recognize the lifetime period specified!")
}

# creating covariates
life_data[is.na(cov_sex_uncontrolled), cov_sex_uncontrolled := ifelse((outcome %like% "partum" | outcome %like% "maternal" | outcome %like% "natal") | confounders %like% "sex", 0, 1)]
life_data[, cov_non_lifetime_recall := ifelse(exposure_recall_type %like% "year" | exposure_recall_type %like% "during pregnancy", 1, 0)]
life_data[, cov_L1 := ifelse(cov_age_uncontrolled == 0 & cov_sex_uncontrolled == 0 & num_confounders > 2, 0, 1)]
life_data[, cov_male_included_es := ifelse(sex %like% "Male", 1, 0)]
life_data[, cov_both_sexes_es := ifelse(sex %like% "Combined", 0, 1)]
life_data[, cov_subpop := ifelse(cov_males_only == 1 | cov_females_only == 1, 1, 0)]

# Remove risk-factor data from the dataset if needed
if(include_riskfactors == F){
  life_data <- life_data[!(type %like% "risk")]
} 

# Remove non-GBD causes from the dataset if we don't want to model with them yet
if(gbd_only == T){
  life_data <- life_data[!(type %like% "non_gbd")]
}

# Remove impairments data from the dataset if needed
if(include_impairments == F){
  life_data <- life_data[!(type %like% "impairment")]
} 

# subset by exposure
if(exposure_breadth == "exact"){
  for(i in exposures){
    if(i == "sexual"){
      temp <- life_data[violence_type == "sexual"]
    } else if(i == "physical"){
      temp <- life_data[violence_type == "physical"]
    } else if(i == "psychological"){
      temp <- life_data[violence_type == "psychological" | violence_type == "emotional" | violence_type %in% c("controlling behavior", "psychological; humiliation", "battered mother", " family violence") | (violence_type %like% "witness" & !(violence_type %like% "; "))]
      
      temp[, cov_component_expdef := ifelse(exposure_definition %like% "fear" | exposure_definition %like% "controlling", 1, 0)] 
      
    } else if(i == "pregnancy"){
      temp <- life_data[violence_type == "pregnancy"]
    } else if(i == "neglect") {
      temp <- life_data[violence_type %like% "neglect" & !(violence_type %like% "; ")]
    } else if(i == "ipv"){
      temp <- life_data[(violence_type %like% "physical" | violence_type %like% "sexual") & sex == "Female" & perpetrator_type %in% c("partner", "partner; former partner")]
      temp <- temp [!(Study_ID %in% intersect(unique(temp[violence_type == 'physical']$Study_ID), 
                                 unique(temp[violence_type == 'sexual']$Study_ID)) 
                                  & !violence_type %in% c('physical', 'sexual'))]
      
      temp[, cov_psych_exp := ifelse(violence_type %like% "psychological" | violence_type %like% "emotional" | violence_type %like% "controlling" | violence_type %like% "economic", 1, 0)]
      
      temp[, original_violence_type := violence_type]
      temp[, violence_type := "ipv"]
      
    } else {
      stop("Violence type hasn't been defined!")
    }
    
    temp[, count_by_outcome := length(unique(Study_ID)), by = c("outcome")]
    temp[, count_by_level4 := length(unique(Study_ID)), by = c("outcome_grouping_level4")]
    temp[, count_by_level3 := length(unique(Study_ID)), by = c("outcome_grouping_level3")]
    temp[, count_by_level2 := length(unique(Study_ID)), by = c("outcome_grouping_level2")]
    
    assign(paste0(i, "_data"), temp)
  }
} else if(exposure_breadth == "inclusive"){
  for(i in exposures){
    if(i == "sexual"){
      temp <- life_data[violence_type %like% "sexual"]
    } else if(i == "physical"){
      temp <- life_data[violence_type %like% "physical"]
    } else if(i == "psychological"){
      temp <- life_data[violence_type %like% "psychological" | violence_type %like% "emotional" | violence_type %in% c("controlling behavior", "psychological; humiliation", "battered mother", " family violence") | violence_type %like% "witness"]
      
      temp[, cov_component_expdef := ifelse(exposure_definition %like% "fear" | exposure_definition %like% "controlling", 1, 0)] # all exposure definitions change to 0
      
    } else if(i == "pregnancy"){
      temp <- life_data[violence_type == "pregnancy"]
    } else if(i == "neglect") {
      temp <- life_data[violence_type %like% "neglect"]
    } else {
      stop("Violence type hasn't been defined!")
    }
    
    temp[, count_by_outcome := length(unique(Study_ID)), by = c("outcome")]
    temp[, count_by_level4 := length(unique(Study_ID)), by = c("outcome_grouping_level4")]
    temp[, count_by_level3 := length(unique(Study_ID)), by = c("outcome_grouping_level3")]
    temp[, count_by_level2 := length(unique(Study_ID)), by = c("outcome_grouping_level2")]
    
    assign(paste0(i, "_data"), temp)
  }
  
} else {
  stop("Breadth of the violence exposure has not been defined!")
}

## IDENTIFY OUTCOMES OF INTEREST FOR EACH EXPOSURE
if(make_figures){
  
  total_counts <- data.table()
  for(i in exposures){
    to_plot <- copy(get(paste0(i, "_data")))
    
    if(make_figures){
      to_plot[, count_by_outcome_strangers := length(unique(Study_ID)), by = c("outcome", "perp_group_strangers")]
      to_plot[, count_by_outcome_partners := length(unique(Study_ID)), by = c("outcome", "perp_group_partners")]
      to_plot[, count_by_outcome_family := length(unique(Study_ID)), by = c("outcome", "perp_group_family")]
      
      to_plot[, count_by_level4_strangers := length(unique(Study_ID)), by = c("outcome_grouping_level4", "perp_group_strangers")]
      to_plot[, count_by_level3_strangers := length(unique(Study_ID)), by = c("outcome_grouping_level3", "perp_group_strangers")]
      to_plot[, count_by_level2_strangers := length(unique(Study_ID)), by = c("outcome_grouping_level2", "perp_group_strangers")]
      
      to_plot[, count_by_level4_partners := length(unique(Study_ID)), by = c("outcome_grouping_level4", "perp_group_partners")]
      to_plot[, count_by_level3_partners := length(unique(Study_ID)), by = c("outcome_grouping_level3", "perp_group_partners")]
      to_plot[, count_by_level2_partners := length(unique(Study_ID)), by = c("outcome_grouping_level2", "perp_group_partners")]
      
      to_plot[, count_by_level4_family := length(unique(Study_ID)), by = c("outcome_grouping_level4", "perp_group_family")]
      to_plot[, count_by_level3_family := length(unique(Study_ID)), by = c("outcome_grouping_level3", "perp_group_family")]
      to_plot[, count_by_level2_family := length(unique(Study_ID)), by = c("outcome_grouping_level2", "perp_group_family")]
      
      setorder(to_plot, count_by_outcome)
      
      to_plot_strangers <- unique(to_plot[, .(outcome, outcome_grouping_level4, outcome_grouping_level3, outcome_grouping_level2, 
                                              count_by_outcome_strangers, count_by_level4_strangers, count_by_level3_strangers, count_by_level2_strangers,
                                              perp_group_strangers)])
      to_plot_partners <- unique(to_plot[, .(outcome, outcome_grouping_level4, outcome_grouping_level3, outcome_grouping_level2, 
                                             count_by_outcome_partners, count_by_level4_partners, count_by_level3_partners, count_by_level2_partners,
                                             perp_group_partners)])
      to_plot_family <- unique(to_plot[, .(outcome, outcome_grouping_level4, outcome_grouping_level3, outcome_grouping_level2, 
                                           count_by_outcome_family, count_by_level4_family, count_by_level3_family, count_by_level2_family,
                                           perp_group_family)])
      to_plot_totals <- unique(to_plot[, .(outcome, outcome_grouping_level4, outcome_grouping_level3, outcome_grouping_level2,
                                           count_by_outcome, count_by_level4, count_by_level3, count_by_level2)])
      
      to_plot_totals[, `:=` (perp_group = "overall", violence_type = i)]
      
      level_order_strangers <- c("overall", "anyone/unspecific", "stranger", "non-stranger")
      level_order_partners <- c("overall", "anyone/unspecific", "partner", "non-partner")
      level_order_family <- c("overall", "anyone/unspecific", "family/household member", "non-family/household member")
      
      total_counts <- rbindlist(list(to_plot_totals, total_counts), fill = T)
      
      fwrite(to_plot_strangers, paste0(output_fp, "exact_", age_of_exposure, "_", i, "_stranger_counts.csv"))
      fwrite(to_plot_partners, paste0(output_fp, "exact_", age_of_exposure, "_", i, "_partners_counts.csv"))
      fwrite(to_plot_family, paste0(output_fp, "exact_", age_of_exposure, "_", i, "_family_counts.csv"))
  
      for(w in "partners"){
        for(m in c("main", "4", "3", "2")){
          to_plot <- copy(get(paste0("to_plot_", w)))
          total_dataset <- copy(to_plot_totals)
          assign(paste0(w, "_", i, "_plot_", m), 
                 outcome_heatmaps(total_dataset = total_dataset, plotting_dataset = to_plot, perp_type = w, outcome_level = m, perp_levels = get(paste0("level_order_", w)), viol_type = i))
        }
      }
    }
  }

  pdf(paste0(output_fp, exposure_breadth,"_", age_of_exposure, "_heatmaps.pdf"), width = 14, height = 12)
  for(m in c("main", "4", "3", "2")){
    version <- copy(total_counts)
    print(outcome_heatmaps(total_dataset = version, total_only = T, perp_type = w, outcome_level = m, perp_levels = exposures, viol_type = i))
  }

  for(i in exposures){
    for(w in "partners"){
      print(get(paste0(w, "_", i, "_plot_main")))
      print(get(paste0(w, "_", i, "_plot_4")))
      print(get(paste0(w, "_", i, "_plot_3")))
      print(get(paste0(w, "_", i, "_plot_2")))
    }
  }
  
  dev.off()
}

### Now we want to select the outcomes with >= 3 studies for each exposure
ro_combos <- data.table()
for(i in exposures){
  exp_data <- get(paste0(i, "_data"))
  outcomes_of_interest <- data.table("outcome" = unique(exp_data[count_by_outcome >= 3 & !(is.na(outcome) | outcome == ""), outcome]), "level" = "main")
  aggregate_outcomes_level4 <- data.table("outcome" = unique(exp_data[(count_by_outcome < 3 | is.na(outcome) | outcome == "") & count_by_level4 >= 3 & !(is.na(outcome_grouping_level4) | outcome_grouping_level4 == ""), outcome_grouping_level4]), "level" = "level4")
  aggregate_outcomes_level3 <- data.table("outcome" = unique(exp_data[(count_by_outcome < 3 | is.na(outcome) | outcome == "") & (count_by_level4 < 3 | is.na(outcome_grouping_level4) | outcome_grouping_level4 == "") & count_by_level3 >= 3 & !(is.na(outcome_grouping_level3) | outcome_grouping_level3 == ""), outcome_grouping_level3]), "level" = "level3")
  aggregate_outcomes_level2 <- data.table("outcome" = unique(exp_data[(count_by_outcome < 3 | is.na(outcome) | outcome == "") & (count_by_level4 < 3 | is.na(outcome_grouping_level4) | outcome_grouping_level4 == "") & (count_by_level3 < 3 | is.na(outcome_grouping_level3) | outcome_grouping_level3 == "") & count_by_level2 >= 3 & !(is.na(outcome_grouping_level2) | outcome_grouping_level2 == ""), outcome_grouping_level2]), "level" = "level2")
  
  exposure_outcome_combos <- rbindlist(list(outcomes_of_interest, aggregate_outcomes_level4, aggregate_outcomes_level3, aggregate_outcomes_level2), use.names = T, fill = T)
  exposure_outcome_combos$exposure <- i
  exposure_outcome_combos <- exposure_outcome_combos[outcome != "" & !is.na(outcome)]
  ro_combos <- rbindlist(list(ro_combos, exposure_outcome_combos), fill = T)
}

# check if there are outcomes that are gbd outcomes but not being selected at a granular level and add them
potential_issue_outcomes <- c("depressive_disorders", "diabetes", "alcohol_use_disorder", "drug_use_disorder", "abortion_miscarriage", "anxiety")
for(i in exposures){
  list_of_outcomes <- unique(ro_combos[level == "main" & exposure == i, outcome])
  
  for(w in potential_issue_outcomes){
    if(w %in% c("abortion_miscarriage")){
      level_option <- "level4"
    } else {
      level_option <- "level3"
    }
    list_to_check <- unique(ro_combos[level ==  level_option& exposure == i, outcome])
    
    if(w %in% list_of_outcomes){
      temp <- data.table("outcome" = w, level = level_option, exposure = i)
      ro_combos <- rbindlist(list(ro_combos, temp), fill = T)
    } else if(w %in% c("alcohol_use_disorder", "drug_use_disorders", "abortion_miscarriage")){
      if(w == "alcohol_use_disorder" & !(w %in% list_to_check) & ("alcohol_dependence" %in% list_of_outcomes | "alcohol_abuse" %in% list_of_outcomes)){
        temp <- data.table("outcome" = w, level = level_option, exposure = i)
        ro_combos <- rbindlist(list(ro_combos, temp), fill = T)
      } else if(w == "drug_use_disorder" & !(w %in% list_to_check) & ("amphetamine_dependence" %in% list_of_outcomes | 
                                             "amphetamine_use" %in% list_of_outcomes | 
                                             "cannabis_abuse" %in% list_of_outcomes | 
                                             "cannabis_dependence" %in% list_of_outcomes | 
                                             "cocaine_use" %in% list_of_outcomes | 
                                             "drug_dependence" %in% list_of_outcomes | 
                                             "drug_use" %in% list_of_outcomes | 
                                             "opioid_abuse" %in% list_of_outcomes | 
                                             "opioid_dependence" %in% list_of_outcomes | 
                                             "opioid_use" %in% list_of_outcomes | 
                                             "other_drug_use" %in% list_of_outcomes)){
        temp <- data.table("outcome" = w, level = level_option, exposure = i)
        ro_combos <- rbindlist(list(ro_combos, temp), fill = T)
      } else if(w == "abortion_miscarriage" & !(w %in% list_to_check) & ("abortion" %in% list_of_outcomes & "miscarriage" %in% list_of_outcomes)){
        temp <- data.table("outcome" = w, level = level_option, exposure = i)
        ro_combos <- rbindlist(list(ro_combos, temp), fill = T)
      }
    }
  }
}

ro_combos <- subset(ro_combos, level != "level2")
ro_combos <- ro_combos[!(outcome %in% c("maternal_disorders", "other", "cardiovascular_disorders", "substance_use_disorders", "substance_use"))]
ro_combos <- unique(ro_combos)

# subset by outcome
models_to_run <- c()
for(i in exposures){
  outcomes_of_interest <- ro_combos[exposure == i]
  exp_data <- get(paste0(i, "_data"))
  
  for(m in unique(outcomes_of_interest$outcome)){
    print(paste0("Exposure: ", i, "; Outcome: ", m))
    # resolve exclusions
    relevant_exclusions <- exclusions[lifetime == age_of_exposure & exposure == i & outcome == m]
    if(nrow(relevant_exclusions) > 0){
      full_exclusions <- relevant_exclusions[other == "all", Study_ID]
      relevant_exclusions <- relevant_exclusions[other != "all" | is.na(other)]
    }
    
    for(level in outcomes_of_interest[outcome == m, level]){
      if(level == "main"){
        model_data <- exp_data[outcome == m]
      } else {
        model_data <- exp_data[get(paste0("outcome_grouping_", level)) == m]
      }
  
      # creating cause-specific covariates
      if(m == "depressive_disorders"){
        model_data[, cov_aggregate_outcomedef := ifelse(cause %in% c("depressive disorders, bipolar disorder", "depressive disorders", 'depressive or anxiety disorder',
                                                                     "depressive disorders, anxiety disorders", "major depressive disorder, bipolar disorder"), 1, 0)]
      } else if(m %in% c("abortion_miscarriage", "abortion", "miscarriage")){
        model_data[, cov_abortion := ifelse(outcome == "abortion", 1, 0)]
        
      } else if(m == "anxiety"){
        model_data[, cov_ptsd := ifelse(outcome == 'ptsd', 1, 0)]

      } else if(m == "diabetes"){
        model_data[, cov_outcome_def:=ifelse(outcome == "diabetes", 0, 1)]
        
      } else if(m == "hiv_sti"){
        model_data[, cov_aggregate_outcomedef := ifelse(outcome == "hiv", 1, 0)]
        model_data[is.na(outcome), cov_aggregate_outcomedef := 0]
        
      } else if(m %in% c("substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
        model_data[, cov_outcome_def := ifelse(outcome %like% "dependence", 0, 1)]

      } else if(m == "selfharm"){
        model_data[, cov_poisoning := ifelse(cause == 'self-poisoning', 1, 0)]
        model_data[, cov_dep_confounding := ifelse(Study_ID %in% 'bandara 2022' | Study_ID %in% 'kaslow 2000', 1, 0)]
        
      } 
      
      model_data <- model_data[, `:=` (outcome_grouping_level2 = NULL, outcome_grouping_level3 = NULL, outcome_grouping_level4 = NULL, outcome = NULL, 
                                       count_by_level4 = NULL, count_by_level2 = NULL, count_by_level3 = NULL, count_by_outcome = NULL)]
      model_data <- unique(model_data)
      
      if(nrow(exclusions[lifetime == age_of_exposure & exposure == i & outcome == m]) > 0){
        if(length(full_exclusions) > 0){
          model_data <- model_data[!(Study_ID %in% full_exclusions)]
        }
      }
      
      # saving model elswehere for sensitivity analyses
      assign(paste0(i, "_", m, "_data"), model_data)
      
      # apply manual exclusions
      if(nrow(relevant_exclusions) > 0){
        for(ex in 1:nrow(relevant_exclusions)){
          # print(relevant_exclusions[ex])
          print(paste0("Pre-exclusions: ", nrow(model_data)))
          model_data <- model_data[!(Study_ID == relevant_exclusions[ex]$Study_ID & get(relevant_exclusions[ex]$variable_of_interest) == relevant_exclusions[ex]$value_of_interest)]
          print(paste0("Post-exclusions: ", nrow(model_data)))
        }
      }
      
      #### NOW WE WANT TO ADD IN THE DATA POINT SELECTION ALGORITHM HERE
      model_data[, observations_by_study := .N, by = "Study_ID"]
      no_selection_needed <- model_data[observations_by_study == 1]
      selection_needed <- model_data[observations_by_study != 1]
      
      if(nrow(selection_needed) > 0){
        # Select by exposure type
        if(i != "ipv"){
          selection_needed[, exposure_options := ifelse(violence_type == i & violence_type_combination == "only", 0, 
                                                        ifelse(violence_type == i & violence_type_combination %in% c("not only", "and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 1, 
                                                               ifelse(violence_type == i, 3, 
                                                                      ifelse(violence_type %like% i & violence_type_combination %in% c("not only", "and"), 4, 
                                                                             ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 5, 
                                                                                    ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown"), 6, 
                                                                                           ifelse(violence_type %like% i, 7, 8)))))))]
          if(i == "sexual" & m == "asthma"){
            selection_needed[Study_ID == "coogan 2013" & violence_type_combination == "not only", exposure_options := 0]
          }
        } else {
          selection_needed[, exposure_options := 0]
        }
          selection_needed[, use := min(exposure_options), by = "Study_ID"]
          
          print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
          selection_needed <- selection_needed[use == exposure_options]
          print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
      
        
        # Select by exposure temporality
        if(age_of_exposure == "adulthood"){
          selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_age_lower >= 18, 0, 
                                                     ifelse(exposure_recall_type == "lifetime", 1, 
                                                            ifelse(exposure_recall_type != "during pregnancy", 2, 3)))]
        } else if(age_of_exposure == "childhood") {
          selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime", 0, 1)]
        } else {
          selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_timing == "Adulthood violence exposure", 0, 
                                                     ifelse(exposure_recall_type == "lifetime", 1, 2))]
        }
  
        selection_needed[, use := min(exposure_time), by = "Study_ID"]
        
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        selection_needed <- selection_needed[use == exposure_time]
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        
        # Select by outcome
        if(include_aggregateoutcomes){
          if(m %in% c("depressive_disorders", "hiv_sti")){
            selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 0, 0, 
                                                         ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 1, 1, 
                                                                ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
          } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
            selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 0, 0, 
                                                         ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 1, 1, 
                                                                ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
          } else {
            selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 , 0, 1)]
          }
        } else {
          if(m %in% c("depressive_disorders", "hiv_sti")){
            selection_needed[, outcome_options := ifelse(cov_aggregate_outcomedef == 0, 0, 1)]
          } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
            selection_needed[, outcome_options := ifelse(cov_outcome_def == 0, 0, 1)]
          } else {
            selection_needed[, outcome_options := 0]
          }
        }
        
        selection_needed[, use := min(outcome_options), by = "Study_ID"]
        
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        selection_needed <- selection_needed[use == outcome_options]
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        
        # Select by perpetrator
        if(i != "ipv"){
          # We want the broadest category of perpetrator for this model
          selection_needed[, perp_options := ifelse(mapped_perps == "anyone/unspecified", 0, 
                                                    ifelse(perp_group_strangers == "anyone/unspecific" & perp_group_family == "anyone/unspecific" & perp_group_family == "anyone/unspecific", 1, 2))]
        } else {
          selection_needed[, perp_options := 0]
        }
       
        
        selection_needed[, use := min(perp_options), by = "Study_ID"]
        
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        selection_needed <- selection_needed[use == perp_options]
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        
        # remove sub-groups
        selection_needed[, subgroup_options := ifelse(subgroup_analysis == "no" | is.na(subgroup_analysis), 0, 1)]
        
        selection_needed[, use := min(subgroup_options), by = "Study_ID"]
        
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        selection_needed <- selection_needed[use == subgroup_options]
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        
        # Select by sex
        selection_needed[, sex_options := ifelse(sex == "Combined Male and Female" | is.na(sex), 0, 1)]

        selection_needed[, use := min(sex_options), by = "Study_ID"]
        
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        selection_needed <- selection_needed[use == sex_options]
        print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
        
        if("demakakos 2020" %in% unique(selection_needed$Study_ID) & m == "abortion_miscarriage"){
          selection_needed[Study_ID == "demakakos 2020", subgroup_analysis_free_text := outcome_definition]
        }
        if("coogan 2013" %in% unique(selection_needed$Study_ID) & m == "asthma"){
          selection_needed[Study_ID == "coogan 2013", subgroup_analysis_free_text := violence_type_combination]
        }
        
        ## Inflate the standard errors as needed
        selection_needed[, inflation_factors := .N, by = c("Study_ID", "sex", "subgroup_analysis_free_text", "age_upper", "exposed_level")]
        selection_needed[inflation_factors != 1, ln_rr_se := ln_rr_se*sqrt(inflation_factors)]
        
        # Combine selected and data points that do not need selection
        model_data <- rbindlist(list(no_selection_needed, selection_needed), fill = T)
      }
      
      if(age_of_exposure == "childhood" & i == "psychological" & m %like% "alcohol"){
        if("fenton 2013" %in% unique(model_data$Study_ID) & "laflair 2013" %in% unique(model_data$Study_ID)){
          model_data[Study_ID == "fenton 2013" | Study_ID == "laflair 2013", `:=` (inflation_factors = 2, ln_rr_se = ln_rr_se*sqrt(2))]
        }
      }
      
      if(length(unique(model_data$study_id)) < 3){
        print(paste0("There is insufficient data to run a model for ", i, "-", m))
        next
      } else if(length(unique(model_data$study_id)) >= 3){
        # Removing duplicates that do not meet our criteria
        model_data <- removing_untestable_covariates(model_data, verbose = F)
        
        model_data[, `:=` (cause_id = 123, cause_name = m)]   
        model_data$seq <- 0:(nrow(model_data)-1)
      
        if(age_of_exposure == "childhood"){
          model_data[, exposure_definition := str_replace_all(iconv(exposure_definition, to = "ASCII//TRANSLIT"), "'", "")]
          model_data[, outcome_definition := str_replace_all(iconv(outcome_definition, to = "ASCII//TRANSLIT"), "'", "")]
          model_data[, Exposure_assessment_instrument := str_replace_all(iconv(Exposure_assessment_instrument, to = "ASCII//TRANSLIT"), "'", "")]
          model_data[, Outcome_assessment_method := str_replace_all(iconv(Outcome_assessment_method, to = "ASCII//TRANSLIT"), "'", "")]
        }
        
        write.csv(model_data, paste0(output_fp, description, "/data/", age_of_exposure, "_", i, "-", m, "_", level, ".csv"), row.names = F, fileEncoding = "UTF-8")
        
        models_to_run <- c(models_to_run, paste0(age_of_exposure, "_", i, "-", m, "_", level))
      }
    }
  }
}

models_to_run_text <- paste(models_to_run, collapse = " ")
models_to_run_text <- paste0("dichotomous_pipeline -i ./data -o ./main_models_trim -p ", models_to_run_text)
setwd(paste0(output_fp, description, "/"))
write.table(models_to_run_text, file = "models_to_run.txt", sep = " ", row.names = FALSE, col.names = FALSE)
write.table(models_to_run, file = "plots_to_run.txt", sep = ", ", row.names = FALSE, col.names = FALSE)
write.table(readme_description, file = "README.txt", sep = " ", row.names = FALSE, col.names = FALSE)



#### Sensitivity analyses
if(run_sensitivity == T){
  for(options in sensitivity_analyses){
    models_to_run <- c()
    female_models_to_run <- c()
    male_models_to_run <- c()
    for(i in exposures){
      outcomes_of_interest <- ro_combos[exposure == i]
      for(m in unique(outcomes_of_interest$outcome)){
        print(paste0("Exposure: ", i, "; Outcome: ", m))
        # resolve exclusions
        relevant_exclusions <- exclusions[lifetime == age_of_exposure & exposure == i & outcome == m]
        if(nrow(relevant_exclusions) > 0){
          full_exclusions <- relevant_exclusions[other == "all", Study_ID]
          relevant_exclusions <- relevant_exclusions[other != "all" | is.na(other)]
        }
        
        for(level in outcomes_of_interest[outcome == m, level]){
          if(level == "main"){
            model_data <- exp_data[outcome == m]
          } else {
            model_data <- exp_data[get(paste0("outcome_grouping_", level)) == m]
          }
        }
        
        # get model data for sensitivity analyses
        model_data <- copy(get(paste0(i, "_", m, "_data")))
        
        # apply manual exclusions
        if(nrow(relevant_exclusions) > 0){
          for(ex in 1:nrow(relevant_exclusions)){
            print(paste0("Pre-exclusions: ", nrow(model_data)))
            model_data <- model_data[!(Study_ID == relevant_exclusions[ex]$Study_ID & get(relevant_exclusions[ex]$variable_of_interest) == relevant_exclusions[ex]$value_of_interest)]
            print(paste0("Post-exclusions: ", nrow(model_data)))
          }
        }
      
        if(options == "noadjustment"){
          model_data[, observations_by_study := .N, by = "Study_ID"]
          no_selection_needed <- model_data[observations_by_study == 1]
          selection_needed <- model_data[observations_by_study != 1]
          
          if(nrow(selection_needed) > 0){
            # Select by exposure type
            selection_needed[, exposure_options := ifelse(violence_type == i & violence_type_combination == "only", 0, 
                                                          ifelse(violence_type == i & violence_type_combination %in% c("not only", "and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 1, 
                                                                 ifelse(violence_type == i, 3, 
                                                                        ifelse(violence_type %like% i & violence_type_combination %in% c("not only", "and"), 4, 
                                                                               ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 5, 
                                                                                      ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown"), 6, 
                                                                                             ifelse(violence_type %like% i, 7, 8)))))))]
            selection_needed[, use := min(exposure_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by exposure temporality
            if(age_of_exposure == "adulthood"){
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_age_lower >= 18, 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 
                                                                ifelse(exposure_recall_type != "during pregnancy", 2, 3)))]
            } else if(age_of_exposure == "childhood") {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime", 0, 1)]
            } else {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_timing == "Adulthood violence exposure", 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 2))]
            }
            
            selection_needed[, use := min(exposure_time), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_time]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by outcome
            if(include_aggregateoutcomes){
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else {
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 , 0, 1)]
              }
            } else {
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_aggregate_outcomedef == 0, 0, 1)]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_outcome_def == 0, 0, 1)]
              } else {
                selection_needed[, outcome_options := 0]
              }
            }
            
            selection_needed[, use := min(outcome_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == outcome_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by perpetrator
            # We want the broadest category of perpetrator for this model
            selection_needed[, perp_options := ifelse(mapped_perps == "anyone/unspecified", 0, 
                                                      ifelse(perp_group_strangers == "anyone/unspecific" & perp_group_family == "anyone/unspecific" & perp_group_family == "anyone/unspecific", 1, 2))]
            
            
            selection_needed[, use := min(perp_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == perp_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # remove sub-groups
            selection_needed[, subgroup_options := ifelse(subgroup_analysis == "no" | is.na(subgroup_analysis), 0, 1)]
            
            selection_needed[, use := min(subgroup_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == subgroup_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by sex
            selection_needed[, sex_options := ifelse(sex == "Combined Male and Female" | is.na(sex), 0, 1)]
            
            selection_needed[, use := min(sex_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == sex_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            ## Inflation factors are not used
            selection_needed[, inflation_factors := .N, by = c("Study_ID", "sex", "subgroup_analysis_free_text", "age_upper", "exposed_level")]

            # Combine selected and data points that do not need selection
            model_data <- rbindlist(list(no_selection_needed, selection_needed), fill = T)
          }
          
          if(length(unique(model_data$study_id)) < 3){
            print(paste0("There is insufficient data to run a model for ", i, "-", m))
            next
          } else if(length(unique(model_data$study_id)) >= 3){
            # Removing duplicates that do not meet our criteria
            model_data <- removing_untestable_covariates(model_data, verbose = F)
            
            model_data[, `:=` (cause_id = 123, cause_name = m)]   
            model_data$seq <- 0:(nrow(model_data)-1)
            
            if(age_of_exposure == "childhood"){
              model_data[, exposure_definition := str_replace_all(iconv(exposure_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, outcome_definition := str_replace_all(iconv(outcome_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Exposure_assessment_instrument := str_replace_all(iconv(Exposure_assessment_instrument, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Outcome_assessment_method := str_replace_all(iconv(Outcome_assessment_method, to = "ASCII//TRANSLIT"), "'", "")]
            }
            
            write.csv(model_data, paste0(output_fp, description, "/", options,"_data/", age_of_exposure, "_", i, "-", m, "_", level, ".csv"), row.names = F, fileEncoding = "UTF-8")
            
            models_to_run <- c(models_to_run, paste0(age_of_exposure, "_", i, "-", m, "_", level))
          }
        } else if(options %in% c("maleonly", "femaleonly")){
          # We want to only select sex-specific data points
          if(options == "maleonly"){
            model_data <- model_data[sex == "Male"]
            
          } else {
            model_data <- model_data[sex == "Female"]
          }

          model_data[, observations_by_study := .N, by = c("Study_ID", "sex")]
          no_selection_needed <- model_data[observations_by_study == 1]
          selection_needed <- model_data[observations_by_study != 1]
          
          if(nrow(selection_needed) > 0){
            # Select by exposure type
            selection_needed[, exposure_options := ifelse(violence_type == i & violence_type_combination == "only", 0, 
                                                          ifelse(violence_type == i & violence_type_combination %in% c("not only", "and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 1, 
                                                                 ifelse(violence_type == i, 3, 
                                                                        ifelse(violence_type %like% i & violence_type_combination %in% c("not only", "and"), 4, 
                                                                               ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 5, 
                                                                                      ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown"), 6, 
                                                                                             ifelse(violence_type %like% i, 7, 8)))))))]
            selection_needed[, use := min(exposure_options), by = c("Study_ID", "sex")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by exposure temporality
            if(age_of_exposure == "adulthood"){
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_age_lower >= 18, 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 
                                                                ifelse(exposure_recall_type != "during pregnancy", 2, 3)))]
            } else if(age_of_exposure == "childhood") {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime", 0, 1)]
            } else {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_timing == "Adulthood violence exposure", 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 2))]
            }
            
            selection_needed[, use := min(exposure_time), by = c("Study_ID", "sex")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_time]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by outcome
            if(include_aggregateoutcomes){
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else {
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 , 0, 1)]
              }
            } else {
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_aggregate_outcomedef == 0, 0, 1)]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_outcome_def == 0, 0, 1)]
              } else {
                selection_needed[, outcome_options := 0]
              }
            }
            
            selection_needed[, use := min(outcome_options), by = c("Study_ID", "sex")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == outcome_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by perpetrator
            # We want the broadest category of perpetrator for this model
            selection_needed[, perp_options := ifelse(mapped_perps == "anyone/unspecified", 0, 
                                                      ifelse(perp_group_strangers == "anyone/unspecific" & perp_group_family == "anyone/unspecific" & perp_group_family == "anyone/unspecific", 1, 2))]
            
            
            selection_needed[, use := min(perp_options), by = c("Study_ID", "sex")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == perp_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # remove sub-groups
            selection_needed[, subgroup_options := ifelse(subgroup_analysis == "no" | is.na(subgroup_analysis), 0, 1)]
            
            selection_needed[, use := min(subgroup_options), by = c("Study_ID", "sex")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == subgroup_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            selection_needed[, inflation_factors := .N, by = c("Study_ID", "sex", "subgroup_analysis_free_text", "age_upper", "exposed_level")]
            selection_needed[inflation_factors != 1, ln_rr_se := ln_rr_se*sqrt(inflation_factors)]
            
            # Combine selected and data points that do not need selection
            model_data <- rbindlist(list(no_selection_needed, selection_needed), fill = T)
          }
          
          if(length(unique(model_data$study_id)) < 3){
            print(paste0("There is insufficient data to run a model for ", i, "-", m))
            next
          } else if(length(unique(model_data$study_id)) >= 3){
            # Removing duplicates that do not meet our criteria
            model_data <- removing_untestable_covariates(model_data, verbose = F)
            
            model_data[, `:=` (cause_id = 123, cause_name = m)]   
            model_data$seq <- 0:(nrow(model_data)-1)
            
            if(age_of_exposure == "childhood"){
              model_data[, exposure_definition := str_replace_all(iconv(exposure_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, outcome_definition := str_replace_all(iconv(outcome_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Exposure_assessment_instrument := str_replace_all(iconv(Exposure_assessment_instrument, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Outcome_assessment_method := str_replace_all(iconv(Outcome_assessment_method, to = "ASCII//TRANSLIT"), "'", "")]
            }
            
            write.csv(model_data, paste0(output_fp, description, "/", options,"_data/", age_of_exposure, "_", i, "-", m, "_", level, ".csv"), row.names = F, fileEncoding = "UTF-8")
            
            models_to_run <- c(models_to_run, paste0(age_of_exposure, "_", i, "-", m, "_", level))
          }
          
        } else if(options %in% c("nopregnancyrecall", "onlypregnancy")){
          # We are dropping all pregnancy-related recall here
          if(options == "nopregnancyrecall"){
            model_data <- model_data[!(exposure_recall_type %like% "during pregnancy")]
          } else if(options == "onlypregnancy"){
            model_data <- model_data[exposure_recall_type %like% "during pregnancy"]
          }
          
          
          model_data[, observations_by_study := .N, by = "Study_ID"]
          no_selection_needed <- model_data[observations_by_study == 1]
          selection_needed <- model_data[observations_by_study != 1]
          
          if(nrow(selection_needed) > 0){
            # Select by exposure type
            selection_needed[, exposure_options := ifelse(violence_type == i & violence_type_combination == "only", 0, 
                                                          ifelse(violence_type == i & violence_type_combination %in% c("not only", "and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 1, 
                                                                 ifelse(violence_type == i, 3, 
                                                                        ifelse(violence_type %like% i & violence_type_combination %in% c("not only", "and"), 4, 
                                                                               ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 5, 
                                                                                      ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown"), 6, 
                                                                                             ifelse(violence_type %like% i, 7, 8)))))))]
            selection_needed[, use := min(exposure_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by exposure temporality
            if(age_of_exposure == "adulthood"){
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_age_lower >= 18, 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 2))]
            } else if(age_of_exposure == "childhood") {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime", 0, 1)]
            } else {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_timing == "Adulthood violence exposure", 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 2))]
            }
            
            selection_needed[, use := min(exposure_time), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_time]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by outcome
            if(include_aggregateoutcomes){
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else {
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 , 0, 1)]
              }
            } else {
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_aggregate_outcomedef == 0, 0, 1)]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_outcome_def == 0, 0, 1)]
              } else {
                selection_needed[, outcome_options := 0]
              }
            }
            
            selection_needed[, use := min(outcome_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == outcome_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by perpetrator
            # We want the broadest category of perpetrator for this model
            selection_needed[, perp_options := ifelse(mapped_perps == "anyone/unspecified", 0, 
                                                      ifelse(perp_group_strangers == "anyone/unspecific" & perp_group_family == "anyone/unspecific" & perp_group_family == "anyone/unspecific", 1, 2))]
            
            
            selection_needed[, use := min(perp_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == perp_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # remove sub-groups
            selection_needed[, subgroup_options := ifelse(subgroup_analysis == "no" | is.na(subgroup_analysis), 0, 1)]
            
            selection_needed[, use := min(subgroup_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == subgroup_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by sex
            selection_needed[, sex_options := ifelse(sex == "Combined Male and Female" | is.na(sex), 0, 1)]
            
            selection_needed[, use := min(sex_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == sex_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            selection_needed[, inflation_factors := .N, by = c("Study_ID", "sex", "subgroup_analysis_free_text", "age_upper", "exposed_level")]
            selection_needed[inflation_factors != 1, ln_rr_se := ln_rr_se*sqrt(inflation_factors)]
            
            # Combine selected and data points that do not need selection
            model_data <- rbindlist(list(no_selection_needed, selection_needed), fill = T)
          }
          
          if(length(unique(model_data$study_id)) < 3){
            print(paste0("There is insufficient data to run a model for ", i, "-", m))
            next
          } else if(length(unique(model_data$study_id)) >= 3){
            # Removing duplicates that do not meet our criteria
            model_data <- removing_untestable_covariates(model_data, verbose = F)
            
            model_data[, `:=` (cause_id = 123, cause_name = m)]   
            model_data$seq <- 0:(nrow(model_data)-1)
            
            if(age_of_exposure == "childhood"){
              model_data[, exposure_definition := str_replace_all(iconv(exposure_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, outcome_definition := str_replace_all(iconv(outcome_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Exposure_assessment_instrument := str_replace_all(iconv(Exposure_assessment_instrument, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Outcome_assessment_method := str_replace_all(iconv(Outcome_assessment_method, to = "ASCII//TRANSLIT"), "'", "")]
            }
            
            write.csv(model_data, paste0(output_fp, description, "/", options,"_data/", age_of_exposure, "_", i, "-", m, "_", level, ".csv"), row.names = F, fileEncoding = "UTF-8")
            
            models_to_run <- c(models_to_run, paste0(age_of_exposure, "_", i, "-", m, "_", level))
          }
          
        } else if(options %in% c("perp_specific", "nonperp_specific", "anyperp_specific")){
          if(options == "perp_specific" & age_of_exposure == "childhood"){
            model_data <- model_data[perp_group_family == "family/household member"]
          } else if (options == "perp_specific" & age_of_exposure == "adulthood")(
            model_data <- model_data[perp_group_partners == "partner"]
          ) else if(options == "nonperp_specific" & age_of_exposure == "childhood"){
            model_data <- model_data[perp_group_family == "non-family/household member"]
            
          } else if(options == "nonperp_specific" & age_of_exposure == "adulthood"){
            model_data <- model_data[perp_group_partners == "non-partner"]
            
          } else if(options == "anyperp_specific" & age_of_exposure == "childhood"){
            model_data <- model_data[perp_group_family == "anyone/unspecific"]
            
          } else if(options == "anyperp_specific" & age_of_exposure == "adulthood"){
            model_data <- model_data[perp_group_partners == "anyone/unspecific"]
            
          }
          
          model_data[, observations_by_study := .N, by = c("Study_ID")]
          no_selection_needed <- model_data[observations_by_study == 1]
          selection_needed <- model_data[observations_by_study != 1]
          
          if(nrow(selection_needed) > 0){
            # Select by exposure type
            selection_needed[, exposure_options := ifelse(violence_type == i & violence_type_combination == "only", 0, 
                                                          ifelse(violence_type == i & violence_type_combination %in% c("not only", "and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 1, 
                                                                 ifelse(violence_type == i, 3, 
                                                                        ifelse(violence_type %like% i & violence_type_combination %in% c("not only", "and"), 4, 
                                                                               ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown") & ((confounders %like% "abuse" & !(confounders %like% "substance abuse" | confounders %like% "drug abuse" | confounders %like% "alcohol abuse")) | confounders %like% "violence" | confounders %like% "neglect" | confounders %like% "maltreatment"), 5, 
                                                                                      ifelse(violence_type %like% i & violence_type_combination %in% c("and/or", "unknown"), 6, 
                                                                                             ifelse(violence_type %like% i, 7, 8)))))))]
            selection_needed[, use := min(exposure_options), by = c("Study_ID")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by exposure temporality
            if(age_of_exposure == "adulthood"){
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_age_lower >= 18, 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 
                                                                ifelse(exposure_recall_type != "during pregnancy", 2, 3)))]
            } else if(age_of_exposure == "childhood") {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime", 0, 1)]
            } else {
              selection_needed[, exposure_time := ifelse(exposure_recall_type == "lifetime" & exposure_timing == "Adulthood violence exposure", 0, 
                                                         ifelse(exposure_recall_type == "lifetime", 1, 2))]
            }
            
            selection_needed[, use := min(exposure_time), by = c("Study_ID")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == exposure_time]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by outcome
            if(include_aggregateoutcomes){
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_aggregate_outcomedef == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 0, 0, 
                                                             ifelse(cov_includes_other_outcomes == 0 & cov_outcome_def == 1, 1, 
                                                                    ifelse(cov_includes_other_outcomes == 1, 2, 3)))]
              } else {
                selection_needed[, outcome_options := ifelse(cov_includes_other_outcomes == 0 , 0, 1)]
              }
            } else {
              if(m %in% c("depressive_disorders", "hiv_sti")){
                selection_needed[, outcome_options := ifelse(cov_aggregate_outcomedef == 0, 0, 1)]
              } else if(m %in% c("diabetes", "substance_use_disorder", "drug_use_disorder", "cannabis_use", "cocaine_use", "amphetamine_use", "opioid_use", "alcohol_use_disorder")){
                selection_needed[, outcome_options := ifelse(cov_outcome_def == 0, 0, 1)]
              } else {
                selection_needed[, outcome_options := 0]
              }
            }
            
            selection_needed[, use := min(outcome_options), by = c("Study_ID")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == outcome_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # remove sub-groups
            selection_needed[, subgroup_options := ifelse(subgroup_analysis == "no" | is.na(subgroup_analysis), 0, 1)]
            
            selection_needed[, use := min(subgroup_options), by = c("Study_ID")]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == subgroup_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            # Select by sex
            selection_needed[, sex_options := ifelse(sex == "Combined Male and Female" | is.na(sex), 0, 1)]
            
            selection_needed[, use := min(sex_options), by = "Study_ID"]
            
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            selection_needed <- selection_needed[use == sex_options]
            print(paste0("Study count: ", length(unique(selection_needed$Study_ID)), "; Rows: ", nrow(selection_needed)))
            
            selection_needed[, inflation_factors := .N, by = c("Study_ID", "sex", "subgroup_analysis_free_text", "age_upper", "exposed_level")]
            selection_needed[inflation_factors != 1, ln_rr_se := ln_rr_se*sqrt(inflation_factors)]
            
            # Combine selected and data points that do not need selection
            model_data <- rbindlist(list(no_selection_needed, selection_needed), fill = T)
          }
          
          if(length(unique(model_data$study_id)) < 3){
            print(paste0("There is insufficient data to run a model for ", i, "-", m))
            next
          } else if(length(unique(model_data$study_id)) >= 3){
            # Removing duplicates that do not meet our criteria
            model_data <- removing_untestable_covariates(model_data, verbose = F)
            
            model_data[, `:=` (cause_id = 123, cause_name = m)]   
            model_data$seq <- 0:(nrow(model_data)-1)
            
            if(age_of_exposure == "childhood"){
              model_data[, exposure_definition := str_replace_all(iconv(exposure_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, outcome_definition := str_replace_all(iconv(outcome_definition, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Exposure_assessment_instrument := str_replace_all(iconv(Exposure_assessment_instrument, to = "ASCII//TRANSLIT"), "'", "")]
              model_data[, Outcome_assessment_method := str_replace_all(iconv(Outcome_assessment_method, to = "ASCII//TRANSLIT"), "'", "")]
            }
            
            write.csv(model_data, paste0(output_fp, description, "/", options,"_data/", age_of_exposure, "_", i, "-", m, "_", level, ".csv"), row.names = F, fileEncoding = "UTF-8")
            
            models_to_run <- c(models_to_run, paste0(age_of_exposure, "_", i, "-", m, "_", level))
          }

        } else {
          stop("Sensitivity analysis type is unknown!")
        }
        
          models_to_run_text <- paste(models_to_run, collapse = " ")
          setwd(paste0(output_fp, description, "/"))
          write.table(models_to_run_text, file = paste0(options, "_models_to_run.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

      }
    }
  }
}






















