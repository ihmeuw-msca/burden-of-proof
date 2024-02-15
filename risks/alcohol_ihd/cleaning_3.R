##########################################################################################
# Title: Alcohol consumption and ischemic heart disease: a Burden of Proof study
# Purpose: Cleaning data (part III)
##########################################################################################

rm(list = ls())

library(dplyr)
library(mrbrt001, lib.loc = "filepath")
library(data.table)
library(ggplot2)
library(zoo)

cause <- "ihd"
model <- "observational"

# read in 2022 data set
path <- "filepath"
df_orig <- fread(path)
df_orig <- df_orig[custom_alcohol_type == "Any alcoholic beverages"]

df_orig <- df_orig[outcome == "Ischemic heart disease"]

# map outcomes (split between including mi and non-mi)
df_orig[, outcome_def_map := as.character(outcome_def_map)]
df_orig[
  outcome_def == "Myocardial infarction is defined according to the universal definition of MI as rise and/or fall of cardiac troponin with at least above the 99th percentile of the upper reference limit in a clinical setting consistent with myocardial ischemia, and with at least one of the following: symptoms of ischemia, New ST elevation in ECG, new horizontal or down-sloping ST depression in ECG, New left bundle brunch block on ECG, Development of pathological Q waves on ECG, Imaging evidence of new loss of viable myocardium or new regional wall, Identification of an intracoronary thrombus by angiography or autopsy.",
  outcome_def_map := "mi"
]
df_orig[outcome_def == "NA", outcome_def_map := NA]
df_orig[is.na(outcome_def_map), outcome_def_map := "ihd"]

# create log effect sizes and UIs
df_orig <- copy(df_orig) %>%
  .[, ln_effect := log(mean)] %>%
  .[, log_lower := log(lower)] %>%
  .[, log_upper := log(upper)] %>%
  .[, ln_se := (log_upper - log_lower) / 3.92]

df_orig$sex <- "Both"
df_orig[percent_male == 0, sex := "Female"]
df_orig[percent_male == 1, sex := "Male"]

df_orig$is_outlier <- 0
df_orig[exp_level_upper > 300, is_outlier := 1]

# indicator variable if effect is harmful
df_orig$harmful <- 0
df_orig[mean >= 1, harmful := 1]
df_orig[upper >= 1, harmful := 1]

# indicator variable for non-drinkers
df_orig$cv_non_drinker <- 0
df_orig[unexposed_group %in% c("Non-drinkers", "Non- or light drinkers", "No- or former drinkers"), cv_non_drinker := 1]

# create sample size variable
df_orig$sample_size <- df_orig$cc_cases_total + df_orig$cc_controls_total
df_orig[is.na(sample_size), sample_size := cohort_sample_size_total]

df_orig <- as.data.frame(df_orig)
df_orig$id <- 1:nrow(df_orig)

cols <- names(df_orig)[grep("confounders", names(df_orig))]
cols <- c("id", cols)

# create indicator variable for each additional confounder in the model per effect size
confounder_master <- df_orig[, cols]
confounder_master$confounders_other_ind <- NULL
confounder_master <- confounder_master %>% tidyr::separate(confounders_other, paste0("cov_", 1:50), ";")
confounder_master[is.na(confounder_master)] <- 0
confounder_master[, -1][confounder_master[, -1] != 0] <- 1
confounder_master <- sapply(confounder_master, as.numeric)
total_adj <- rowSums(confounder_master[, c(2:ncol(confounder_master))])
confounder_master <- cbind(confounder_master, total_adj)

df_orig <- merge(df_orig, confounder_master[, c("id", "total_adj")], by = "id", all = T)
df_orig <- as.data.table(df_orig)

# create indicator variable for level of control for confounding
df_orig[confounders_age == 0 & confounders_sex == 0, cv_adjusted := 0]
df_orig[confounders_age == "age" | confounders_sex == "sex", cv_adjusted := 1]
df_orig[confounders_age == "age" & confounders_sex == "sex", cv_adjusted := 2]
df_orig[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking", cv_adjusted := 3]
df_orig[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking" & total_adj <= 4, cv_adjusted := 4]
df_orig[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking" & total_adj > 4, cv_adjusted := 5]

# create indicator variables for each level of control for confounding
df_orig[cv_adjusted == 0, `:=`(cv_adjusted_0 = 1, cv_adjusted_1 = 1, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
df_orig[cv_adjusted == 1, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 1, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
df_orig[cv_adjusted == 2, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
df_orig[cv_adjusted == 3, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
df_orig[cv_adjusted == 4, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 0, cv_adjusted_4 = 1)]
df_orig[cv_adjusted == 5, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 0, cv_adjusted_4 = 0)]

# create indicator variables for sample with disease
df_orig <- setnames(df_orig, old = "rep_prevalent_disease", new = "cv_rep_prevalent_disease")
df_orig[is.na(cv_exposure_study), cv_exposure_study := 0]

# age split indicator variable
df_orig$cv_older <- 0
df_orig[age_start >= 50, cv_older := 1]

# rescale every indicator variable so that we can print out for cv = 0
df_orig[, cv_mortality := (cv_mortality - 1) * -1]
df_orig[, cv_incidence := (cv_incidence - 1) * -1]
df_orig[, cv_older := (cv_older - 1) * -1]
df_orig[, cv_sick_quitters := (cv_sick_quitters - 1) * -1]
df_orig[, cv_exposure_study := (cv_exposure_study - 1) * -1]
df_orig[, cv_rep_geography := (cv_rep_geography - 1) * -1]
df_orig[, cv_rep_prevalent_disease := (cv_rep_prevalent_disease - 1) * -1]
df_orig[, cv_subpopulation := (cv_subpopulation - 1) * -1]
df_orig[confounders_bmi == "bmi", confounders_bmi := 1]
df_orig[
  confounders_hypertension == "hypertens" |
    confounders_other %like% "blood pressure" |
    confounders_other %like% "sbp",
  cv_blood_pressure := 1
]
df_orig[confounders_hypercholesterolemia == "hyperchol" |
  confounders_other %like% "cholesterol" |
  confounders_other %like% "LDL", cv_cholesterol := 1]
df_orig[confounders_other %like% "apoAI" |
  confounders_other %like% "apolipoprotein" |
  confounders_other %like% "apolipoproteins" |
  confounders_other %like% "Apo A-I", cv_apolipoprotein := 1]
df_orig[confounders_other %like% "fibrinogen", cv_fibrinogen := 1]
df_orig[confounders_other %like% "adiponectin", cv_adiponectin := 1]

df_orig <- setnames(df_orig,
  old = c("confounders_bmi", "confounders_diabetes"),
  new = c("cv_bmi", "cv_diabetes")
)

# drop redundant columns
df_orig$confounders_hypertension <- NULL
df_orig$confounders_hypercholesterolemia <- NULL
df_orig$cv_diabetes <- NULL

df_orig[is.na(cv_blood_pressure), cv_blood_pressure := 0]
df_orig[is.na(cv_cholesterol), cv_cholesterol := 0]
df_orig[, cv_apolipoprotein := 0]
df_orig[is.na(cv_fibrinogen), cv_fibrinogen := 0]
df_orig[is.na(cv_adiponectin), cv_adiponectin := 0]

cvs <- c(
  "cv_adjusted_0",
  "cv_adjusted_1",
  "cv_adjusted_2",
  "cv_adjusted_3",
  "cv_adjusted_4",
  "cv_rep_prevalent_disease",
  "cv_rep_geography",
  "cv_subpopulation",
  "cv_exposure_study",
  "cv_exposure_selfreport",
  "cv_outcome_selfreport",
  "cv_sick_quitters",
  "cv_non_drinker",
  "percent_male",
  "cv_older",
  "cv_incidence",
  "cv_mortality"
)

if (model == "observational") {
  drop_cvs <- c("cv_adjusted_3", "cv_subpopulation")
  df_orig$cv_ihd <- 1
  df_orig[outcome_def_map == "ihd", cv_ihd := 0]

  # MI
  df_orig$cv_mi <- 0
  df_orig[outcome_def_map != "mi", cv_mi := 1]

  cvs <- c(
    cvs,
    "cv_bmi",
    "cv_blood_pressure",
    "cv_cholesterol",
    "cv_apolipoprotein",
    "cv_fibrinogen",
    "cv_adiponectin",
    "cv_ihd",
    "cv_mi"
  )
} else if (model %in% c("MR_2SLS", "MR_MVMR", "MR_IVW")) {
  drop_cvs <- c("cv_sick_quitters", "cv_adjusted_3", "cv_subpopulation")
  df_orig[nid == 501650, cv_exposure_selfreport := 1]
  df_orig[nid == 428498, cv_exposure_selfreport := 0]
  df_orig[nid == 501522, cv_exposure_selfreport := 1]
  df_orig[nid == 501158 & model == "MR_2SLS", cv_exposure_selfreport := 0]
  df_orig[nid == 501158 & model == "MR_IVW", cv_exposure_selfreport := 1]
  df_orig[nid == 501158 & model == "MR_MVMR", cv_exposure_selfreport := 1]

  df_orig[, cv_exposure_selfreport := ifelse(cv_exposure_selfreport == T, 1, 0)]

  df_orig$cv_ihd <- 1
  df_orig[outcome_def_map == "ihd", cv_ihd := 0]

  # MI
  df_orig$cv_mi <- 0
  df_orig[outcome_def_map != "mi", cv_mi := 1]

  cvs <- c(
    cvs,
    "cv_bmi",
    "cv_blood_pressure",
    "cv_cholesterol",
    "cv_apolipoprotein",
    "cv_fibrinogen",
    "cv_adiponectin",
    "cv_ihd",
    "cv_mi"
  )
}


# subset data and rename columns for burden of proof pipeline

# add indicator column for MR study or non-MR study
df_orig[note_sr %like% "MR study", mr := 1]
df_orig[!(note_sr %like% "MR study"), mr := 0]

# create lower and upper bounds for exposure mean value
df_orig[is.na(exp_level_lower), exp_level_lower := exp_level_value]
df_orig[is.na(exp_level_upper), exp_level_upper := exp_level_value]

# subset to variables relevant for modeling
df_orig <- df_orig[, c("mr", "nid", "ln_effect", "ln_se", "exp_level_lower", "exp_level_upper", "unexp_level_lower", "unexp_level_upper", ..cvs)]
df_orig[, drop_cvs] <- NULL

df_orig <- setnames(df_orig,
  old = c("exp_level_lower", "exp_level_upper", "unexp_level_lower", "unexp_level_upper"),
  new = c("b_0", "b_1", "a_0", "a_1")
)

# subset data to MR or observational studies and save data
final_obs_total <- df_orig[mr == 0]
final_obs_total$mr <- NULL

final_mr_total <- df_orig[mr == 1]
final_mr_total$mr <- NULL

if (model == "observational") {
  # combine 2016, 2020, and 2022 data
  old_data <- fread("filepath")
  final <- rbind(final_obs_total, old_data, fill = T)

  # study-design specific
  nid_study_design <- fread("filepath")
  final <- merge(final, nid_study_design, by = "nid")
  final$cv_design <- 1
  final[design %like% "cohort", cv_design := 0]
  final$design <- NULL

  # data set for all observational data
  write.csv(final, "filepath", row.names = F)

  # subset observational data by sex
  data_females <- final[percent_male == "0", ]
  data_females$percent_male <- NULL

  data_males <- final[percent_male == "1", ]
  data_males$percent_male <- NULL

  # subset observational data by endpoint
  data_incidence <- final[cv_incidence == 0 & cv_mortality == 1, ]
  data_incidence$cv_incidence <- NULL
  data_incidence$cv_mortality <- NULL

  data_mortality <- final[cv_incidence == 1 & cv_mortality == 0, ]
  data_mortality$cv_incidence <- NULL
  data_mortality$cv_mortality <- NULL

  write.csv(data_incidence, "filepath", row.names = F)
  write.csv(data_mortality, "filepath", row.names = F)

  # subset observational data by study design
  data_cohort <- final[cv_design == "0", ]
  data_cohort$cv_design <- NULL

  date_case_control <- final[cv_design == "1", ]
  date_case_control$cv_design <- NULL

  write.csv(data_cohort, "filepath", row.names = F)
  write.csv(date_case_control, "filepathv", row.names = F)

  # sex-endpoint combination
  data_incidence_prelim <- data_incidence[(data_incidence$nid %in% data_incidence[data_incidence$percent_male == 0, ]$nid) & (data_incidence$nid %in% data_incidence[data_incidence$percent_male == 1, ]$nid), ]
  data_females_incidence <- data_incidence_prelim[percent_male == "0", ]
  data_females_incidence$percent_male <- NULL

  data_males_incidence <- data_incidence_prelim[percent_male == "1", ]
  data_males_incidence$percent_male <- NULL

  write.csv(data_females_incidence, "filepath", row.names = F)
  write.csv(data_males_incidence, "filepath", row.names = F)

  data_mortality_prelim <- data_mortality[(data_mortality$nid %in% data_mortality[data_mortality$percent_male == 0, ]$nid) & (data_mortality$nid %in% data_mortality[data_mortality$percent_male == 1, ]$nid), ]
  data_females_mortality <- data_mortality_prelim[percent_male == "0", ]
  data_females_mortality$percent_male <- NULL

  data_males_mortality <- data_mortality_prelim[percent_male == "1", ]
  data_males_mortality$percent_male <- NULL

  write.csv(data_females_mortality, "filepath", row.names = F)
  write.csv(data_males_mortality, "filepath", row.names = F)

  # myocardial infarction only
  data_mi <- final[cv_mi == 0, ]
  write.csv(data_mi, "filepath", row.names = F)
  data_mi_incidence <- data_mi[cv_incidence == 0 & cv_mortality == 1, ]
  write.csv(data_mi_incidence, "filepath", row.names = F)
  data_mi_mortality <- data_mi[cv_incidence == 1 & cv_mortality == 0, ]
  write.csv(data_mi_mortality, "filepath", row.names = F)
}


# MR specific data
if (model %in% c("MR_2SLS", "MR_MVMR", "MR_IVW")) {
  write.csv(final_mr_total, paste0("filepath"), row.names = F)
}

# subset data for cohort and case-control studies to location only assessed in MR studies
MR <- fread("filepath")
location <- fread("filepath")
location <- location[, .(nid, location_name)]
MR_locations <- merge(MR, location, keep.x = T)

cohort <- fread("filepath")
cohort <- merge(cohort, location, keep.x = T)
cohort_location <- cohort[location_name %like% "United Kingdom" | location_name %like% "Republic of Korea" | location_name %like% "China"]

# drop NIDs that are multi-country
cohort_location <- cohort_location[!(nid == "428293" | nid == "431103")]
cohort_location$location_name <- NULL
write.csv(cohort_location, "filepath", row.names = F)

# subset data to conventional estimates of MR
MR <- fread("filepath")
final <- fread("filepath")
nids <- unique(MR$nid)
MR_conventional <- final[nid %in% nids, ]
write.csv(MR_conventional, "filepath", row.names = F)