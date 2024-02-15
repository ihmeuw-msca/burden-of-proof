##########################################################################################
# Title: Alcohol consumption and ischemic heart disease: a Burden of Proof study
# Purpose: Cleaning data (part II)
##########################################################################################

rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(zoo)

date <- Sys.Date()
replicate <- F

model <- "observational" # or MR

# read in 2016 data set
if (T) {
  path <- "filepath"
  old_data <- fread(paste0(path, cause, ".csv"))

  nrow(old_data)
  length(unique(old_data$nid))
  table((old_data$is_outlier))

  length(unique(old_data[is.na(custom_exp_level_upper) & !is.na(custom_exp_level_lower)]$nid))
  old_data[is.na(custom_exp_level_upper), custom_exp_level_upper := custom_exp_level_lower * 1.5]

  nrow(old_data[is.na(custom_exp_level_lower) & is.na(custom_exp_level_upper)])
  length(unique(old_data[is.na(custom_exp_level_lower) & is.na(custom_exp_level_upper)]$nid))
  test <- old_data[is.na(custom_exp_level_lower) & is.na(custom_exp_level_upper)]
  table(test$is_outlier)
  old_data[is.na(custom_exp_level_lower) & is.na(custom_exp_level_upper), is_outlier := 1]

  old_data[nid == 298262, is_outlier := 1]
  old_data[custom_exp_level_upper > 300, is_outlier := 1]

  table(old_data$is_outlier)
  length(unique(old_data[is_outlier == 1]$nid))

  old_data$harmful <- 0
  old_data[effect_size >= 1, harmful := 1]
  table(old_data$harmful, old_data$sex)
  prop.table(table(old_data$harmful))

  old_data[upper >= 1, harmful := 1]
  table(old_data$harmful, old_data$sex)
  prop.table(table(old_data$harmful))

  length(unique(old_data$nid))
  table(old_data$sex)
  length(unique(old_data[sex == "Male"]$nid))
  length(unique(old_data[sex == "Female"]$nid))
  length(unique(old_data[sex == "Both"]$nid))

  table(old_data$former_drinkers)
  length(unique(old_data[former_drinkers == 1]$nid))

  df_orig_sub <- copy(old_data) %>%
    .[, ln_effect := log(effect_size)] %>%
    .[, log_lower := log(lower)] %>%
    .[, log_upper := log(upper)] %>%
    .[, ln_se := (log_upper - log_lower) / 3.92] %>%
    .[ln_se > 0] # %>%


  names(df_orig_sub)
  df_orig_sub$custom_unexp_level_lower <- 0
  df_orig_sub$custom_unexp_level_upper <- 0
  df_orig_sub[sex == "Male", percent_male := 1]
  df_orig_sub[sex == "Female", percent_male := 0]
  df_orig_sub[sex == "Both", percent_male := 0.5]
  df_orig_sub$cv_non_drinker <- df_orig_sub$former_drinkers

  df_orig_sub <- df_orig_sub[, c(
    "nid", "ln_effect", "ln_se", "custom_exp_level_lower", "custom_exp_level_upper", "custom_unexp_level_lower", "custom_unexp_level_upper",
    "adjusted_for",
    "percent_male", "sex",
    "cv_non_drinker",
    "measure"
  )]

  write.csv(df_orig_sub, paste0("filepath"), row.names = F)
}


# read in 2020 data set
if (T) {
  df_orig <- fread("filepath")
  df_orig <- df_orig[custom_alcohol_type == "Any alcoholic beverages"]

  df_orig <- df_orig[outcome == "Ischemic heart disease"]
  df_orig[
    outcome_def == "\"\"Definite\"\" myocardial infarction was indicated by typical chest pain (lasting for 30 minutes or longer) with the appearance of abnormal and persistent Q or QS waves, changes in cardiac enzyme activity, or both. \"\"Suspect\"\" myocardial infarction was indicated by typical chest pain without a positive electrocardiogram or enzyme activity findings. In this report, definite and suspect myocardial infarctions were combined and presented as myocardial infarction because the relation with alcohol intake was similar.",
    outcome_def_map := "mi"
  ]
  df_orig <- df_orig[!(nid %in% c(309660, 309662, 309668, 309701, 298283, 309710, 309728, 437478, 438565, 429454, 428941, 309677) & outcome_def_map != "ihd")]
  df_orig <- df_orig[!(nid %in% c(432186) & outcome_def_map != "ihd without mi")]


  df_orig[nid %in% c(309664) & note_sr == "these rows are Mantel-Haenszel estimates" & custom_exp_level_lower < 10 & upper == 1.1, note_sr := "these rows are multiple logistic regression estimates"]
  df_orig <- df_orig[!(nid %in% c(309664) & note_sr != "these rows are multiple logistic regression estimates")]
  df_orig <- df_orig[!(nid %in% c(437478) & exp_type != "Baseline")]

  df_orig <- df_orig[!(nid %in% c(309723, 309723) & !percent_male %in% c(0, 1))]
  df_orig <- df_orig[!(nid == 309729 & !(age_start == 17 & age_end == 79))]
  df_orig <- df_orig[!(nid == 338008 & !(age_start == 35 & age_end == 79))]

  df_orig[nid == 298281, new_upper := c(rollmean(custom_exp_level_lower, 2), max(custom_exp_level_upper))]
  df_orig[nid == 298281, new_lower := c(rollmean(c(0, custom_exp_level_lower), 2))]
  df_orig[nid == 298281, custom_exp_level_lower := new_lower]
  df_orig[nid == 298281, custom_exp_level_upper := new_upper]

  model_key <- fread("filepath")
  model_key[, c("Notes", "a", "b", "c", "adjusted_for", "2020", "2016")] <- NULL
  model_key <- model_key[to_use_2020 == 1]

  df_orig <- merge(df_orig, model_key, by = c("nid", "confounder_label", "confounders_other"), all.x = T)
  length(unique(df_orig$nid))
  length(unique(df_orig[to_use_2020 == 1]$nid))

  df_orig <- df_orig[to_use_2020 == 1]

  df_orig <- copy(df_orig) %>%
    .[, ln_effect := log(effect_size)] %>%
    .[, log_lower := log(lower)] %>%
    .[, log_upper := log(upper)] %>%
    .[, ln_se := (log_upper - log_lower) / 3.92]

  unique(df_orig[is.na(lower), "nid"])
  nrow(df_orig[is.na(lower)])
  nrow(df_orig[is.na(ln_se)])

  if (T) {
    df_orig_sub$version <- 2016
    df_orig$version <- 2020

    df_orig <- df_orig[custom_alcohol_type == "Any alcoholic beverages"]
    total <- rbind(df_orig_sub, df_orig, fill = T)
    total$flag <- 0
    total[nid %in% c(173897, 298337, 298338), flag := 1]
    total[nid %in% c(298336) | (nid == 324111 & custom_exp_level_lower == 40.00) | (nid == 298335 & is.na(custom_exp_level_lower)), flag := 2]
    total[is.na(nid), nid := 999999]

    total[is.na(confounder_label), confounder_label := adjusted_for]
    total[is.na(confounder_label), confounder_label := "none"]

    sub_sub_sub <- unique(total[, c(
      "nid", "confounder_label", "confounders_other", "adjusted_for",
      "version"
    )])
    sub_sub_sub$indc <- 1
    sub_sub_sub <- dcast(sub_sub_sub, nid + confounder_label + confounders_other + adjusted_for
    ~ version)
  }

  nrow(df_orig)
  length(unique(df_orig$nid))
  table((df_orig$is_outlier))

  if (T) {
    nrow(df_orig)
    length(unique(df_orig$nid))
    nrow(df_orig[custom_most_adj_model == 1])

    nrow(df_orig[percent_male == 1])
    nrow(df_orig[percent_male == 0])
    nrow(df_orig[!percent_male %in% c(0, 1)])

    length(unique(df_orig[percent_male == 1]$nid))
    length(unique(df_orig[percent_male == 0]$nid))
    length(unique(df_orig[!percent_male %in% c(0, 1)]$nid))

    nrow(df_orig[custom_most_adj_model == 1 & percent_male == 1])
    nrow(df_orig[custom_most_adj_model == 1 & percent_male == 0])
    nrow(df_orig[custom_most_adj_model == 1 & !percent_male %in% c(0, 1)])

    table(df_orig$outcome, is.na(df_orig$custom_sick_alc_quitters))
  }

  df_orig$sex <- "Both"
  df_orig[percent_male == 0, sex := "Female"]
  df_orig[percent_male == 1, sex := "Male"]

  df_orig$is_outlier <- 0
  df_orig[custom_exp_level_upper > 300, is_outlier := 1]

  df_orig$harmful <- 0
  df_orig[effect_size >= 1, harmful := 1]
  table(df_orig$harmful, df_orig$sex)
  prop.table(table(df_orig$harmful))

  df_orig[upper >= 1, harmful := 1]
  table(df_orig$harmful, df_orig$sex)
  prop.table(table(df_orig$harmful))
  prop.table(table(df_orig[custom_exp_level_lower > 50]$harmful))

  unique(df_orig$exposed_group)
  unique(df_orig$unexposed_group)
  table(df_orig$exposed_group, df_orig$unexposed_group)

  df_orig[, c("exposed_group", "unexposed_group", "mean_exp", "mean_unexp", "nid")]

  df_orig$cv_non_drinker <- 0
  df_orig[ # exposed_group == "Non-drinkers" |
    unexposed_group %in% c("Non-drinkers", "Non- or light drinkers"), cv_non_drinker := 1
  ]

  table(df_orig$cv_non_drinker)
  length(unique(df_orig[cv_non_drinker == 1]$nid))

  df_orig_sub <- copy(df_orig)

  summary(df_orig_sub$cohort_sample_size_exp)
  summary(df_orig_sub$cohort_sample_size_unexp)
  summary(df_orig_sub$cohort_sample_size_total)

  summary(df_orig_sub$cc_cases)
  summary(df_orig_sub$cc_control)

  df_orig_sub$sample_size <- df_orig_sub$cc_cases + df_orig_sub$cc_control
  df_orig_sub[is.na(sample_size), sample_size := cohort_sample_size_total]
  nrow(df_orig_sub[is.na(sample_size)])

  df_orig_sub <- as.data.frame(df_orig_sub)
  df_orig_sub$id <- 1:nrow(df_orig_sub)

  cols <- names(df_orig_sub)[grep("confounders", names(df_orig_sub))]
  cols <- c("id", cols)

  confounder_master <- df_orig_sub[, cols]
  confounder_master$confounders_other_ind <- NULL
  confounder_master <- confounder_master %>% tidyr::separate(confounders_other, paste0("cov_", 1:50), ";")
  confounder_master[is.na(confounder_master)] <- 0
  confounder_master[, -1][confounder_master[, -1] != 0] <- 1
  confounder_master <- sapply(confounder_master, as.numeric)
  total_adj <- rowSums(confounder_master[, c(2:ncol(confounder_master))])
  confounder_master <- cbind(confounder_master, total_adj)

  df_orig_sub <- merge(df_orig_sub, confounder_master[, c("id", "total_adj")], by = "id", all = T)
  df_orig2 <- copy(df_orig_sub)
  df_orig2 <- as.data.table(df_orig2)

  table(df_orig2$confounder_label)
  df_orig2[confounders_age == 0 & confounders_sex == 0, cv_adjusted := 0]
  df_orig2[confounders_age == "age" | confounders_sex == "sex", cv_adjusted := 1]
  df_orig2[confounders_age == "age" & confounders_sex == "sex", cv_adjusted := 2]
  df_orig2[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking", cv_adjusted := 3]
  df_orig2[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking" & total_adj <= 4, cv_adjusted := 4]
  df_orig2[confounders_age == "age" & confounders_sex == "sex" & confounders_smoking == "smoking" & total_adj > 4, cv_adjusted := 5]

  table(df_orig2$cv_adjusted)
  table(is.na(df_orig2$cv_adjusted))

  df_orig2[cv_adjusted == 0, `:=`(cv_adjusted_0 = 1, cv_adjusted_1 = 1, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
  df_orig2[cv_adjusted == 1, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 1, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
  df_orig2[cv_adjusted == 2, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 1, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
  df_orig2[cv_adjusted == 3, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 1, cv_adjusted_4 = 1)]
  df_orig2[cv_adjusted == 4, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 0, cv_adjusted_4 = 1)]
  df_orig2[cv_adjusted == 5, `:=`(cv_adjusted_0 = 0, cv_adjusted_1 = 0, cv_adjusted_2 = 0, cv_adjusted_3 = 0, cv_adjusted_4 = 0)]

  df_orig2 <- setnames(df_orig2, old = "rep_prevalent_disease", new = "cv_rep_prevalent_disease")
  df_orig2[is.na(cv_exposure_study), cv_exposure_study := 0]
  df_orig2[, cv_exposure_selfreport := ifelse(exp_method_1 == "Self-report (human/environment)" & is.na(exp_method_2), 1, 0)]
  df_orig2[, cv_selection_bias := ifelse(cohort_dropout_rate <= 1 - 0.95, 0,
    ifelse(cohort_dropout_rate <= 1 - 0.85 & cohort_dropout_rate >= 1 - .95, 1,
      ifelse(cohort_dropout_rate > 1 - 0.85, 2, NA)
    )
  )]
  df_orig2$cv_older <- 0
  df_orig2[age_start >= 50, cv_older := 1]

  # rescale so can print out for cv = 0
  df_orig2[, cv_mortality := (cv_mortality - 1) * -1]
  df_orig2[, cv_incidence := (cv_incidence - 1) * -1]
  df_orig2[, cv_older := (cv_older - 1) * -1]
  df_orig2[, cv_sick_quitters := (cv_sick_quitters - 1) * -1]
  df_orig2[, cv_exposure_study := (cv_exposure_study - 1) * -1]
  df_orig2[, cv_exposure_selfreport := (cv_exposure_selfreport - 1) * -1]
  df_orig2[, cv_rep_geography := (cv_rep_geography - 1) * -1]
  df_orig2[, cv_rep_prevalent_disease := (cv_rep_prevalent_disease - 1) * -1]
  df_orig2[, cv_subpopulation := (cv_subpopulation - 1) * -1]
  df_orig2[, cv_outcome_selfreport := (cv_outcome_selfreport - 1) * -1]
  df_orig2[confounders_hypercholesterolemia == "hyperchol", confounders_hypercholesterolemia := 1]
  df_orig2[confounders_hypertension == "hypertens", confounders_hypertension := 1]
  df_orig2 <- setnames(df_orig2,
    old = c("confounders_hypercholesterolemia", "confounders_hypertension"),
    new = c("cv_hypercholesterolemia", "cv_hypertension")
  )
  df_orig2$cv_hypercholesterolemia <- as.numeric(df_orig2$cv_hypercholesterolemia)
  df_orig2$cv_hypertension <- as.numeric(df_orig2$cv_hypertension)
  df_orig2[, cv_hypertension := (cv_hypertension - 1) * -1]
  df_orig2[, cv_hypercholesterolemia := (cv_hypercholesterolemia - 1) * -1]

  cvs <- c(
    "cv_adjusted_0", "cv_adjusted_1", "cv_adjusted_2", "cv_adjusted_3", "cv_adjusted_4", "cv_rep_prevalent_disease", "cv_rep_geography", "cv_subpopulation", "cv_exposure_study", "cv_outcome_selfreport",
    "cv_sick_quitters", "cv_non_drinker", "percent_male", "cv_older", "cv_incidence", "cv_mortality"
  )


  drop_cvs <- c("cv_subpopulation")
  df_orig2$cv_ihd <- 0
  df_orig2[outcome_def_map == "ihd", cv_ihd := 1]
  cvs <- c(cvs, "cv_hypertension", "cv_hypercholesterolemia", "cv_ihd")

  write.csv(df_orig2, paste0("filepath"), row.names = F)
}

# combine 2016 and 2020 data
data_2016 <- fread(paste0("filepath"))
data_2020 <- fread(paste0("filepath"))
final <- rbind(data_2016, data_2020, fill = T)
write.csv(final, paste0("filepath"), row.names = F)