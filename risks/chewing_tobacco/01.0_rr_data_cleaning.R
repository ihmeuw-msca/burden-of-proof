## HEADER #################################################################
# Purpose: Cleaning RR Data
#          

## SET-UP #################################################################
# Clear memory
rm(list=ls())

# Runtime configuration ========================================================================
date <- format(Sys.Date(), "%d_%b_%y")

# Load packages, and install if missing ========================================================================
packages <- c("data.table","tidyverse","ggplot2", "readxl", "openxlsx")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Load files and source code =========================================================================
extraction_fp <- paste0("filepath")
file_names <- c("ihd_filename", "stroke_filename", "hnc_filename")
core_fp <- paste0("filepath")
map_fp <- paste0(core_fp, "filepath/")
output_fp <- paste0("filepath")

## SCRIPT ##################################################################
data <- data.table()
for(i in file_names){
  data_temp <- fread(paste0(extraction_fp, i))
  data_temp$review <- str_to_upper(unlist(str_split(i, "_"))[2])
  data <- rbindlist(list(data, data_temp), use.names = T)
  rm(data_temp)
  rm(i)
}

# Dropping all columns with no details
data[data == ""] <- NA
data <- data[, which(unlist(lapply(data, function(x)!all(is.na(x))))), with=F]
data[, dup_NID := .N, by = c("NID", "review")]
if(length(unique(data$dup_NID)) > 1){
  View(data[dup_NID >= 2])
  stop("There is a duplicated row. Check!")
} else {
  print("No duplicates in dataset. Continue!")
  data$dup_NID <- NULL
}
# Fixing a column naming issue
if("ES 8: If this effect size is for a group of the study sample, what is the subgroup?" %in% names(data)){
  setnames(data, "ES 8: If this effect size is for a group of the study sample, what is the subgroup?", "ES 8: If this effect size is for a subgroup of the study sample, what is the subgroup?")
}
names(data) <- str_trim(names(data))

### Changing data format to have one row per effect size, instead of one row per study
es_col_names <- names(data)[which(names(data) %like% "ES " | names(data) %like% "Effect size ")]
study_col_names <- names(data[which(!(names(data) %in% es_col_names))])

es_traits <- str_split(es_col_names, "(:\\s)|Effect size [:alnum:]\\s|Effect size [:alnum:][:alnum:]\\s")
es_traits <- unique(sapply(es_traits, "[", 2)[which(!is.na(sapply(es_traits, "[", 2)))])

es_to_melt <- lapply(es_traits, function(x){if(x %like% "subgroup of the study"){es_col_names[es_col_names %like% x]} else{es_col_names[str_which(es_col_names, paste0(x, "$"))]}})

data <- melt.data.table(data, measure = es_to_melt, value.name = es_traits, variable.name = "Effect size version")

if("Title.1" %in% names(data)){
  data$Title.1 <- NULL
}

## Remove rows that are created with no effect sizes
data <- data[!is.na(`Effect size`)]

## Clean up columns with "Other: " included during consensus
temp <- data %>%
  select_if(~any(.x %like% "Other: ")) %>%
  filter_all(any_vars(. %like% "Other: ")) #%>% 
  # select(!c('Study-level Outcome', 'Outcome', 'Confounders controlled for'))
setDT(temp)
issue_names <- names(temp)
data[, (issue_names) := lapply(.SD, function(x) gsub(pattern = "^Other:", replacement = "", x)), .SDcols = issue_names]
rm(temp)

## Clean up trailing spaces 
columns <- names(data)
data[, (columns) := lapply(.SD, str_trim), .SDcols = columns]

if(F){
  fwrite(data, paste0(extraction_fp, "wide_data_", date, ".csv"))
}

# Clean up missing "percent males"
print(paste0(nrow(data[is.na(`Percentage male`)]), " data sources are missing percentage male details."))
# All of the data 
data[is.na(`Percentage male`), `Percentage male` := 0.5]

# Report some basic summary details
by_review <- copy(data) %>% .[, `ES count` := .N, by = review]
by_review <- unique(by_review[, .(review, NID, `ES count`)])
by_review[, `NID count` := .N, by = review]
by_review <- unique(by_review[, .(review, `NID count`, `ES count`)])
by_review$`Total NIDs`<- length(unique(data$NID))
if(F){
  View(by_review)
  fwrite(by_review, paste0(extraction_fp, "filename", date, ".csv"))
}

## Now lets clean up the exposure and outcome definitions
exposures <- unique(data[, .(NID, `Exposed group definition`, `Exposure temporality`, `Unexposed group definition`, 
                      `Temporality of unexposed group`, `Risk mapping`, `Study-level risk definition`)])
exposures_mapped <- fread(paste0(map_fp, "filename.csv"))
exposures[, NID := as.integer(NID)]
exposures_mapped[exposures_mapped == ""] <- NA
exposures[exposures == ""] <- NA
exp_to_map <- merge(exposures, exposures_mapped, all = T)[is.na(`Effect risk mapping`) | `Effect risk mapping` == ""]
if(nrow(exp_to_map) > 0){
  fwrite(exp_to_map, paste0(map_fp, "filename.csv"))
  stop("You have some exposures that need mapping!")
} else {
  print("Yay! All of your exposures are mapped. Continue.")
}

data[, NID := as.integer(NID)]
data[data == ""] <- NA
data <- merge(data, exposures_mapped, all.x = T)
if(nrow(data[(is.na(`Effect risk mapping`) | `Effect risk mapping` == "")]) > 0){
  stop("You have some exposures that need mapping or there is an issue with merging the correct exposures.")
}  else {
  print("Yay! All of your correct exposures have been merged on!")
}

setDT(data)
outcomes <- unique(data[, .(NID, `Study-level Outcome`, `Study-level outcome definition`, `ICD codes (if available)`, 
                            `Outcome`, `Outcome definition`, `ICD-10 codes`)])
outcomes_mapped <- fread(paste0(map_fp, "filename.csv"))
acause_mapped <- fread(paste0(map_fp, "filename.csv"))
outcomes_mapped[outcomes_mapped == ""] <- NA
outcomes_mapped$cause_id <- NULL
outcomes_mapped$cause <- NULL

if(nrow(outcomes_mapped[is.na(acause)]) > 0){
  outcomes_mapped <- merge(outcomes_mapped, acause_mapped, all.x = T)
  outcomes_mapped$temp <- ""
  for(i in unique(acause_mapped$cause)){
    label <- unique(acause_mapped[cause == i, acause])
    outcomes_mapped[is.na(acause) & Outcome %like% i, temp := ifelse(temp == "", paste0(label), paste0(label, "; ", temp))]
    rm(label)
  }
  outcomes_mapped[is.na(acause), acause := temp]
  outcomes_mapped[is.na(acause)|acause == "" & NID == 512094, acause := "neo_nasal"]
  outcomes_mapped$temp <- NULL
  fwrite(outcomes_mapped, paste0(map_fp, "filename.csv"))
} 

outcomes[, NID := as.integer(NID)]
outcomes_mapped[outcomes_mapped == ""] <- NA
outcomes[outcomes == ""] <- NA
outcomes_to_map <- merge(outcomes, outcomes_mapped, all = T)[is.na(acause)| acause == ""]
if(nrow(outcomes_to_map) > 0){
  fwrite(outcomes_to_map, paste0(map_fp, "filename.csv"))
  stop("You have some outcomes that need mapping!")
} else {
  print("Yay! All of your outcomes are mapped. Continue.")
}
outcomes_mapped[outcomes_mapped == ""] <- NA
data <- merge(data, outcomes_mapped, all.x = T, by = c("NID", "Study-level Outcome", "Outcome", "Outcome definition"))
data[, `Multi-outcome` := ifelse(acause %like% "; ", "Yes", "No")]
if(nrow(data[(is.na(acause) | acause == "")]) > 0){
  stop("You have some outcomes that need mapping or there is an issue with merging the correct outcomes.")
}  else {
  print("Yay! All of your correct outcomes have been merged on!")
}

## Now to clean up covariates
if("If this effect size is for a subgroup of the study sample, what is the subgroup?" %in% names(data)){
  covariates <- unique(data[, .(NID, `Confounders controlled for`, `If this effect size is for a subgroup of the study sample, what is the subgroup?`, `Smoking status`, `Percentage male`)])
} else {
  covariates <- unique(data[, .(NID, `Confounders controlled for`, `Smoking status`, `Percentage male`)])
  covariates$`If this effect size is for a subgroup of the study sample, what is the subgroup?` <- NA
}
covariates_mapped <- fread(paste0(map_fp, "filename.csv"))

covariates[, row_marker := 1:.N]
covariates[, `:=` (standard_covariates = str_split_fixed(`Confounders controlled for`, "Other: ", 2)[1], 
                   `Other covariates` = str_split_fixed(`Confounders controlled for`, "Other: ", 2)[2]), by = row_marker]
covariates$row_marker <- NULL
covariates <- merge(covariates, covariates_mapped, all.x = T)
covariates[covariates == ""] <- NA
if(nrow(covariates[!is.na(`Other covariates`) & is.na(`Covs mapped to`)]) > 0){
  fwrite(covariates[!is.na(`Other covariates`) & is.na(`Covs mapped to`)], paste0(map_fp, "filename.csv"))
  stop("You have some covariates that need mapping!")
} else {
  print("Yay! All of your covariates are mapped. Continue.")
}

covariates[, Confounders := ifelse(!is.na(standard_covariates) & !is.na(`Covs mapped to`), paste0(standard_covariates, `Covs mapped to`), 
                                   ifelse(!is.na(standard_covariates), standard_covariates, `Covs mapped to`))]
covariates$standard_covariates <- NULL
covariates$`Other covariates` <- NULL
covariates$`Covs mapped to` <- NULL
covariates <- unique(covariates)

data <- merge(data, covariates, all.x = T, by = names(covariates)[which(names(covariates) %in% names(data))])
if(nrow(data[(is.na(Confounders) | Confounders == "") & (!is.na(`Confounders controlled for`) | `Confounders controlled for` == "")]) > 0){
  stop("You have some confounders that need mapping or there is an issue with merging the correct confounders")
}  else {
  print("Yay! All of your correct confounders have been merged on!")
}

if(F){
  fwrite(data, paste0(extraction_fp, "filename_", date, ".csv"))
} 

## check if things have been marked for controlled for sex
# We are doing this because the controlling for sex measure is only used to flag data points that are "maximally controlled" and a study that is only among men will still be considered maximally controlled if it controls for age and smoking but does not explicitly control for sex since the entire sample is of one sex.
if(nrow(data[!(Confounders %like% "Sex") & `Percentage male` %in% c(0, 1)]) > 0){
  print("There are some rows that are only Male or only Female and not marked as controlled for Sex. Fixing now!")
  data[!(Confounders %like% "Sex") & `Percentage male` %in% c(0, 1), Confounders := ifelse(is.na(Confounders) | Confounders == "", "Sex", paste0(Confounders, "; Sex"))]
}

## check if things have been marked for controlled for smoking
data[`Smoking status` == "Unknown", `Smoking status` := "Any smoking status"]

# We are doing this because the controlling for smoking measure is only used to flag data points that are "maximally controlled" and a study that is only among smokers will still be considered maximally controlled if it controls for age and sex but does not explicitly control for smoking since the entire sample is of one smoking status.
if(nrow(data[!(Confounders %like% "Smoking") & `Smoking status` != "Any smoking status"])){
  print("There are some rows with specific smoking populations and not marked as controlled for smoking. Fixing now!")
  data[!(Confounders %like% "Smoking") & `Smoking status` != "Any smoking status", Confounders := ifelse(is.na(Confounders) | Confounders == "", "Smoking", paste0(Confounders, "; Smoking"))]
}

##### Making modeling dataset
modelset <- copy(data)
modelset <- modelset[`Effect risk mapping` != "khaini"]

## Turning it long
modelset <- separate_rows(modelset, acause, sep = ";\\s")
setDT(modelset)

##
confounders <- unique(modelset$Confounders)
confounders <- unique(unlist(str_split(confounders, pattern = ";\\s")))
confounders <- confounders[which(!is.na(confounders))]
modelset[, adjustment_count := ifelse(is.na(Confounders), 0, str_count(Confounders, ";\\s")+1)]

modelset[, `:=` (cv_age = ifelse(Confounders %like% "Age", 1, 0), 
                 cv_sex = ifelse(Confounders %like% "Sex", 1, 0), 
                 cv_income = ifelse(Confounders %like% "Income", 1, 0), 
                 cv_smoking = ifelse(Confounders %like% "Smoking", 1, 0), 
                 cv_alcohol = ifelse(Confounders %like% "Alcohol use", 1, 0), 
                 cv_religion = ifelse(Confounders %like% "Religion", 1, 0), 
                 cv_geography = ifelse(Confounders %like% "Geographic region", 1, 0), 
                 cv_other = ifelse(Confounders %like% "Other", 1, 0), 
                 cv_other_smokeless_tob = ifelse(Confounders %like% "Use of other smokeless", 1, 0), 
                 cv_occupation = ifelse(Confounders %like% "Occupation", 1, 0), 
                 cv_oral_hygiene = ifelse(Confounders %like% "Oral hygiene", 1, 0), 
                 cv_race = ifelse(Confounders %like% "Race or ethnicity", 1, 0), 
                 cv_time = ifelse(Confounders %like% "Time", 1, 0), 
                 cv_urban = ifelse(Confounders %like% "Urbanicity", 1, 0), 
                 cv_nontob_chewing = ifelse(Confounders %like% "Use of non-tobacco chewing", 1, 0), 
                 cv_diet = ifelse(Confounders %like% "Dietary components", 1, 0), 
                 cv_genotype = ifelse(Confounders %like% "Genotype", 1, 0), 
                 cv_ses = ifelse(Confounders %like% "Socioeconomic status", 1, 0), 
                 cv_medical_history = ifelse(Confounders %like% "Medical history", 1, 0),
                 cv_language = ifelse(Confounders %like% "Language", 1, 0),
                 cv_shs = ifelse(Confounders %like% "Secondhand smoke exposure", 1, 0),
                 cv_stress = ifelse(Confounders %like% "Stress", 1, 0),
                 cv_aspirin = ifelse(Confounders %like% "Aspirin use", 1, 0),
                 cv_obesity = ifelse(Confounders %like% "Obesity", 1, 0),
                 cv_bmi = ifelse(Confounders %like% "BMI", 1, 0),
                 cv_hypertension = ifelse(Confounders %like% "Hypertension", 1, 0),
                 cv_lpa = ifelse(Confounders %like% "Physical activity", 1, 0),
                 cv_diabetes = ifelse(Confounders %like% "Diabetes", 1, 0))]

# Flipping where the exposed group is non-users
modelset[, `:=` (`Upper CI` = as.numeric(`Upper CI`), `Lower CI` = as.numeric(`Lower CI`), `Effect size` = as.numeric(`Effect size`))]
modelset[`Exposure temporality` == "Never", `:=` (`Upper CI` = 1/`Lower CI`, `Lower CI` = 1/`Upper CI`, `Effect size` = 1/`Effect size`, `Exposure temporality` = `Temporality of unexposed group`, `Temporality of unexposed group` = "Never" )]

## Calculate standard error
if(!("standard_error" %in% names(modelset))){
  modelset[, standard_error := NA_real_]
}

modelset[is.na(standard_error) & `CI type` == "95", `:=` (se_calc = "uses confidence interval")]
modelset[is.na(standard_error) & `CI type` == "90", `:=` (se_calc = "uses confidence interval")]
modelset[is.na(standard_error) & `CI type` == "99", `:=` (se_calc = "uses confidence interval")]

print(paste0("Imputing standard errors for ", nrow(modelset[is.na(standard_error) & is.na(se_calc)]), " rows of data."))

modelset[is.na(standard_error) & is.na(se_calc), `:=` (se_calc = "uses 98th percentile of observed standard errors")]

if(nrow(modelset[is.na(se_calc)]) > 0){
  stop("You have some standard error flags that are still missing!")
} else {
  print("Standard error flags are all filled in!")
}

## Create bias covariates
setnames(modelset, c("Were the study participants geographically representative of the IHME location(s)?", "Were the study participants representative of the study's geographic scope?", "Is the study looking at mortality among people who have already developed the outcome?", "If this effect size is for a subgroup of the study sample, what is the subgroup?"), 
         c("rep_geography", "rep_participants", "rep_prev_disease", "subgroup"))

print("Saving data file now!")
fwrite(modelset, paste0(output_fp, "filename_", date, ".csv"))

