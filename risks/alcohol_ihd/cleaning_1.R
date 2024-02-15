##########################################################################################
# Title: Alcohol consumption and ischemic heart disease: a Burden of Proof study
# Purpose: Cleaning data (part I)
##########################################################################################

rm(list = ls())

library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
source("filepath")

locs <- get_location_metadata(location_set_id = 22, gbd_round_id = 6)
type_map <- fread("filepath")
unit_map <- fread("filepath")
path <- "filepath"
input <- paste0(path, "input/")

model <- "observational"

# read data
data <- as.data.table(read_xlsx(paste0("filepath"), sheet = "extraction"))

# Clean up columns
data <- unique(data) %>%
  .[!is.na(extractor) & extractor != "UW NetID of person who extracted the data"] %>%
  .[, ":="(mean = as.numeric(mean),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    nid = str_trim(nid))]
data <- data[!is.na(nid), ]

type_map[, ":="(alcohol_type = custom_alcohol_type,
  alcohol_type_new = custom_alcohol_type_new)]
data <- data %>%
  merge(., type_map, by = "alcohol_type", all.x = T) %>%
  .[!is.na(alcohol_type_new), alcohol_type := alcohol_type_new] %>%
  .[, alcohol_type_new := NULL]

# drop all empty columns
data <- data[rowSums(is.na(data)) != ncol(data), ]

# make these columns numeric
changeCols <- c(
  "exp_freq_lower", "exp_freq_upper", "unexp_freq_lower", "unexp_freq_upper",
  "confounders_age", "confounders_sex", "confounders_education", "confounders_income", "confounders_smoking", "confounders_alcohol_use", "confounders_physical_activity",
  "confounders_dietary_components", "confounders_bmi", "confounders_hypertension", "confounders_diabetes", "confounders_hypercholesterolemia",
  "reverse_causation"
)
data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]

data[exp_unit_def == "NA", exp_unit_def := NA]
data[unexp_unit_def == "NA", unexp_unit_def := NA]

## sex-specific areas columns
changeCols <- c(
  "percent_male", "male_exp_level_lower", "male_exp_level_upper", "male_unexp_level_lower",
  "male_unexp_level_upper", "female_exp_level_lower", "female_exp_level_upper",
  "female_unexp_level_lower", "female_unexp_level_upper", "exp_level_lower", "exp_level_upper",
  "unexp_level_lower", "unexp_level_upper"
)
data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]

# fix exposure upper level signs
data <- data %>%
  .[(!is.na(exp_level_lower)) & is.na(exp_level_lower_sign), exp_level_lower_sign := "="] %>%
  .[(!is.na(exp_level_upper)) & is.na(exp_level_upper_sign), exp_level_upper_sign := "="] %>%
  .[(!is.na(male_exp_level_lower)) & is.na(male_exp_level_lower_sign), male_exp_level_lower_sign := "="] %>%
  .[(!is.na(male_exp_level_upper)) & is.na(male_exp_level_upper_sign), male_exp_level_upper_sign := "="] %>%
  .[(!is.na(female_exp_level_lower)) & is.na(female_exp_level_lower_sign), emale_exp_level_lower_sign := "="] %>%
  .[(!is.na(female_exp_level_upper)) & is.na(female_exp_level_upper_sign), female_exp_level_upper_sign := "="] %>%
  .[(!is.na(unexp_level_lower)) & is.na(unexp_level_lower_sign), unexp_level_lower_sign := "="] %>%
  .[(!is.na(unexp_level_upper)) & is.na(unexp_level_upper_sign), unexp_level_upper_sign := "="] %>%
  .[(!is.na(male_unexp_level_lower)) & is.na(male_unexp_level_lower_sign), male_unexp_level_lower_sign := "="] %>%
  .[(!is.na(male_unexp_level_upper)) & is.na(male_unexp_level_upper_sign), male_unexp_level_upper_sign := "="] %>%
  .[(!is.na(female_unexp_level_lower)) & is.na(female_unexp_level_lower_sign), female_unexp_level_lower_sign := "="] %>%
  .[(!is.na(female_unexp_level_upper)) & is.na(female_unexp_level_upper_sign), female_unexp_level_upper_sign := "="]

data[nid == "500504.0", percent_male := 0.5]

# lower and upper bounds
data <- data %>%
  .[!is.na(male_exp_level_lower) & !is.na(female_exp_level_lower) & is.na(exp_level_lower), exp_level_lower := (male_exp_level_lower * percent_male) + (female_exp_level_lower * (1 - percent_male))] %>%
  .[!is.na(male_exp_level_upper) & !is.na(female_exp_level_upper) & is.na(exp_level_upper), exp_level_upper := (male_exp_level_upper * percent_male) + (female_exp_level_upper * (1 - percent_male))] %>%
  .[!is.na(male_unexp_level_lower) & !is.na(female_unexp_level_lower) & is.na(unexp_level_lower), unexp_level_lower := (male_unexp_level_lower * percent_male) + (female_unexp_level_lower * (1 - percent_male))] %>%
  .[!is.na(male_unexp_level_upper) & !is.na(female_unexp_level_upper) & is.na(unexp_level_upper), unexp_level_upper := (male_unexp_level_upper * percent_male) + (female_unexp_level_upper * (1 - percent_male))] %>%
  .[!is.na(male_exp_level_lower_sign) & !is.na(female_exp_level_lower_sign) & is.na(exp_level_lower_sign) & male_exp_level_lower_sign == female_exp_level_lower_sign, exp_level_lower_sign := female_exp_level_lower_sign] %>%
  .[!is.na(male_exp_level_upper_sign) & !is.na(female_exp_level_upper_sign) & is.na(exp_level_upper_sign) & male_exp_level_upper_sign == female_exp_level_upper_sign, exp_level_upper_sign := female_exp_level_upper_sign] %>%
  .[!is.na(male_unexp_level_lower_sign) & !is.na(female_unexp_level_lower_sign) & is.na(unexp_level_lower_sign) & male_unexp_level_lower_sign == female_unexp_level_lower_sign, unexp_level_lower_sign := female_unexp_level_lower_sign] %>%
  .[!is.na(male_unexp_level_upper_sign) & !is.na(female_unexp_level_upper_sign) & is.na(unexp_level_upper_sign) & male_unexp_level_upper_sign == female_unexp_level_upper_sign, unexp_level_upper_sign := female_unexp_level_upper_sign]

# for mean values
data <- data %>%
  .[!is.na(male_exp_level_value) & !is.na(female_exp_level_value) & is.na(exp_level_value), exp_level_lower := (male_exp_level_value * percent_male) + (female_exp_level_value * (1 - percent_male))] %>%
  .[!is.na(male_unexp_level_value) & !is.na(female_unexp_level_value) & is.na(unexp_level_value), unexp_level_lower := (male_unexp_level_value * percent_male) + (female_unexp_level_value * (1 - percent_male))]


# if data only entered for males and percent_male is 1
# lower and upper bounds
data <- data %>%
  .[!is.na(male_exp_level_lower) & is.na(female_exp_level_lower) & is.na(exp_level_lower) & percent_male == 1, exp_level_lower := male_exp_level_lower] %>%
  .[!is.na(male_exp_level_upper) & is.na(female_exp_level_upper) & is.na(exp_level_upper) & percent_male == 1, exp_level_upper := male_exp_level_upper] %>%
  .[!is.na(male_exp_level_lower_sign) & is.na(female_exp_level_lower_sign) & is.na(exp_level_lower_sign) & percent_male == 1, exp_level_lower_sign := male_exp_level_lower_sign] %>%
  .[!is.na(male_exp_level_upper_sign) & is.na(female_exp_level_upper_sign) & is.na(exp_level_upper_sign) & percent_male == 1, exp_level_upper_sign := male_exp_level_upper_sign] %>%
  .[!is.na(male_unexp_level_lower) & is.na(female_unexp_level_lower) & is.na(unexp_level_lower) & percent_male == 1, unexp_level_lower := male_unexp_level_lower] %>%
  .[!is.na(male_unexp_level_upper) & is.na(female_unexp_level_upper) & is.na(unexp_level_upper) & percent_male == 1, unexp_level_upper := male_unexp_level_upper] %>%
  .[!is.na(male_unexp_level_lower_sign) & is.na(female_unexp_level_lower_sign) & is.na(unexp_level_lower_sign) & percent_male == 1, unexp_level_lower_sign := male_unexp_level_lower_sign] %>%
  .[!is.na(male_unexp_level_upper_sign) & is.na(female_unexp_level_upper_sign) & is.na(unexp_level_upper_sign) & percent_male == 1, unexp_level_upper_sign := male_unexp_level_upper_sign]

# for mean values
data <- data %>%
  .[!is.na(male_exp_level_value) & is.na(female_exp_level_value) & is.na(exp_level_value) & percent_male == 1, exp_level_value := male_exp_level_value] %>%
  .[!is.na(male_unexp_level_value) & is.na(female_unexp_level_value) & is.na(unexp_level_value) & percent_male == 1, unexp_level_value := male_unexp_level_value]


# if data is only entered for females and percent_male is 0
# lower and upper bounds
data <- data %>%
  .[is.na(male_exp_level_lower) & !is.na(female_exp_level_lower) & is.na(exp_level_lower) & percent_male == 0, exp_level_lower := female_exp_level_lower] %>%
  .[is.na(male_exp_level_upper) & !is.na(female_exp_level_upper) & is.na(exp_level_upper) & percent_male == 0, exp_level_upper := female_exp_level_upper] %>%
  .[is.na(male_exp_level_lower_sign) & !is.na(female_exp_level_lower_sign) & is.na(exp_level_lower_sign) & percent_male == 0, exp_level_lower_sign := female_exp_level_lower_sign] %>%
  .[is.na(male_exp_level_upper_sign) & !is.na(female_exp_level_upper_sign) & is.na(exp_level_upper_sign) & percent_male == 0, exp_level_upper_sign := female_exp_level_upper_sign] %>%
  .[is.na(male_unexp_level_lower) & !is.na(female_unexp_level_lower) & is.na(unexp_level_lower) & percent_male == 0, unexp_level_lower := female_unexp_level_lower] %>%
  .[is.na(male_unexp_level_upper) & !is.na(female_unexp_level_upper) & is.na(unexp_level_upper) & percent_male == 0, unexp_level_upper := female_unexp_level_upper] %>%
  .[is.na(male_unexp_level_lower_sign) & !is.na(female_unexp_level_lower_sign) & is.na(unexp_level_lower_sign) & percent_male == 0, unexp_level_lower_sign := female_unexp_level_lower_sign] %>%
  .[is.na(male_unexp_level_upper_sign) & !is.na(female_unexp_level_upper_sign) & is.na(unexp_level_upper_sign) & percent_male == 0, unexp_level_upper_sign := female_unexp_level_upper_sign]

# for mean values
data <- data %>%
  .[is.na(male_exp_level_value) & !is.na(female_exp_level_value) & is.na(exp_level_value) & percent_male == 0, exp_level_value := female_exp_level_value] %>%
  .[is.na(male_unexp_level_value) & !is.na(female_unexp_level_value) & is.na(unexp_level_value) & percent_male == 0, unexp_level_value := female_unexp_level_value]


# clean definitions
data[, exp_unit := str_trim(exp_unit), by = exp_unit]
data[, exp_unit_def := str_trim(exp_unit_def), by = exp_unit_def]

## map units
unit_map[, c("num", "denom", "notes", "beer", "wine", "spirits")] <- NULL
unit_map <- unique(unit_map)

# add multiplier to convert to gram per frequency measure
data <- merge(data, unit_map, by.x = c("exp_unit", "exp_unit_def"), by.y = c("given_unit", "def"), all.x = T)
temp <- data[, .(extractor, exp_unit, multiplier)]

data[is.na(multiplier) & exp_unit == "ml soju/day", multiplier := 0.19731]
data[is.na(multiplier) & exp_unit == "drinks/day", multiplier := 10]
data[is.na(multiplier) & exp_unit == "units/week", multiplier := 13 / 7]
data[is.na(multiplier) & exp_unit == "drinks/year", multiplier := 8 / 365]
data[is.na(multiplier) & exp_unit == "drinks/month", multiplier := 8 / 30]
data[is.na(multiplier) & exp_unit == "drinks/week", multiplier := 8 / 7]

test <- unique(data[is.na(multiplier), c("nid", "outcome", "exp_unit", "exp_unit_def", "extractor")])
test <- unique(data[is.na(multiplier), c("exp_unit", "exp_unit_def")])
test <- test[!(exp_unit %in% c(
  "DU (drink units) years:  (DU per week x 52 x years of drinking)", "drinks per day * years of drinking",
  "g-ethanol/kg-body weight/day", "g/day-years", "ml-years"
))]

data <- data[!is.na(multiplier)]
data <- setnames(data, old = c("multiplier"), new = c("exp_multiplier"))

data <- merge(data, unit_map, by.x = c("unexp_unit", "unexp_unit_def"), by.y = c("given_unit", "def"), all.x = T)
data[is.na(multiplier) & unexp_unit == "drinks/day", multiplier := 10]
data[is.na(multiplier) & unexp_unit == "units/week", multiplier := 13 / 7]
data[is.na(multiplier) & unexp_unit == "drinks/week", multiplier := 11 / 7]
data[is.na(multiplier) & unexp_unit == "glass/day", multiplier := 10 / 1]
data[is.na(multiplier) & unexp_unit == "shot/day", multiplier := 10 / 1]
data[is.na(multiplier) & unexp_unit == "bottle/day", multiplier := 10 / 1]
data[is.na(multiplier) & unexp_unit == "l/week", multiplier := (789.24 * 0.05) / 7]

data <- setnames(data, old = c("multiplier"), new = c("unexp_multiplier"))

## exposed
data[!is.na(exp_freq_lower) & is.na(exp_freq_upper) & exp_freq_unit %in% c("times/week"), exp_freq_upper := 7]
data[!is.na(exp_freq_lower) & is.na(exp_freq_upper) & exp_freq_unit %in% c("days/month heavy episodic drinking"), exp_freq_upper := 30]

data[exp_freq_unit %in% c("days/week", "times/week"), exp_freq_denom := 7]
data[exp_freq_unit %in% c("days/month", "days/month heavy episodic drinking"), exp_freq_denom := 30]
data[exp_freq_unit %in% c("days/year"), exp_freq_denom := 365]

## unexposed
data[!is.na(unexp_freq_lower) & is.na(unexp_freq_upper) & unexp_freq_unit %in% c("times/week"), unexp_freq_upper := 7]

data[unexp_freq_unit %in% c("days/week", "times/week"), unexp_freq_denom := 7]
data[unexp_freq_unit %in% c("days/month", "days/month heavy episodic drinking"), unexp_freq_denom := 30]
data[unexp_freq_unit %in% c("days/year"), unexp_freq_denom := 365]

data[, freq_adjustment := ((unexp_freq_lower + unexp_freq_upper) / 2) / unexp_freq_denom]

data[!is.na(freq_adjustment), unexp_multiplier := unexp_multiplier * freq_adjustment]

# convert to gram per day
changeCols <- c("exp_level_lower", "exp_level_upper", "exp_level_value", "upper", "lower", "mean")
data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]
data[, c("exp_level_lower", "exp_level_upper", "exp_level_value") :=
  list(
    exp_level_lower * exp_multiplier,
    exp_level_upper * exp_multiplier,
    exp_level_value * exp_multiplier
  )]
data[, exp_map := "grams_day"]

changeCols <- c("unexp_level_lower", "unexp_level_upper", "unexp_level_value")
data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]
data[, c("unexp_level_lower", "unexp_level_upper", "level_rr") :=
  list(
    unexp_level_lower * unexp_multiplier,
    unexp_level_upper * unexp_multiplier,
    unexp_level_value * unexp_multiplier
  )]
data[, unexp_map := "grams_day"]

# fill in the lower values when missing
data[!is.na(exp_level_upper) & is.na(exp_level_lower), exp_level_lower := 0]
data[!is.na(unexp_level_upper) & is.na(unexp_level_lower), unexp_level_lower := 0]

# fill in the upper values when missing
data[exp_level_upper >= 499 | is.na(exp_level_upper), exp_level_upper := exp_level_lower * 1.5]
data[exp_level_upper > 500, exp_level_upper := 500]

data[unexp_level_upper >= 499 | is.na(unexp_level_upper), unexp_level_upper := unexp_level_lower * 1.5]
data[unexp_level_upper > 500, unexp_level_upper := 500]

# calculate mean exposure
data[is.na(exp_level_lower) & is.na(exp_level_upper) & !is.na(exp_level_value), mean_exp := exp_level_value]
data[is.na(unexp_level_lower) & is.na(unexp_level_upper) & !is.na(unexp_level_value), mean_unexp := unexp_level_value]

# exposure ranges
data[!is.na(exp_level_lower) & !is.na(exp_level_upper), mean_exp := (exp_level_upper + exp_level_lower) / 2]
data[!is.na(unexp_level_lower) & !is.na(unexp_level_upper), mean_unexp := (unexp_level_upper + unexp_level_lower) / 2]


# generate standard errors for effect size
data[!is.na(lower) & !is.na(upper), se_effect := (upper - lower) / (2 * 1.96)]
impute_se_effect <- quantile(x = data$se_effect, probs = 0.95, na.rm = TRUE)[[1]]
data[is.na(se_effect), se_effect := impute_se_effect]


# create recall variable
data[, exp_recall_period_value := as.numeric(exp_recall_period_value)]

if (T) {
  data[exp_recall_period_value == 1 & exp_recall_period == "Period: months", recall := "1_month"]
  data[exp_recall_period_value == 30 & exp_recall_period == "Period: days", recall := "1_month"]
  data[exp_recall_period_value == 12 & exp_recall_period == "Period: months", recall := "1_year"]
  data[exp_recall_period_value == 24 & exp_recall_period == "Period: months", recall := "2_years"]

  data[exp_recall_period_value == 1 & exp_recall_period == "Period: years", recall := "1_year"]
  data[exp_recall_period_value == 2 & exp_recall_period == "Period: years", recall := "2_years"]
  data[exp_recall_period_value == 10 & exp_recall_period == "Period: years", recall := "10_years"]
  data[exp_recall_period_value == 20 & exp_recall_period == "Period: years", recall := "20_years"]

  data[exp_recall_period_value == 3 & exp_recall_period == "Period: months", recall := "3_months"]
  data[exp_recall_period_value == 6 & exp_recall_period == "Period: months", recall := "6_months"]
  data[exp_recall_period_value == 5 & exp_recall_period == "Period: years", recall := "5_years"]

  data[exp_recall_period_value == 1 & exp_recall_period == "Period: days", recall := "1_days"]
  data[exp_recall_period_value == 7 & exp_recall_period == "Period: days", recall := "1_weeks"]
  data[exp_recall_period_value == 1 & exp_recall_period == "Period: weeks", recall := "1_weeks"]
  data[exp_recall_period_value == 14 & exp_recall_period == "Period: days", recall := "2_weeks"]
  data[exp_recall_period_value == 2 & exp_recall_period == "Period: weeks", recall := "2_weeks"]

  data[exp_recall_period == "Lifetime", recall := "lifetime"]
  data[exp_recall_period_other == "Several months", recall := "3_months"]
  data[exp_recall_period_other == "Point", recall := "point"]
  data[exp_recall_period_other == "point", recall := "point"]
  data[exp_recall_period_other == "consumption 5 years prior to baseline", recall := "other"]
  data[exp_recall_period_other == "consumption at age 21", recall := "other"]
  data[exp_recall_period_other == "Consumption in 1990's", recall := "other"]
  data[exp_recall_period_other == "Recorded at least 2 years before index date randomly generated for each control or date of diagnosis/reported tumor for cases", recall := "other"]
  data[exp_recall_period_other == "subjects were asked about their usual intake of beer, wine, and liquor from the age at which they\r\nstarted drinking at least one alcoholic beverage per month until\r\n1 year before the interview", recall := "other"]
  data[exp_recall_period_other == "not specified", recall := "not_specified"]
  data[exp_recall_period_other == "point recall", recall := "point"]
  data[exp_recall_period_other == "I year prior to cancer diagnosis or hospital admission", recall := "other"]
  data[exp_recall_period_other == "among participants who reported having had at least four alcoholic beverages in any 1 year, we elicited separate age-defined periods in which alcohol consumption (e.g., frequency of drinking, number of drinks) was relatively unchanged", recall := "other"]
  data[exp_recall_period_other == "drinking habits before getting sick", recall := "other"]
  data[exp_recall_period_other == "Consumption at 10 years prior to baseline (includes subjects who were drinkers then but quit before baseline)", recall := "other"]

  data[exp_recall_period_other == "not specified", recall := "not_specified"]
  data[exp_recall_period_other == "unspecified", recall := "not_specified"]
  data[exp_recall_period_other == "Not defined, probably point", recall := "not_specified"]
  data[exp_recall_period_other == "not specified, probably point", recall := "not_specified"]
}

# map exposure definitions
if (T) {
  def_map <- fread("filepath")
  def_map <- unique(def_map)

  data[, ":="(exp_def = exp_group_def,
    unexp_def = unexp_reference_group_def)]

  data <- merge(data, def_map, by.x = "exp_def", by.y = "given", all.x = T)
  data <- setnames(data, old = "mapped", new = "exposed_group")

  data <- merge(data, def_map, by.x = "unexp_def", by.y = "given", all.x = T)
  data <- setnames(data, old = "mapped", new = "unexposed_group")

  data[exp_def == "Middle drinkers", exposed_group := "Current drinkers"]
  data[exp_def == "Current unhealthy alcohol use", exposed_group := "Heavy drinkers"]
  data[exp_def == "Light drinkers or moderate drinkers", exposed_group := "Current drinkers"]
  data[exp_def == "Medium-heavy drinkers", exposed_group := "Current drinkers"]
  data[unexp_def == "Irregular drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Light or non-drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Very light drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Rare or non-drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Rare or never-drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Non-drinkers and occasional drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Low, non-weekly, or never-drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Current healthy alcohol use or non-drinkers", unexposed_group := "Non- or light drinkers"]
  data[unexp_def == "Abstainers & former users", unexposed_group := "No- or former drinkers"]
  data[exp_def == "Binge drank alcohol in past 30 days", exposed_group := "Binge drinkers"]
  data[unexp_def == "Did not binge drink in the past 30 days", unexposed_group := "Non- binge drinkers"]
  data[exp_def == "Former drinkers", exposed_group := "Former drinkers"]
  data[exp_def == "Current drinkers (age <= 63.5 y/o; excluding heavy users)", exposed_group := "Current drinkers"]
  data[exp_def == "Current drinkers (age > 63.5 y/o; excluding heavy users)", exposed_group := "Current drinkers"]
}

# map outcome definitions
outcome_def <- fread("filepath")
data <- merge(data, outcome_def, by = c("outcome", "outcome_def"), all.x = T)

data[(exposed_group %in% c("Never drinkers", "Non- or light drinkers")) |
  (exposed_group %in% c("Non-drinkers") & unexposed_group != "Never drinkers") |
  (mean_unexp > mean_exp), `:=`(
  mean = 1 / mean, lower = 1 / upper, upper = 1 / lower,
  exp_level_lower = unexp_level_lower,
  exp_level_upper = unexp_level_upper,
  unexp_level_lower = exp_level_lower,
  unexp_level_upper = exp_level_upper,
  exp_level_value = unexp_level_value,
  unexp_level_value = exp_level_value,
  exposed_group = unexposed_group,
  unexposed_group = exposed_group,
  mean_exp = mean_unexp,
  mean_unexp = mean_exp
)]

# make confounder indicators numeric
changeCols <- c(
  "confounders_age", "confounders_sex", "confounders_education", "confounders_income", "confounders_smoking", "confounders_alcohol_use", "confounders_physical_activity",
  "confounders_dietary_components", "confounders_bmi", "confounders_hypertension", "confounders_diabetes", "confounders_hypercholesterolemia"
)

data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]

# create indicator variables: control for confounding
if (T) {
  data$cv_adjusted_0 <- 0
  data$cv_adjusted_1 <- 0
  data$cv_adjusted_2 <- 0
  data[confounders_sex != 1 | confounders_age != 1, cv_adjusted_0 := 1]
  data[confounders_sex == 1 & confounders_age == 1 & confounders_smoking == 0, cv_adjusted_1 := 1]
  data[cv_adjusted_0 == 0 & cv_adjusted_1 == 0, cv_adjusted_2 := 1]

  data$cv_incidence <- 0
  data$cv_mortality <- 0
  table(data$outcome_type)
  data[outcome_type == "Incidence", cv_incidence := 1]
  data[outcome_type == "Mortality", cv_mortality := 1]
  data[outcome_type == "Incidence & Mortality", cv_incidence := 1]
  data[outcome_type == "Incidence & Mortality", cv_mortality := 1]
  table(data$cv_incidence)
  table(data$cv_mortality)

  data[, cv_sick_quitters := 1]
  data[reverse_causation == 0, cv_sick_quitters := 0]
  data[, rep_prevalent_disease := 0]
}

# create indicator variable for number of other confounders adjusted for
if (T) {
  data$id <- 1:nrow(data)
  confounder_col_names <- grep("^extract|confounder|id|^nid|^outcome$",
    colnames(data),
    value = TRUE
  )
  confounders <- data[, ..confounder_col_names]

  confounders[, confounders_other := gsub("(?<!log)\\([^\\)]+(\\)|\\s)",
    "",
    confounders_other,
    perl = TRUE
  )][, confounders_other := gsub("\\.$", "", confounders_other, perl = TRUE)]

  other_confounders <- confounders[, c(tstrsplit(confounders_other,
    split = ",( and)?|;",
    perl  = TRUE
  ))]

  other_confounders <- other_confounders[
    grep("^Race and study center$|^Caffeine intake and gender$|^Use of NSAIDs and aspirin$|^Social class and parity$",
      V1,
      perl = TRUE
    ),
    c("V1", "V2") := tstrsplit(V1, split = " and ", perl = TRUE)
  ]
  other_confounders <- as.matrix(other_confounders) %>%
    trimws() %>%
    as.data.table()

  other_confounders$id <- 1:nrow(other_confounders)

  write.csv("filepath", row.names = F)
  other_covs <- read.csv("filepath")
  other_covs$covariates <- other_covs$value
  other_covs$value <- NULL
  other_covs <- unique(other_covs)

  other_confounders_long <- melt(other_confounders, id.vars = "id", value.name = "covariates")
  other_confounders_long$variable <- NULL
  other_confounders_long <- unique(other_confounders_long)
  other_confounders_long <- other_confounders_long[!is.na(covariates)]

  other_confounders_long <- merge(other_confounders_long, other_covs, by = "covariates", all.x = T)
  other_confounders_long$covariates <- NULL
  other_confounders_long$value <- 1

  other_confounders_wide <- dcast(other_confounders_long, id ~ value)
  names(other_confounders_wide)[names(other_confounders_wide) == "1"] <- "number_other_confounders"

  data <- merge(data, other_confounders_wide, by = "id", all = T)

  # create number of data points contributed per study
  data[, field_citation_value_substr := substr(field_citation_value, 1, 15)]
  data[, study_size := length(field_citation_value), by = .(field_citation_value_substr, location_name)]

  test <- data %>%
    group_by(nid) %>%
    summarise(count = n())
  hist(test$count)

  unique(data$underlying_nid)

  data_sub <- data[outcome == "Ischemic heart disease"]

  covs_summary <- data_sub[, .SD, .SDcols = c(names(data)[names(data) %like% "confounders" & names(data) != "confounders_other"], "nid")] %>% unique()
  covs_string <- names(data)[names(data) %like% "confounders" & names(data) != "confounders_other"]
  covs_summary[, (covs_string) := lapply(.SD, as.numeric), .SDcols = covs_string]
  covs_counts <- data.table(reshape2::melt(covs_summary[, colSums(.SD, na.rm = T) / nrow(.SD), .SDcols = names(data)[names(data) %like% "confounders" & names(data) != "confounders_other"]]))
  covs_counts[, confounder := names(data)[names(data) %like% "confounders" & names(data) != "confounders_other"]]
  setnames(covs_counts, c("value"), c("prevalence"))

  nid_counts <- melt(covs_summary, id.vars = "nid", variable.name = "confounder")
  nid_counts <- nid_counts[value == 1]
  nid_counts <- unique(nid_counts)
  nid_counts <- nid_counts %>%
    group_by(confounder) %>%
    summarise(nid_count = n())

  covs_counts <- merge(covs_counts, nid_counts, by = "confounder", all.x = T)

  covs_counts[, confounder := factor(confounder, levels = covs_counts[order(prevalence)]$confounder)]
}

# create other variables
if (T) {
  data$age_start <- as.numeric(data$age_start)
  data$age_end <- as.numeric(data$age_end)
  data$age_mean <- as.numeric(data$age_mean)

  data[, age_mid := (age_start + age_end) / 2]

  data[!is.na(age_mid), age := age_mid]
  data[!is.na(age_mean), age := age_mean]

  # categorize age
  data[age >= 15 & age < 20, age_cat := "15-19"]
  data[age >= 20 & age < 25, age_cat := "20-24"]
  data[age >= 25 & age < 30, age_cat := "25-29"]
  data[age >= 30 & age < 35, age_cat := "30-34"]
  data[age >= 35 & age < 40, age_cat := "35-39"]
  data[age >= 40 & age < 45, age_cat := "40-44"]
  data[age >= 45 & age < 50, age_cat := "45-49"]
  data[age >= 50 & age < 55, age_cat := "50-54"]
  data[age >= 55 & age < 60, age_cat := "55-59"]
  data[age >= 60 & age < 65, age_cat := "60-64"]
  data[age >= 65 & age < 70, age_cat := "65-69"]
  data[age >= 70 & age < 75, age_cat := "70-74"]
  data[age >= 75 & age < 80, age_cat := "75-79"]
  data[age >= 80 & age < 85, age_cat := "80-84"]
  data[is.na(age_cat), age_cat := "No age entered"]

  # sample size calculations
  data$cohort_sample_size_total <- as.numeric(data$cohort_sample_size_total)
  data$cohort_sample_size_exp <- as.numeric(data$cohort_sample_size_exp)
  data$cohort_sample_size_unexp <- as.numeric(data$cohort_sample_size_unexp)
  data$age_sd <- as.numeric(data$age_sd)

  data[is.na(cohort_sample_size_total), cohort_sample_size_total := cohort_sample_size_exp + cohort_sample_size_unexp]
  data[is.na(age_start), age_lower := age_mean - 1.96 * (age_sd / sqrt(cohort_sample_size_total))]
  data[is.na(age_start), age_upper := age_mean + 1.96 * (age_sd / sqrt(cohort_sample_size_total))]

  # exposure and outcome assessment indicator variables
  data[exp_method_1 == "Self-report (human/environment)", cv_exposure_selfreport := 1]
  data[outcome_assess_1 == "Self-report" & is.na(outcome_assess_2) & is.na(outcome_assess_3), cv_outcome_selfreport := 1]
  data[!(outcome_assess_1 == "Self-report" & is.na(outcome_assess_2) & is.na(outcome_assess_3)), cv_outcome_selfreport := 0]
  data[!is.na(exp_assess_num), cv_exposure_study := 1]
  data[is.na(exp_assess_num), cv_exposure_study := 0]

  data <- data[!is.na(exp_level_lower) | !is.na(exp_level_value)]

  data[is.na(unexp_level_lower) & unexposed_group == "Non-drinkers", unexp_level_lower := 0]
  data[is.na(unexp_level_upper) & unexposed_group == "Non-drinkers", unexp_level_upper := 0]

  data[percent_male %in% c(0, 1), confounders_sex := 1]

  # sublocation and representativeness indicator variables
  data <- setnames(data, old = c("smaller_site_unit", "representativeness"), new = c("cv_rep_geography", "cv_subpopulation"))

  # uncertainty intervals / standard errors
  data[mean > lower, standard_error := (mean - lower) / 1.96]
  data[, standard_error := as.numeric(standard_error)]
  data[, upper_se := quantile(standard_error, 0.975, na.rm = T)]
  data[is.na(lower) & is.na(upper), lower := mean - 1.96 * upper_se]
  data[is.na(lower) & is.na(upper), upper := mean + 1.96 * upper_se]
}

# create confounder character variables
names(data)[grep("confounder", names(data))]

changeCols <- c(
  "confounders_age", "confounders_sex", "confounders_education", "confounders_income",
  "confounders_smoking", "confounders_alcohol_use", "confounders_physical_activity",
  "confounders_dietary_components", "confounders_bmi", "confounders_hypertension", "confounders_diabetes",
  "confounders_hypercholesterolemia"
)
data[, (changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]

if (T) {
  data[, confounders_age := ifelse(confounders_age == 1, "age", 0)]
  data[, confounders_sex := ifelse(confounders_sex == 1, "sex", 0)]
  data[, confounders_education := ifelse(confounders_education == 1, "education", 0)]
  data[, confounders_income := ifelse(confounders_income == 1, "income", 0)]
  data[, confounders_smoking := ifelse(confounders_smoking == 1, "smoking", 0)]
  data[, confounders_alcohol_use := ifelse(confounders_alcohol_use == 1, "alc_use", 0)]
  data[, confounders_physical_activity := ifelse(confounders_physical_activity == 1, "phys_act", 0)]
  data[, confounders_dietary_components := ifelse(confounders_dietary_components == 1, "diet", 0)]
  data[, confounders_bmi := ifelse(confounders_bmi == 1, "bmi", 0)]
  data[, confounders_hypertension := ifelse(confounders_hypertension == 1, "hypertens", 0)]
  data[, confounders_diabetes := ifelse(confounders_diabetes == 1, "diabetes", 0)]
  data[, confounders_hypercholesterolemia := ifelse(confounders_hypercholesterolemia == 1, "hyperchol", 0)]
  data[!is.na(confounders_other), confounders_other_ind := "other"]

  data[, confounder_label := paste0(
    confounders_age, ", ", confounders_sex, ", ", confounders_education, ", ", confounders_income, ", ",
    confounders_smoking, ", ", confounders_alcohol_use, ", ", confounders_physical_activity, ", ",
    confounders_dietary_components, ", ", confounders_bmi, ", ", confounders_hypertension, ", ",
    confounders_diabetes, ", ", confounders_hypercholesterolemia, ", ", confounders_other_ind
  )]
}

# save data
data <- data[subgroup_analysis == "0"]
write.csv(data, paste0("filepath"), row.names = F)