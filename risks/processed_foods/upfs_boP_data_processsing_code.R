## Upfs Cleaning Script
rm(list = ls())

## Libraries
library(data.table)
library(tidyverse)

# --- 1. Data Loading ---
dt2 <- fread("Filepath")
dt_nid <- fread("Filepath2")

# Clean NIDs
dt_nid <- dt_nid[keep == 1, .(nid, field_citation_value)]
dt_nid <- dt_nid[nid != "" & !is.na(nid)]
dt_nid2 <- unique(dt_nid, by = "field_citation_value")

# Merge and filter
dt <- merge(dt2, dt_nid2, by = "field_citation_value", all.x = TRUE)
dt <- dt[!is.na(nid)]

# Assign Global Variables: Assinged based on specific risk and cause
dt[, `:=`(
  risk_type = risk_type,
  risk_unit = risk_unit,
  acause    = acause
)]

# --- 2. Covariate Cleaning & Consistency ---
dt[, representative := ifelse(rep_geography == 1, 0, 1)]

# Clean inconsistencies in exposure level
dt[exp_level == "At the individual ", exp_level := "At the individual"]
dt[, `:=`(
  cov_exposure_asses_level     = ifelse(exp_level == "At the individual", 0, 1),
  cov_exposure_method_1        = ifelse(exp_method_1 == "Self-report (human/environment)", 0, 1),
  cov_outcome_asses_method_1   = ifelse(outcome_assess_1 == "Self-report", 1, 0),
  outcome_2                    = 0,
  reverse_causation            = 1,
  washout_years                = NA_real_,
  seq                          = 1:.N
)]

# Study design and bias covariates
dt[, `:=`(
  cov_selection_bias = ifelse(notes_selection_bias_follow_up >= 0.85, 0, 1),
  cov_study_design   = ifelse(design %in% c("Prospective cohort", "prospective cohort", "case-cohort", "Nested case-control"), 1, 0),
  cov_incidence_only = ifelse(outcome_type == "Incidence", 0, 1),
  cov_mortality_only = ifelse(outcome_type == "Mortality", 0, 1),
  cov_odds_ratio     = ifelse(effect_size_measure == "Odds ratio (OR)", 1, 0)
)]

# --- 3. Numeric Conversions & Transformations ---
num_cols <- c("effect_size", "upper", "lower", "value_of_duration_fup")
dt[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# Create Log-Relative Risk and SE
dt[, `:=`(
  ln_rr    = log(effect_size),
  ln_rr_se = (log(upper) - log(lower)) / 3.92
)]

# Follow-up and Risk levels
dt[, cov_follow_up := ifelse(value_of_duration_fup > 10, 1, 0)]
dt[, `:=`(
  alt_risk_lower = cohort_exp_level_rr_lower,
  alt_risk_upper = cohort_exp_level_rr_upper,
  ref_risk_lower = cohort_unexp_level_rr_lower,
  ref_risk_upper = cohort_unexp_level_rr_upper
)]

# Final numeric check for BoP columns
columns_to_convert <- c("ln_rr", "ln_rr_se", "ref_risk_lower", "ref_risk_upper", "alt_risk_lower", "alt_risk_upper")
dt[, (columns_to_convert) := lapply(.SD, as.numeric), .SDcols = columns_to_convert]
dt <- dt[!is.na(ln_rr)]

# --- 4. Column Renaming & Dropping ---
# Harmonize covariate names
names(dt) <- gsub("^(cofounder_+|confounder_+|confounders_+|cov_+)", "cov_", names(dt))
names(dt) <- gsub("^(cofounder|confounder|confounders|cov)_*", "cov_", names(dt))


# --- 5. Final Formatting for Upload ---
dt[, measure := "relrisk"]
dt[, sex := fcase(percent_male == 1, "Male", 
                  percent_male == 0, "Female", 
                  default = "Both")]

dt[, `:=`(
  crosswalk_parent_seq = NA,
  is_outlier           = 0,
  input_type           = NA,
  effect_size_unit     = "linear",
  source_type          = "Survey - longitudinal",
  standard_error       = "NA",
  rei                  = "UPF"
)]

# Harmonize string names
dt[design == "Prospective cohort", design := "prospective cohort"]
dt[effect_size_measure == "Hazard ratio (HR)", effect_size_measure := "hazard ratio"]
dt[effect_size_measure == "Relative risk (RR)", effect_size_measure := "relative risk"]
dt[effect_size_measure == "Odds ratio (OR)", effect_size_measure := "odds ratio"]

# Merge Location Metadata
source("FILEPATH/get_location_metadata.R")
loc <- get_location_metadata(location_set_id = 35, release_id = 16)
loc[, location_id := as.numeric(location_id)]
dt[, location_id := as.numeric(location_id)]
dt_final <- merge(dt, loc, by = "location_id", all.x = TRUE)

# --- 6. Automated Bias Covariate Filtering ---
# Remove constant covariates
cov_columns <- grep("^cov_", names(dt_final), value = TRUE)
constant_covs <- sapply(dt_final[, ..cov_columns], function(x) length(unique(na.omit(x))) == 1)
dt_final2 <- dt_final[, (cov_columns[constant_covs]) := NULL]

# Remove redundant/identical covariates
find_identical_columns_unique <- function(data) {
  col_names <- colnames(data)
  identical_columns <- c()
  if(length(col_names) < 2) return(NULL)
  for (i in 1:(length(col_names) - 1)) {
    for (j in (i + 1):length(col_names)) {
      if (all(data[[col_names[i]]] == data[[col_names[j]]], na.rm = TRUE)) {
        identical_columns <- c(identical_columns, col_names[j]) # Keep the first, drop the second
      }
    }
  }
  return(unique(identical_columns))
}

bias_covs <- grep("^cov_", names(dt_final2), value = TRUE)
overlap_to_drop <- find_identical_columns_unique(dt_final2[, ..bias_covs])
if(!is.null(overlap_to_drop)) dt_final2[, (overlap_to_drop) := NULL]

# Final specific drops

fwrite(dt_final2, paste0(save_d, folder, subfolder, "UPF", "-", acause, ".csv"))