#' Perform outcome-specific one-off cleaning
#' 
#' @note Once a given change is made to the raw extractions file, remove it from
#'   the source code below.

outcome_specific_cleaning <- function(dat, outcome) {
  # LBW
  if (outcome == "lbw") {  
    # Clean bio_unit_1 column
    dat[, bio_unit_1 := ifelse(
      bio_unit_1 %in% c("mg/dL", "mg/dl"), "g/dL", bio_unit_1
    )]
    # Fix one-off issue - fill NA with exposure bounds unit
    dat[record_nid == 551792 & is.na(bio_unit_1), bio_unit_1 := "g/dl"]
    
    # Convert kilograms to grams for outcome indication in later step
    dat$out_upper <- as.numeric(dat$out_upper)
    dat <- dat[!(is.na(out_unit)) & out_unit == "kilograms", `:=` (
      out_upper = out_upper * 1000,
      out_lower = out_lower * 1000,
      out_unit = "grams"
    )]
    
  # PTB
  } else if (outcome == "ptb") {
    # Fix to replace implausible bio_unit_1
    dat[record_nid == 530809 & tolower(bio_unit_1) == "mg/dl",
        bio_unit_1 := "g/dl"]
    dat[record_nid == 562850 & tolower(bio_unit_1) == "mg/dl",
        bio_unit_1 := "g/dl"]
    
    # Remove obs with ambiguous exp_timepoint
    dat <- dat[-(189:196), ]
    
    dat$out_upper <- as.numeric(dat$out_upper)
    dat <- dat[!(is.na(out_unit)) & out_unit == "days", `:=`(
      out_upper = out_upper / 7,
      out_lower = out_lower / 7,
      out_unit = "weeks"
    )]
    
  # SGA
  } else if (outcome == "sga") {
    # Clean bio_unit_1 column
    dat[177, bio_unit_1 := "g/dL"]
    dat[179, bio_unit_1 := "g/dL"]
    dat[180, bio_unit_1 := "g/dL"]
    dat[181, bio_unit_1 := "g/dL"]
    
  # Stillbirth or neonatal mortality
  } else if (outcome %in% c("stillbirth", "neo_mort")) {
    dat <- dat[, bio_unit_1 := ifelse(
      (bio_unit_1 %in% c("mg/dL", "mg/dl")),
      "g/dL",
      bio_unit_1
    )]
    
    bio_type_vec <- c("Hemoglobin", "hemoglobin", "Hematocrit", "Packed_Cell_Volume")
    dat <- dat[bio_type_1 %in% bio_type_vec & is.na(bio_type_2)]
    dat <- dat[!(bio_form_1 %in% c("Continuous", "Categorical, quintile"))]
    
  # Maternal mortality or maternal sepsis
  } else if (outcome %in% c("mat_mort", "mat_sepsis", "preecl", "ghtn", "hdop")) {
    location_sheet <- data.table::fread(
      "PATH TO LOCATION DATASET"
    )
    names(dat)[names(dat) == "location_name"] <- "location_id" # so the hb imputation would work
    location_sheet$location_name_placeholder <- paste(location_sheet$location_name, location_sheet$ihme_loc_id, sep = "|")
    match_indices <- match(dat$location_id, location_sheet$location_id)
    dat$location_name <- location_sheet$location_name_placeholder[match_indices]
    
  # Maternal sepsis only
  } else if (outcome == "mat_sepsis") {
    names(dat)[names(dat) == "location_name"] <- "location_id"
    dat$bio_unit_1_hct[dat$redcap_id == 45] <- "%" # missing unit for hct
    dat$bio_unit_1_hgb[dat$redcap_id == 49] <- "g/dl" # though reported in paper as mg/dl, numeric value suggests g/dL
    
  # Maternal hemorrhage
  } else if (outcome == "mat_hem") {
    dat[, bio_unit_1 := ifelse((bio_unit_1 %in% c("mg/dL", "mg/dl")), "g/dL", bio_unit_1)]
  
  # Hypertensive disorders
  } else if (outcome %in% c("hdop", "ghtn", "preecl")) {
    dat[c(25, 26, 27, 68, 102), bio_unit_1_hgb := "g/dl"]
    dat[c(103, 104), bio_unit_1_hgb := "mmol/l"]
    dat <- dat[!(bio_unit_1 %in% c("Serum ferritin", "MCV"))]
    
  # PPD  
  } else if (outcome == "ppd") {
    dat$out_upper <- NA
    dat$out_lower <- NA
    dat$out_unit <- NA
  }
  
  return(dat)
}
