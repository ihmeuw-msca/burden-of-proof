#' Subset data to desired outcome or outcome subgroup, if relevant

indicate_outcome <- function(dat,
                             outcome,
                             outcome_severity = "") {
  box::use(
    data.table
  )
  
  checkmate::assert_subset(x = outcome,
                           choices = c("preecl", "lbw", "lga", "mat_hem",
                                       "mat_mort", "mat_sepsis", "neo_mort",
                                       "neo_sepsis", "ppd", "ptb", "sga", 
                                       "stillbirth", "hdop", "ghtn"),
                           empty.ok = FALSE)
  
  checkmate::assert_subset(x = outcome_severity,
                           choices = c("mptb", "vptb", "eptb", "vlbw", "elbw", 
                                       "early_neonatal", "late_neonatal", 
                                       "antenatal_hem", "prenatal_dep",
                                       "postnatal_dep", ""),
                           empty.ok = TRUE)
  
  if (is.null(outcome_severity) || is.na(outcome_severity)) {
    outcome_severity <- ""
  }
  
  if (outcome == "lbw") {
    dat <- indicate_lbw(dat, outcome, outcome_severity)
  } else if (outcome %in% c("sga",
                            "lga")) {
    dat <- indicate_sga_lga(dat, outcome, outcome_severity)
  } else if (outcome %in% c("mat_hem",
                            "mat_mort",
                            "mat_sepsis",
                            "ppd",
                            "preecl",
                            "ghtn",
                            "hdop")) {
    dat <- indicate_maternal_outcome(dat, outcome, outcome_severity)
  } else if (outcome == "neo_mort") {
    dat <- indicate_neo_mortality(dat, outcome, outcome_severity)
  } else if (outcome == "neo_sepsis") {
    dat <- indicate_neo_sepsis(dat, outcome, outcome_severity)
  } else if (outcome == "ptb") {
    dat <- indicate_ptb(dat, outcome, outcome_severity)
  } else if (outcome == "stillbirth") {
    dat <- indicate_stillbirth(dat, outcome, outcome_severity)
  }
  
  dat <- data.table::setDT(dat)
  return(dat)
}


#' Indicate LBW severity and subset as needed
indicate_lbw <- function(dat, outcome, outcome_severity = "") {
  # Indicate and exclude high birthweight (macrosomia)
  dat[, macrosomia_ind := ifelse((out_unit == "grams" &
                                    !is.na(out_lower) &
                                    out_lower >= 4000), 
                                 1, 
                                 0)]
  dat <- dat[macrosomia_ind == 0]
  
  if (outcome_severity == "vlbw") {
    dat[, vlbw_ind := ifelse((out_unit == "grams" &
                                out_upper >= 1499 &
                                out_upper <= 1500), 1, 0)]
    dat <- dat[vlbw_ind == 1]
  } else if (outcome_severity == "elbw") {
    dat[, elbw_ind := ifelse((out_unit == "grams" & 
                                out_upper >= 999 &
                                out_upper <= 1000), 1, 0)]
    dat <- dat[elbw_ind == 1]
  } else {
    dat[, lbw_ind := ifelse((out_unit == "grams" & 
                               out_upper >= 2499 & 
                               out_upper <= 2500), 1, 0)]
    dat <- dat[lbw_ind == 1]
  }
  return(dat)
}


#' Indicate SGA/LGA
indicate_sga_lga <- function(dat, outcome, outcome_severity = "") {
  if (outcome == "sga") {
    dat <- dat[, sga_ind := ifelse((out_unit == "percentile" &
                                      out_upper >= 9.9 &
                                      out_upper <= 10), 1, 0)] 
    dat <- dat[sga_ind == 1]
  } else if (outcome == "lga") {
    dat <- dat[, lga_ind := ifelse((out_unit == "percentile" & 
                                      out_lower >= 90 &
                                      out_lower <= 90.1), 1, 0)] 
    dat <- dat[lga_ind == 1]
  }
  return(dat)
}


#' Indicate maternal outcomes and subset as needed
indicate_maternal_outcome <- function(dat, outcome, outcome_severity = "") {
  # Maternal mortality
  if (outcome == "mat_mort") {
    dat <- dat[out_type == "Maternal mortality"]
    
    # Maternal hemorrhage
  } else if (outcome == "mat_hem") {
    if (outcome_severity == "antenatal_hem") { 
      dat <- dat[out_type == "Antenatal Hemorrhage"]
    } else {
      dat <- dat[out_type %in% c("Postpartum Hemorrhage - Primary",
                                 "Postpartum Hemorrhage - Unspecified")]
    }
    
    # Maternal sepsis
  } else if (outcome == "mat_sepsis") {
    dat <- dat[out_type =="Maternal sepsis"]
    
    # Preeclampsia
  } else if (outcome == "preecl") {
    dat <- dat[outsub_htn == "Pre-eclampsia"]
    
  } else if (outcome == "hdop") {
    dat <- dat[out_type == "Hypertensive disorders of pregnancy"]
    
  } else if (outcome == "ghtn") {
    dat <- dat[outsub_htn == "Gestational hypertension or pregnancy-induced hypertension"]
    
    # PPD
  } else if (outcome == "ppd") {
    dat <- dat[analysis_exclcrx == 1]
    if (outcome_severity == "prenatal_dep") {
      dat <- dat[prenataldep == 1]
    } else if (outcome_severity == "postnatal_dep") {
      dat <- dat[postnataldep == 1]
    }
  }
  return(dat)
}


#' Indicate neonatal mortality subtype and subset as needed
indicate_neo_mortality <- function(dat, outcome, outcome_severity = "") {
  if (outcome == "neo_mort") {
    dat$out_type_original <- dat$out_type
    dat[out_type_original %in% c("Neonatal mortality", "Perinatal mortality") &
          out_timepoint %in% c("0-6days of life",
                               "0-6 days of life",
                               "within 24 hrs of birth"),
        out_type := "early neonatal"]
    dat[(out_type_original == "Neonatal mortality" &
           out_timepoint == "7-27 days of life"),
        out_type := "late neonatal"]
    dat[out_type_original %in% c("Neonatal mortality", "Perinatal mortality") &
          (out_timepoint %in% c("0-27 days of life",
                                "0-28 days of life",
                                "0-28 days",
                                "unspecified") |
             is.na(out_timepoint)),
        out_type := "total neonatal"]
    
    if (outcome_severity == "early_neonatal") {
      dat <- dat[out_type == "early neonatal"]
    } else if (outcome_severity == "late_neonatal") {
      dat <- dat[out_type == "late neonatal"]
    } else {
      dat <- dat[out_type %in% c("total neonatal", "early neonatal", "late neonatal")]
    }
  }
  return(dat)
}


#' Indicate neonatal sepsis and subset as needed
indicate_neo_sepsis <- function(dat, outcome, outcome_severity = "") {
  if (outcome == "neo_sepsis") {
    dat$out_type_original <- dat$out_type
    message("No neonatal sepsis indicator provided, all observations considered neonatal-sepsis-related")
  }
  return(dat)
}


#' Indicate PTB severity and subset as needed
indicate_ptb <- function(dat, outcome, outcome_severity = "") {
  if (outcome_severity == "vptb") {
    dat[, vptb_ind := ifelse((out_unit == "weeks" & 
                                (out_lower >= 27.9 & out_lower <= 28) &
                                (out_upper >= 31.9 & out_upper <= 32)), 1, 0)]
    dat <- dat[vptb_ind == 1]
  } else if (outcome_severity == "mptb") {
    dat[, mptb_ind := ifelse((out_unit == "weeks" &
                                (out_lower >= 31.9 & out_lower <= 32) &
                                (out_upper >= 36.9 & out_upper <= 37)), 1, 0)]
    dat <- dat[mptb_ind == 1]
  } else if (outcome_severity == "eptb") {
    dat[, eptb_ind := ifelse((out_unit == "weeks" & 
                                out_upper >= 27.9 &
                                out_upper <= 28), 1, 0)]
    dat <- dat[eptb_ind == 1]
  } else {
    dat[, ptb_ind := ifelse((out_unit == "weeks" & 
                               out_upper >= 36.9 & 
                               out_upper <= 37), 1, 0)]
    dat <- dat[ptb_ind == 1]
  }
  return(dat)
}


#' Indicate stillbirth and subset as needed
indicate_stillbirth <- function(dat, outcome, outcome_severity = "") {
  if (outcome == "stillbirth") {
    dat$out_type_original <- dat$out_type
    dat <- dat[(out_type_original == "Stillbirth" |
                 (out_type_original %in% c("Neonatal mortality",
                                           "Perinatal mortality") &
                    out_timepoint == "at delivery")) &
                 out_lower == 20,
               out_type := "stillbirth"]
    dat <- dat[out_type %in% c("stillbirth")]
  }
  return(dat)
}
