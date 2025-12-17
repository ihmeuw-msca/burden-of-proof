#' Create columns for acause, bundle_id, bundle_version_id, and out_sub
#' 
#' @param dat Data object containing extractions data for a given outcome.
#' @param outcome String containing outcome of interest. Must be one of: lbw,
#'   lga, mat_hem, mat_mort, mat_sepsis, neo_mort, neo_sepsis, ppd, ptb, sga,
#'   stillbirth.
#' 
#' @note The out_sub column helps differentiate between sub-levels of acause.

add_outcome_specific_cols <- function(dat, outcome) {
  dat$out_sub <- outcome
  dat$bundle_id <- 11111
  dat$bundle_version_id <- 11111
  
  if (outcome %in% c("lbw", "ptb", "sga", "lga")) {
    dat$cause_id <- 1061
  } else if (outcome == "stillbirth") {
    dat$cause_id <- 744
  } else if (outcome == "ppd") {
    dat$cause_id <- 567
  } else if (outcome == "neo_sepsis") {
    dat$cause_id <- 383
  } else if (outcome == "mat_sepsis") {
    dat$cause_id <- 368
  } else if (outcome == "mat_hem") {
    dat$cause_id <- 367
    dat$bundle_id <- 10313
    dat$bundle_version_id <- 42832
  } else if (outcome == "mat_mort") {
    dat$cause_id <- 366
  } else if (outcome == "neo_mort") {
    dat$cause_id <- 294
  } else if (outcome %in% c("hdop", "ghtn", "preecl")) {
    dat$cause_id <- 369
    dat$bundle_id <- 999
    dat$bundle_version_id <- 999
  }
  
  cause_dt <- data.table::as.data.table(ihme::get_ids("cause"))
  dat$acause <- cause_dt[cause_id == dat$cause_id[1], acause]
  
  return(dat)
}
