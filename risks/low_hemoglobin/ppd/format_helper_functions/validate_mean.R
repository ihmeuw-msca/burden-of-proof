##' *************************************************************************************
##' Title: validate_mean.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# assign mean column ------------------------------------------------------

mean_hot_fix <- function(input_df){
  df <- copy(input_df)
  i_vec <- (is.na(df$cases) & is.na(df$mean)) |
    (
      df$measure == 'prevalence' &
        (df$mean < 0 | df$mean > 1)
    )
  i_vec <- which(!i_vec)
  
  return(df[i_vec, ])
}


assign_mean <- function(mean_vec, measure_vec, cases_vec){
  i_vec <- (is.na(cases_vec) & is.na(mean_vec))  |
    (
      (measure_vec == 'prevalence' & !is.na(mean_vec)) &
        (mean_vec < 0 | mean_vec > 1)
    )
  validation_flag <- any(i_vec)
  
  if(!validation_flag){
    return(mean_vec)
  }else{
    assign("invalid_mean_index", which(i_vec), envir = .GlobalEnv)
    stop("Mean is undefined or or has invalid bounds. Please check and rerun.")
  }
}