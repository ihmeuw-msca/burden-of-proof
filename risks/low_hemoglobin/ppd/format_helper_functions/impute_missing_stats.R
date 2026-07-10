##' *************************************************************************************
##' Title: impute_missing_stats.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# impute missing stats ----------------------------------------------------

impute_missing_statistics <- function(input_df, se_hot_fix){
  df <- copy(input_df)
  
  df <- impute_cases(df)
  df <- impute_standard_error(df, se_hot_fix)
  df <- impute_upper_lower(df)
  # df$variance <- calculate_variance(df$standard_error)
  
  return(df)
}


impute_cases <- function(input_df){
  df <- copy(input_df)
  
  i_vec <- which(
    df$measure == "prevalence" &
      is.na(df$cases)
  )
  
  df$cases[i_vec] <- df$mean[i_vec] * df$sample_size[i_vec] ###TSN note - object of type 'closure' is not subsettable 
  
  return(df)
}

impute_standard_error <- function(input_df, se_hot_fix){
  df <- copy(input_df)
  i_vec <- which(
    df$measure == "prevalence" &
      !(is.na(df$mean)) &
      !(is.na(df$sample_size)) &
      df$sample_size > 0 &
      is.na(df$standard_error)
  )
  
  df$standard_error[i_vec] <- sqrt(
    df$mean[i_vec] * (1 - df$mean[i_vec]) / df$sample_size[i_vec]
  )
  
  i_vec <- which(
    !(is.na(df$upper)) &
      !(is.na(df$lower)) &
      is.na(df$standard_error)
  )
  
  df$standard_error[i_vec] <- (df$upper[i_vec] - df$lower[i_vec]) / (2 * qnorm(0.975))
  
  if(se_hot_fix){
    i_vec <- which(
      is.na(df$standard_error) &
        df$sample_size < 1000
    )
    df$standard_error[i_vec] <- .0025
    
    i_vec <- which(
      is.na(df$standard_error) &
        df$sample_size >= 1000
    )
    df$standard_error[i_vec] <- .00125
  }
  
  return(df)
}

impute_upper_lower <- function(input_df){
  df <- copy(input_df)
  
  bound_vec <- c("lower", "upper")
  for(b in 1:length(bound_vec)){
    i_vec <- which(
      !(is.na(df$mean)) &
        !(is.na(df$standard_error)) &
        is.na(df[[bound_vec[b]]])
    )
    
    df[[bound_vec[b]]][i_vec] <- df$mean[i_vec] + df$standard_error[i_vec] * qnorm(0.975) * ((-1) ^ b)
    df$uncertainty_type_value[i_vec] <- 95
  }
  
  #quick fix - for lowers under 0, code them to 0
  i_vec <- which(
      df$lower < 0
  )
  df$lower[i_vec] <- 0

  
  return(df)
}
