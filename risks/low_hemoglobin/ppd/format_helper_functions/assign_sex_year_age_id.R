##' *************************************************************************************
##' Title: add_sex_year_age_id.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)


# add sex id/sex to df ----------------------------------------------------

add_sex_id <- function(input_df){
  df <- copy(input_df)
  
  df$sex_id<-NA
  
  sex_list <- list(
    Male = 1,
    Female = 2,
    Both = 3
  )
  
  for(i in names(sex_list)){
    i_vec <- which(
      is.na(df$sex) & 
        !(is.na(df$sex_id)) & 
        df$sex_id == sex_list[[i]]
    )
    df$sex[i_vec] <- i
    
    i_vec <- which(
      is.na(df$sex_id) & 
        !(is.na(df$sex)) & 
        df$sex == i
    )
    df$sex_id[i_vec] <- as.integer(sex_list[[i]])
  }
  
  if(any(is.na(df$sex_id) | is.na(df$sex))){
    stop("Sex IDs/Sex columns are not all defined. Please update and rerun.")
  }
  
  return(df)
}



# assign year id ----------------------------------------------------------

validate_year_id <- function(df){
  if(!all(df$year_start <= df$year_end)){
    stop("Not all year_start are less than or equal to year_end.")
  }
}


# assign age group ids ----------------------------------------------------

validate_age_group_ids <- function(df){
  
  validation_flag <- any(
    is.na(df$age_start) |
      is.na(df$age_end) |
      df$age_start > df$age_end
  )
  
  if(validation_flag){
    stop("age_start, age_end not defined or age_start>age_end. Please check and rerun.")
  }
}

