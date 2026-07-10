##' *************************************************************************************
##' Title: fill_source_type.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# update source type ------------------------------------------------------

update_source_type <- function(source_col){
  i_vec <- which(is.na(source_col) | source_col == "")
  source_col[i_vec] <- "Survey - other/unknown"
  return(source_col)
}