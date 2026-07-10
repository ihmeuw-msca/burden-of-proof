##' *************************************************************************************
##' Title: check_measure.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# update measure value ----------------------------------------------------

update_measure <- function(measure_col){
  if(any(is.na(vec))){
    stop("Not all measures defined. Please check and rerun.")
  }
  
  return(vec)
}