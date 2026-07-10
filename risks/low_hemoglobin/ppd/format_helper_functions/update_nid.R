##' *************************************************************************************
##' Title: format_dismod_bundle.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# update NID values -------------------------------------------------------

update_nid <- function(u_nid, nid){
  i_vec <- which(is.na(nid) | nid == "")
  nid[i_vec] <- u_nid[i_vec]
  return(nid)
}