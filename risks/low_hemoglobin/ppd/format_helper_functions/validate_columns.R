
##' *************************************************************************************
##' Title: validate_columns.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks - makes sure that all 
##'         columns required for dismod validations are present 
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# source helper functions ---------------------------------------------------
invisible(sapply(
  list.files(
    path = "repos/hemog/nonfatal_pipeline/format_bundle_data/",
    full.names = TRUE,
    pattern = "*.R"
  ),
  source
))

# validate columns --------------------------------------------------------

validate_columns <- function(cols){
  req_cols <- c(
    "nid", "underlying_nid", "location_id", "sex", "measure",
    "year_start", "year_end", "age_start", "age_end", 
    "mean", "standard_error", "upper", "lower",
    "source_type", "representative_name", "urbanicity_type", "recall_type", 
    "unit_type", "unit_value_as_published", "is_outlier", "seq", "input_type", 
    "design_effect", "recall_type_value", "uncertainty_type", "sampling_type", 
    "effective_sample_size"
  )
  
  if(!(all(req_cols %in% cols))){
    i_vec <- which(!(req_cols %in% cols))
    col_string <- paste(req_cols[i_vec], collapse = ", ")
    stop("The follow required columns are not in your df: ", col_string)
  }
}