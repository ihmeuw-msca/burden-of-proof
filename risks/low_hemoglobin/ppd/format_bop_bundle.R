# source libraries --------------------------------------------------------
library(data.table)
library(stringr)

# source helper functions ---------------------------------------------------
invisible(sapply(
  list.files(
    path = "PATH TO format_helper_functions FOLDER PATH",
    full.names = TRUE,
    pattern = "*\\.R$"
  ),
  source
))

# format bundles for dismod -----------------------------------------------
format_bop_bundle <- function(input_df, release_id, hot_fix = FALSE){
  df <- copy(input_df)
  validate_year_id(df)
  validate_age_group_ids(df)
  
  df <- add_sex_id(input_df = df)
  df$nid <- update_nid(u_nid = df$underlying_nid, nid = df$nid)
  df$is_outlier <- 0 
  
  return(df)
}
