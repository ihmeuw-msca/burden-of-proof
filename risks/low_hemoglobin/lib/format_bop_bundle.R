#' Format bundles for dismod
format_bop_bundle <- function(input_df, release_id, hot_fix = FALSE) {
  box::use(
    transform_data/format_helper_functions/assign_sex_year_age_id,
    transform_data/format_helper_functions/update_nid
  )
  
  df <- data.table::copy(input_df)
  assign_sex_year_age_id$validate_year_id(df)
  assign_sex_year_age_id$validate_age_group_ids(df)
  
  df <- assign_sex_year_age_id$add_sex_id(input_df = df)
  df$nid <- update_nid$update_nid(u_nid = df$underlying_nid, nid = df$nid)
  df$is_outlier <- 0 
  
  return(df)
}
