#' Drop given column or columns

drop_cols <- function(dat, covs_to_drop) {
  covs_to_drop <- unlist(strsplit(covs_to_drop, "\\s*,\\s*|\\s+"))
  
  present_covs <- covs_to_drop[covs_to_drop %in% names(dat)]
  if (length(present_covs) > 0) {
    message("Dropping columns in covs_to_drop: ",
            paste(present_covs, collapse = ", "))
    dat <- dat |>
      dplyr::select(-dplyr::all_of(present_covs))
  }
  
  absent_covs <- setdiff(covs_to_drop, present_covs)
  if (length(absent_covs) > 0) {
    message("Columns not found and not dropped: ", paste(absent_covs, collapse = ", "))
  }
  
  return(dat)
}
