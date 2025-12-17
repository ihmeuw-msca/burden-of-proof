#' Subset by trimester, if specified

subset_by_trimester <- function(dat,
                                trimester) {
  
  if (is.na(trimester)) {
    message("No trimester specified, using all relevant observations")
  } else if (trimester == 1) {
    dat <- dat[exp_timepoint == "first_trimester"]
  } else if (trimester == 2) {
    dat <- dat[exp_timepoint == "second_trimester"]
  } else if (trimester == 3) {
    dat <- dat[exp_timepoint == "third_trimester"]
  } else if (!is.na(trimester)) {
    stop("Invalid value provided for `trimester`")
  }
  
  if (nrow(dat) == 0) {
    stop("Stopping due to zero observations for trimester ", trimester)
  }
  
  return(dat)
}
