#'  Subset by specified value(s) of given covariate column(s) to facilitate
#'  sensitivity analyses

subset_by_given_cov_value <- function(dat,
                                      cov_con_socio = NULL,
                                      cov_rev = NULL) {
  
  if (!is.null(cov_rev) && cov_rev == 0) {
    dat <- dat[dat$cov_rev == 0, ]
    message("Subsetting by cov_rev = 0")
  } else if (!is.null(cov_rev) && cov_rev == 1) {
    dat <- dat[dat$cov_rev == 1, ]
    message("Subsetting by cov_rev = 1")
  } else {
    message("No subsetting by cov_rev indicated")
  }
  
  if (!is.null(cov_con_socio) && cov_con_socio == 0) {
    dat <- dat[dat$cov_con_socio == 0, ]
    message("Subsetting by cov_con_socio = 0")
  } else if (!is.null(cov_con_socio) && cov_con_socio == 1) {
    dat <- dat[dat$cov_con_socio == 1, ]
    message("Subsetting by cov_con_socio = 1")
  } else {
    message("No subsetting by cov_con_socio indicated")
  }

  return(dat)
}
