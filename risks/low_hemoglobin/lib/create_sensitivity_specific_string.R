#' Create outcome- and sensitivity-specific string to be used in paths
create_sensitivity_specific_string <- function(outcome,
                                               outcome_severity = NULL,
                                               trimester = NULL,
                                               cov_con_socio = NULL,
                                               cov_rev = NULL) {
  
  if (!is.null(outcome_severity) && outcome_severity != "") {
    outcome <- paste0(outcome, "-", outcome_severity)
  }
  
  if (!is.null(trimester) && trimester != "") {
    outcome <- paste0(outcome, "-trim", trimester)
  }
  
  if (!is.null(cov_con_socio) &&
      length(cov_con_socio) == 1 &&
      cov_con_socio %in% 0:1) {
    outcome <- paste0(outcome, "-cov_con_socio_", cov_con_socio)
  }
  
  if (!is.null(cov_rev) && cov_rev != "") {
    cov_rev <- gsub("[^0-9]", "", cov_rev)
    if (cov_rev == 1) {
      outcome <- paste0(outcome, "-cov_rev_1")
    } else if (cov_rev == "01" || cov_rev == "10") {
      outcome <- paste0(outcome, "-cov_rev_0_1")
    }
  }
  
  return(outcome)
}
