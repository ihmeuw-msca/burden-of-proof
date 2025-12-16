#' Create title fr MR-BRT plots using values of `trimester` and bias cov columns
create_mrbrt_plot_title <- function(outcome_title,
                                    hb_risk_type,
                                    cov_con_socio = NULL,
                                    cov_rev = NULL,
                                    outcome_severity = NULL,
                                    trimester = NULL,
                                    incl_undef = NULL) {

  if (outcome_severity == "vlbw") {
    outcome_title <- paste("very", outcome_title)
  } else if (outcome_severity == "elbw") {
    outcome_title <- paste("extremely", outcome_title)
  } else if (outcome_severity == "mptb") {
    outcome_title <- paste("moderate", outcome_title)
  } else if (outcome_severity == "vptb") {
    outcome_title <- paste("very", outcome_title)
  } else if (outcome_severity == "eptb") {
    outcome_title <- paste("extremely", outcome_title)
  } else if (outcome_severity == "postnatal_dep") {
    outcome_title <- "Postnatal Depression"
  } else if (outcome_severity == "prenatal_dep") {
    outcome_title <- "Prenatal Depression"
  }
  
  forest_plot_title <- glue::glue(
    "Meta-analysis of {stringr::str_to_lower(hb_risk_type)} hemoglobin as a risk factor for {outcome_title}"
  )
  
  if (!is.null(trimester) && !is.na(trimester) && trimester != "") {
    forest_plot_title <- glue::glue(
      "{forest_plot_title} (trimester {trimester})"
    )
  }
  
  if (!is.null(cov_con_socio) && length(cov_con_socio) == 1 && cov_con_socio == 0) {
    forest_plot_title <- glue::glue(
      "{forest_plot_title} (restricted to observations with adjustment for socioeconomic confounders)"
    )
  }
  
  if (is.null(cov_rev) || (!is.null(cov_rev) && cov_rev != 0 && cov_rev != 1)) {
    forest_plot_title <- glue::glue(
      "{forest_plot_title} (agnostic to reverse causality)"
    )
  }

  if (!is.null(incl_undef) && length(incl_undef) == 1 && incl_undef == 1) {
    forest_plot_title <- glue::glue(
      "{forest_plot_title} (including observations with undefined outcomes)"
    )
  }
  
  return(forest_plot_title)
}
