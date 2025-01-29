#######################################################################    
# Title: Create table for defined MR-BRT parameters for a BoP study
# Author: Sinclair Carr
#######################################################################    
rm(list = ls())

## create customer r library folder
user <- Sys.getenv("USER")
user_rlibs <- file.path("/homes", user, "rlibs")

if (!dir.exists(user_rlibs)) {
  dir.create(user_rlibs)
} 

## install/load packages
packages <- c("yaml", "data.table", "stringr", "officer", "flextable")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, lib = user_rlibs)
    library(p, lib.loc = user_rlibs, character.only = TRUE)
  }
}

## load cc functions
source("/FILEPATH/r/get_ids.R") 

## main function
#' This function creates a table for parameter specifications for MR-BRT to be shown in a Burden of Proof study
#'
#' @param input is the path to the folder where the results of the BoP pipeline are stored
#' @param output is the path to the folder where the table should be saved
#' @param heading optional heading of the table; "Table SX. MR-BRT spline and prior specifications" by default.
#'
#' @return
#' @export
#'
#' @examples
#' parameter_table(input = "/homes/shcarr/bop_pipeline/results/dichotomous/",
#'                 output = "/homes/shcarr/bop_pipeline/tables/", 
#'                 heading = "Table S7. MR-BRT spline and prior specifications")
#'                 
#'                 


risk_folder = "processed_meat"
trimming = "no_trimming/"
exposure = paste0(risk_folder, "consumption")

input = paste0("/FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/", trimming, "/risk_outcome_pair/")
output = paste0("/FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/", trimming, "/summary_table/")
heading = paste0("Table:", exposure ,"  MR-BRT spline and prior specifications")


parameter_table <- function(input, output, heading = "Table SX. MR-BRT spline and prior specifications"){
  
  if (!dir.exists(input)) stop("Please provide a valid input path")
  if (!dir.exists(output)) stop("Please provide a valid output path")
  
  # get cause_ids
  cause_ids <- get_ids("cause")
  
  results <- data.table()
  
  # loop through all risk-outcome pairs
  for (ro_pair in list.files(input)) {
    ro_path <- file.path(input, ro_pair)
    
    # load files
    settings <- yaml.load_file(file.path(ro_path, "settings.yaml"))
    summary <- yaml.load_file(file.path(ro_path, "summary.yaml"))
    
    if(summary$risk_type == "continuous"){
      
      # Spline degree and interior knots
      spline_degree <- switch(settings$fit_signal_model$cov_model$spline_degree, 
                              "linear", "quadratic", "cubic", "N/A")
      num_knots <- paste(
        length(settings$fit_signal_model$cov_model$spline_knots) - 2, 
        "interior knots")
      
      degree_knots <- paste(spline_degree, num_knots, sep = ", ")
      
      # Gaussian prior
      gauss_prior_mean <- settings$fit_signal_model$cov_model$prior_spline_maxder_gaussian[[1]]
      gauss_prior_sd <- settings$fit_signal_model$cov_model$prior_spline_maxder_gaussian[[2]]
      gauss_prior <- paste0("Gaussian max derivative prior on the right tail (", 
                            gauss_prior_mean[length(gauss_prior_mean)], ", ", 
                            gauss_prior_sd[length(gauss_prior_sd)], ")")
      
      # monotonicity constraint
      prior_spline_monotonicity <- settings$fit_signal_model$cov_model$prior_spline_monotonicity
      
      monotonicity <- ifelse(is.null(prior_spline_monotonicity), NA,
                             ifelse(prior_spline_monotonicity == "decreasing", 
                                    "monotonic decreasing", "monotonic increasing"))
      
      # linear tails
      r_linear_tail <- settings$fit_signal_model$cov_model$spline_r_linear
      l_linear_tail <- settings$fit_signal_model$cov_model$spline_l_linear
      
      if (!(r_linear_tail) & !(l_linear_tail)) {
        linear_tail <- NA
      } else if (r_linear_tail & !(l_linear_tail)) {
        linear_tail <- "right linear tail"
      } else if (!(r_linear_tail) & l_linear_tail) {
        linear_tail <- "left linear tail"
      } else {
        linear_tail <- "left and right linear tail"
      }
      
      priors_constraints <- paste0(ifelse(is.na(monotonicity), "", paste0(monotonicity, ", ")), 
                                   linear_tail, ", ", gauss_prior)
      
    } else {
      
      degree_knots <- "N/A"
      priors_constraints <- "N/A"
      
    }
    
    # exposure
    risk <- stringr::str_split(ro_pair, "-")[[1]][1]
    
    # health outcome
    cause <- stringr::str_split(ro_pair, "-")[[1]][2]
    health_outcome <- ifelse(cause %in% cause_ids$acause, cause_ids$cause_name[cause_ids$acause == cause], cause)
    
    temp <- data.frame(
      cbind(
        risk, 
        health_outcome,
        degree_knots,
        priors_constraints
      )
    )
    results <- rbind(results, temp)
  }
  
  # make first letters uppercase
  #results <- data.table(
   # sapply(results, 
         #  function(x) paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep = "")
   # )
  #) 
  
  # rename columns
  colnames(results) <- c("Risk", "Health outcome", "Spline degree, number of interior knots", "Priors & constraints")
  
  #check if multiple exposures
  if (length(unique(results$`Risk`)) == 1) {
    results <- results[, -c("Risk")]
  }
  
  
  print(ncol(results))
  
  
  # create table
  flextable(results) %>%
    add_header_lines(values = heading) %>%
    add_footer_lines(ifelse(sum(results == "N/A") > 0, "N/A, not available", "")) %>%
    style(part = "header", pr_p = fp_par(text.align = "center")) %>%
    style(j = ifelse(ncol(results) == 4, 1:2, 1), part = "header", pr_p = fp_par(text.align = "left")) %>%
    style(part = "body", pr_p = fp_par(text.align = "center")) %>%
    style(j = ifelse(ncol(results) == 4, 1:2, 1), part = "body", pr_p = fp_par(text.align = "left")) %>%
    width(width = 1.7) %>%
    width(j = ifelse(ncol(results) == 4, 4, 3), width = 2.4) %>%
    font(fontname = "Calibri", part = "all") %>%
    fontsize(size = 10) %>%
    fontsize(i = 1, size = 12, part = "header") %>%
    line_spacing(space = 1.15) %>%
    save_as_docx(path = file.path(output, "parameter_table.docx"), align = "center")
    
  }


##

parameter_table(input, output, heading)
