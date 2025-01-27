#######################################################################    
# Title: Create table for selected bias covariates and gamma solution for a BoP study
# Author: 
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
source("FILEPATH/r/get_ids.R")

## main function
#' This function creates a table for selected bias covariates and gamma solution to be shown in a Burden of Proof study
#'
#' @param input is the path to the folder where the results of the BoP pipeline are stored
#' @param output is the path to the folder where the table should be saved
#' @param heading optional heading of the table; "Table SX. Selected bias covariates and gamma solution" by default.
#' @param footnote optional footnote of the table; vector(mode = "character") by default.
#'
#' @return
#' @export
#'
#' @examples 
#' bc_gamma_table(input = "/homes/shcarr/bop_pipeline/results/dichotomous/",
#'                output = "/homes/shcarr/bop_pipeline/tables/", 
#'                heading = "Table S4. Selected bias covariates and gamma solution")


risk_folder = "processed_meat"
trimming = "no_trimming/"
exposure = paste0(risk_folder, "consumption")

input = paste0("FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/", trimming, "/risk_outcome_pair/")
output = paste0("FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/", trimming, "/summary_table/")
heading = paste0("Table:", exposure ,"  Selected bias covariates and gamma solution")


bc_gamma_table <- function(input, output, heading = "Table SX. Selected bias covariates and gamma solution", footnote = vector(mode = "character")){
  
  if (!dir.exists(input)) stop("Please provide a valid input path")
  if (!dir.exists(output)) stop("Please provide a valid output path")
  
  # load cause_ids
  cause_ids <- get_ids("cause")
  
  results <- data.table()
  
  # loop through all risk-outcome pairs
  for (ro_pair in list.files(input)) {
    ro_path <- file.path(input, ro_pair)
    
    # load files
    summary <- yaml.load_file(file.path(ro_path, "summary.yaml"))
    cov_finder_result <- yaml.load_file(file.path(ro_path, "cov_finder_result.yaml"))
    
    # bias covariate(s)
    selected_bc_covs <- cov_finder_result$selected_covs
    selected_bc_covs <- gsub("cov_", "", selected_bc_covs)
    selected_bc_covs <- ifelse(length(selected_bc_covs) > 0, paste(selected_bc_covs, collapse = ", "), "None")
    
    # gamma solution
    gamma <- paste0(round(summary$gamma[1], 10), " (", round(summary$gamma[2], 5), ")")
    
    # exposure
    risk <- stringr::str_split(ro_pair, "-")[[1]][1]
    
    # health outcome
    cause <- stringr::str_split(ro_pair, "-")[[1]][2]
    health_outcome <- ifelse(cause %in% cause_ids$acause, cause_ids$cause_name[cause_ids$acause == cause], cause)
    
    # compile results
    temp <- data.frame(
      cbind(
        risk,
        health_outcome,
        selected_bc_covs,
        gamma
      )
    )
    results <- rbind(results, temp)
  }
  
  # rename columns
  colnames(results) <- c("Risk", "Health outcome", "Selected bias covariates", "Gamma solution (mean and sd)")
  
  # check if multiple exposures
  if (length(unique(results$`Risk`)) == 1) {
    results <- results[, -c("Risk")]
  } 
  
  # define table settings
  if ("Risk" %in% colnames(results)) {
    cols <- c(1:4)
    width <- c(1.8, 1.8, 1.8, 1.3)
  } else {
    cols <- c(1:3)
    width <- c(1.8, 1.8, 1.3)
  }
  
  # create table
  flextable(results) %>%
    add_header_lines(values = heading) %>%
    add_footer_lines(footnote) %>%
    width(j = cols, width = width) %>%
    style(part = "header", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "header", pr_p = fp_par(text.align = "left")) %>%
    style(part = "body", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "body", pr_p = fp_par(text.align = "left")) %>%
    font(fontname = "Calibri", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(i = 1, size = 12, part = "header") %>%
    line_spacing(space = 1.15) %>%
    save_as_docx(path = file.path(output, "bc_gamma_table.docx"), align = "center")
}

## run the function


bc_gamma_table(input,output, heading)

