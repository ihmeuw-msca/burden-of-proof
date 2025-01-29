#######################################################################    
# Title: BoP Summary Table
#######################################################################    
rm(list = ls())

## create customer r library folder
user <- Sys.getenv("USER")
user_rlibs <- file.path("/homes", user, "rlibs")

if (!dir.exists(user_rlibs)) {
  dir.create(user_rlibs)
} 

## install/load packages
packages <- c("yaml", "data.table", "stringr", "flextable", "officer", "english")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, lib = user_rlibs)
    library(p, lib.loc = user_rlibs, character.only = TRUE)
  }
}

## install/load emo package
if (!require(emo, lib.loc = user_rlibs)) {
  devtools::install_github("hadley/emo", lib = user_rlibs)
  library(emo, lib.loc = user_rlibs)
}

## load cc functions
source("FILEPATH/r/get_ids.R")

## main function
#' This function creates a table for key results to be shown in a Burden of Proof study
#'
#' @param input is the path to the folder where the results of the BoP pipeline are stored
#' @param output is the path to the folder where the table should be saved
#' @param exposure is the name of the exposure of interest
#' @param exp_levels optional vector of exposure levels for continuous exposure for which relative risks are to be calculated; vector(mode = "double") by default.
#' @param footnote optional footnote of the table; standard footnote by default.
#'
#' @return
#' @export
#'
#' @examples
#' summary_table(input = "/homes/user/bop_pipeline/results/dichotomous/",
#'              output = "/homes/user/bop_pipeline/tables/", 
#'             exposure = "alcohol consumption")
 
 risk_folder = "SSB"
 trimming = "no_trimming/"
 input =  paste0("FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/all_data/", trimming, "/risk_outcome_pair/")
 output = paste0("FILEPATH/", risk_folder, "/mrbrt_run/BoP_summary_table/all_data/", trimming, "/summary_table/")
 exposure = paste0(risk_folder, "consumption")
 

summary_table <- function(input, output, exposure, exp_levels = vector(mode = "double"), footnote = paste0('The reported relative risk (RR) and its 95% uncertainty interval (UI) reflect the risk an individual who has been exposed to ', exposure, ' has of developing the outcome of interest relative to that of someone who has not been exposed to ', ifelse('and' %in% exposure, 'these risk factors', exposure), '. Gamma (γ) is the estimated between-study heterogeneity. We report the 95% UI when not incorporating between-study heterogeneity (γ) \U2212 "95% UI without γ" \U2212 and when accounting for between-study heterogeneity \U2212 "95% UI with γ." The Burden of Proof Risk Function (BPRF) is calculated for risk-outcome pairs that were found to have significant relationships at an 0.05 level of significance when not incorporating between-study heterogeneity (i.e., the lower bound of the 95% UI without γ does not cross the null RR value of one). The BPRF corresponds to the 5th quantile estimate of relative risk accounting for between-study heterogeneity closest to the null for each risk–outcome pair, and it reflects the most conservative estimate of excess risk associated with ', exposure, ' that is consistent with the available data. Since we define ', exposure,' exposure as dichotomous risk factors, i.e., an individual either has been exposed or has not, the risk-outcome score (ROS) is calculated as the signed value of log(BPRF) divided by two. Negative ROSs indicate that the evidence of the association is very weak and inconsistent. For ease of interpretation, we have transformed the ROS and BPRF into a star rating (0-5) with a higher rating representing a larger effect with stronger evidence. The potential existence of publication bias, which, if present, would affect the validity of the results, was tested using Egger', "\'", 's Regression. Included studies represent all available relevant data identified through our systematic reviews from January 1970 through January 2023. The selected bias covariates were chosen for inclusion in the model using an algorithm that systematically detects bias covariates that correspond to significant sources of bias in the observations included. If selected, the observations were adjusted to better reflect the gold standard values of the covariate. See the Supplementary Information for more information about the candidate bias covariates that were selected for in each model.')) {
  
  if (!dir.exists(input)) stop("Please provide a valid input path")
  if (!dir.exists(output)) stop("Please provide a valid output path")
  if (!(is.character(exposure) && nchar(exposure) > 0)) {
    stop("Please provide a valid exposure")
  }
  
  # load cause_ids
  cause_ids <- get_ids("cause")
  
  results <- data.table()
  
  # loop through all risk-outcome pairs
  for (ro_pair in list.files(input)) {
    
    ro_path <- file.path(input, ro_pair)
    
    # load files
    summary <- yaml.load_file(file.path(ro_path, "summary.yaml"))
    cov_finder_result <- yaml.load_file(file.path(ro_path, "cov_finder_result.yaml"))
    inner_ui <- fread(file.path(ro_path, "inner_quantiles.csv"), header = T)
    outer_ui <- fread(file.path(ro_path, "outer_quantiles.csv"), header = T)
    raw_data <- fread(file.path(ro_path, paste0("/raw-", ro_pair, ".csv")))
    
    if (!all(c(0.025, 0.05, 0.5, 0.95, 0.975) %in% colnames(inner_ui)) &
        !all(c(0.025, 0.05, 0.5, 0.95, 0.975) %in% colnames(outer_ui))) {
      stop("The quantiles were not defined as [0.025, 0.05, 0.5, 0.95, 0.975]")
    }
    
    # mean RR and 95% UI with between-study heterogeneity
    data <- data.table(
      rr = exp(outer_ui[['0.5']]),
      lb_i = exp(inner_ui[['0.025']]),
      ub_i = exp(inner_ui[['0.975']]),
      lb_o = exp(outer_ui[['0.025']]),
      ub_o = exp(outer_ui[['0.975']])
    )
    
    # determine harmful or protective of the ro_pair
    harmful <- mean(data$rr) > 1
    
    # burden of proof risk function
    if (harmful) {
      data$bprf <- exp(outer_ui[['0.05']])
    } else {
      data$bprf <- exp(outer_ui[['0.95']])
    }
    
    # get mean RR and 95% UI without gamma
    create_cell <- function(rr, lb, ub) {
      paste0(rr, " (", lb, ", ", ub, ")")

    }
    
    create_cell_i <- function(rr,lbi,ubi) {
      paste0(rr, " (", lbi, ", ", ubi, ")")
      
    }
    
    # RR at select exposure levels
    if(summary$risk_type == "continuous"){
      
      data[, exposure := outer_ui$risk]
      
      # at 85th percentile of exposure
      risk_lvl_85 <- summary$risk_score_bounds[2]
      rr_fun <- approxfun(data$exposure, data$rr)
      lb_fun <- approxfun(data$exposure, data$lb_o)
      ub_fun <- approxfun(data$exposure, data$ub_o)
      
      lbi_fun <- approxfun(data$exposure, data$lb_i)
      ubi_fun <- approxfun(data$exposure, data$ub_i)
      
      rr_fun_round <- function(x) round(rr_fun(x), 2)
      lb_fun_round <- function(x) round(lb_fun(x), 2)
      ub_fun_round <- function(x) round(ub_fun(x), 2)
      
      lbi_fun_round <- function(x) round(lbi_fun(x), 2)
      ubi_fun_round <- function(x) round(ubi_fun(x), 2)
      
      
      rr_85 <- rr_fun_round(risk_lvl_85)
      lb_85 <- lb_fun_round(risk_lvl_85)
      ub_85 <- ub_fun_round(risk_lvl_85)
      
      lbi_85 <- lbi_fun_round(risk_lvl_85)
      ubi_85 <- ubi_fun_round(risk_lvl_85)
      
      cell_85 <- create_cell(rr_85, lb_85, ub_85)
      cell_85_i <- create_cell_i(rr_85, lbi_85, ubi_85)
      risk_lvl_85 <- round(risk_lvl_85, 2)
      
      # at defined exposure values
      select_rr <- data.table()
      for(i in exp_levels){
        select_rr[, as.character(i) := create_cell(rr_fun_round(i), 
                                                   lb_fun_round(i), 
                                                   ub_fun_round(i),
                                                   lbi_fun_round(i), 
                                                   ubi_fun_round(i),
                                                   
                                                   )]
      }
      
      # average bprf between 15th and 85th exposure level
      average_bprf <- mean(data[
        (exposure >= summary$risk_score_bounds[1]) & (exposure <= summary$risk_score_bounds[2]),
        bprf
      ])
      average_bprf <- round(average_bprf, 2)
      
    } else {
      
      risk_unit <- "dichotomous"
      
      risk_lvl_85 <- "N/A"
      cell_85 <- "N/A"
      
      select_rr <- data.table()
      for(i in exp_levels){
        select_rr[[i]] <- "N/A"
      }
      
      # burden of proof risk function
      average_bprf <- round(data[, bprf], 2)
      
    }
    
    # average RR
    if(summary$risk_type == "dichotomous"){
      cell_mean_inner <- create_cell(round(data[, rr], 2), round(data[, lb_i], 2), round(data[, ub_i], 2))
      cell_mean_outer <- create_cell(round(data[, rr], 2), round(data[, lb_o], 2), round(data[, ub_o], 2))
    } else {
      cell_mean_inner <- NA
      cell_mean_outer <- NA
    }
    
    # exposure
    risk <- stringr::str_split(ro_pair, "-")[[1]][1]
    
    # health outcome
    cause <- stringr::str_split(ro_pair, "-")[[1]][2]
    health_outcome <- ifelse(cause %in% cause_ids$acause, 
                             cause_ids$cause_name[cause_ids$acause == cause], 
                             cause)
    
    # (average) change in risk
    if (!is.na(summary$score) & summary$score > 0) {
      average_risk <- paste0(
        round(ifelse(harmful, (average_bprf - 1), (1 - average_bprf)) * 100, 2), "%"
      )
    } else {
      average_risk <- "N/A"
    }
    
    # risk-outcome score
    alt_round <- function(x) {
      max(abs(round(x, 2)), abs(signif(x, 1))) * sign(x)
    }
    ros <- alt_round(summary$score)   
    
    # star rating
    thresholds <- c(-Inf, log(1 + c(0.0, 0.15, 0.5, 0.85)), Inf)
    stars <- c(emo::ji("star"),
               paste(rep(emo::ji("star"), 2), collapse = ""),
               paste(rep(emo::ji("star"), 3), collapse = ""),
               paste(rep(emo::ji("star"), 4), collapse = ""),
               paste(rep(emo::ji("star"), 5), collapse = ""))
    star_rating <- cut(summary$score, thresholds, right = FALSE,
                       labels = stars, include.lowest = TRUE)
    
    # when the fixed effects uncertainty includes the null
    beta_qs <- qnorm(c(0.025, 0.975), mean = summary$beta[1], sd = summary$beta[2])
    if (min(sign(summary$beta[1]) * beta_qs) < 0) {
      average_bprf <- "N/A"
      ros <- NA
      star_rating <- ""
    }
    
    # publication/reporting bias
    pub_bias <- ifelse(summary$pub_bias == 1, "Yes", "No")
    
    # number of studies
    n_studies <- length(unique(raw_data$study_id))
    
    # bias covariate(s)
    selected_bc_covs <- cov_finder_result$selected_covs
    selected_bc_covs <- gsub("bc_", "", selected_bc_covs)
    selected_bc_covs <- ifelse(length(selected_bc_covs) > 0, paste(selected_bc_covs, collapse = ", "), "None")
    
    # compile results
    temp <- as.data.table(
      cbind(
        risk, 
        health_outcome,
        risk_type = summary$risk_type,
        cell_mean_inner,
        cell_mean_outer,
        select_rr,
        risk_lvl_85,
        cell_85,
        cell_85_i,
        average_bprf,
        average_risk,
        ros,
        star_rating,
        pub_bias,
        n_studies,
        selected_bc_covs
      )
    )
    results <- rbind(results, temp)
  }
  
  # order by ros
  results <- results[order(-ros, na.last = TRUE), ]
  results[, ros := as.character(ros)]
  results[is.na(ros), ros := "N/A"]
  

  if(unique(results$risk_type) == "continuous") {
    results$risk_type <- NULL
    results$cell_mean_inner <- NULL
    results$cell_mean_outer <- NULL
    
    # rename columns
    colnames(results) <- c("Risk", 
                           "Health outcome", 
                           exp_levels, 
                           "85th percentile risk level", 
                           "RR (95% UI with γ)", 
                           "RR (95% UI without γ)", 
                           "Exposure-averaged BPRF", 
                           "Conservative interpretation of the average risk increase/decrease", 
                           "ROS", 
                           "Star rating", 
                           "Pub. bias", 
                           "No. of studies",
                           "Selected bias covariates")
       #dhw
    
    
  } else {
    
    # drop irrelevant columns
    results <- results[, c("risk", 
                           "health_outcome", 
                           "cell_mean_inner",
                           "cell_mean_outer",
                           "average_bprf", 
                           "ros", 
                           "star_rating", 
                           "pub_bias", 
                           "n_studies",
                           "selected_bc_covs"), 
                       with = FALSE]
    
    # rename columns
    colnames(results) <- c("Risk", 
                           "Health outcome", 
                           "RR (95% UI without γ)", 
                           "RR (95% UI with γ)", 
                           "BPRF", 
                           "ROS", 
                           "Star rating", 
                           "Pub. bias", 
                           "No. of studies",
                           "Selected bias covariates")
    
  }
  
  # check if multiple exposures
  if (length(unique(results$`Risk`)) == 1) {
    results <- results[, -c("Risk")]
  }  
  
  
  # create table
  
  ## Doc margin,
  
  sect_properties <- prop_section(
    page_size = page_size(orient = "landscape",
                          width = 8.3, height = 11.7),
    type = "continuous",
    page_margins = page_mar()
  )
  
  flextable(results) %>%
    add_header_lines(values = 
                       paste0("Table 2. Strength of the evidence for the relationship between ", 
                              exposure, " and the ", 
                              ifelse(nrow(results) > 1, 
                                     paste0(as.character(english(nrow(results))), " health outcomes analyzed"), 
                                     " health outcome analyzed."))) %>%
    add_footer_lines(footnote) %>%
    style(part = "header", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "header", pr_p = fp_par(text.align = "left")) %>%
    style(part = "body", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "body", pr_p = fp_par(text.align = "left")) %>%
    width(width = 0.8) %>%
    width(j = "Health outcome", width = 1.0) %>%
    font(fontname = "Calibri", part = "all") %>%
    fontsize(size = 8, part = "all") %>%
    fontsize(size = 10, part = "header") %>%
    line_spacing(space = 1.15) %>%
    
    save_as_docx(path = file.path(output, "summary_table.docx"), align = "center",  pr_section = sect_properties)
}


## Run the Function: Added by Demewoz


summary_table(input = input, output = output,exposure = exposure)

