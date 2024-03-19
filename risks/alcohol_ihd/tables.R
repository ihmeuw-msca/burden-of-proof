##########################################################################################
# Title: A burden of proof study on alcohol consumption and ischemic heart disease
# Purpose: Tables for manuscript (after running burden of proof pipeline)
##########################################################################################

if (F) {
  rm(list = ls())

  ## create customer r library folder
  user <- Sys.getenv("USER")
  user_rlibs <- file.path("filepath")

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
  source("filepath")

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
  #' summary_table(input = "filepath",
  #'               output = "filepath",
  #'               exposure = "alcohol consumption")
  summary_table <- function(input, output, analysis_name, exposure, exp_levels = vector(mode = "double"), footnote = "insert manually") {
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
      inner_ui <- fread(file.path(ro_path, "inner_quantiles.csv"), header = T)
      outer_ui <- fread(file.path(ro_path, "outer_quantiles.csv"), header = T)
      raw_data <- fread(file.path(ro_path, paste0("/raw-", ro_pair, ".csv")))

      if (!all(c(0.025, 0.05, 0.5, 0.95, 0.975) %in% colnames(inner_ui)) &
        !all(c(0.025, 0.05, 0.5, 0.95, 0.975) %in% colnames(outer_ui))) {
        stop("The quantiles were not defined as [0.025, 0.05, 0.5, 0.95, 0.975]")
      }

      # mean RR and 95% UI with between-study heterogeneity
      data <- data.table(
        rr = exp(outer_ui[["0.5"]]),
        lb_i = exp(inner_ui[["0.025"]]),
        ub_i = exp(inner_ui[["0.975"]]),
        lb_o = exp(outer_ui[["0.025"]]),
        ub_o = exp(outer_ui[["0.975"]])
      )

      # determine harmful or protective of the ro_pair
      harmful <- mean(data$rr, na.rm = T) > 1

      # burden of proof risk function
      if (harmful) {
        data$bprf <- exp(outer_ui[["0.05"]])
      } else {
        data$bprf <- exp(outer_ui[["0.95"]])
      }

      # get mean RR and 95% UI
      create_cell <- function(rr, lb, ub) {
        paste0(rr, " (", lb, ", ", ub, ")")
      }

      # RR at select exposure levels
      if (summary$risk_type == "continuous") {
        data[, exposure := outer_ui$risk]

        # at 85th percentile of exposure
        risk_lvl_85 <- round(summary$risk_score_bounds[2])
        rr_fun <- approxfun(data$exposure, data$rr)
        lb_fun <- approxfun(data$exposure, data$lb_o)
        ub_fun <- approxfun(data$exposure, data$ub_o)
        rr_fun_round <- function(x) format(round(rr_fun(x), 2), nsmall = 2)
        lb_fun_round <- function(x) format(round(lb_fun(x), 2), nsmall = 2)
        ub_fun_round <- function(x) format(round(ub_fun(x), 2), nsmall = 2)

        rr_85 <- rr_fun_round(risk_lvl_85)
        lb_85 <- lb_fun_round(risk_lvl_85)
        ub_85 <- ub_fun_round(risk_lvl_85)
        cell_85 <- create_cell(rr_85, lb_85, ub_85)
        risk_lvl_85 <- paste(round(risk_lvl_85, 2), "g/day")

        rr_nadir <- min(data$rr, na.rm = T)
        lb_nadir <- format(round(data[rr == rr_nadir, lb_o], 2), nsmall = 2)
        ub_nadir <- format(round(data[rr == rr_nadir, ub_o], 2), nsmall = 2)
        cell_nadir <- create_cell(format(round(rr_nadir, 2), nsmall = 2), lb_nadir, ub_nadir)
        nadir <- paste(round(data[rr == rr_nadir, exposure]), "g/day")

        # at defined exposure values
        select_rr <- data.table()
        for (i in exp_levels) {
          select_rr[, as.character(i) := create_cell(
            rr_fun_round(i),
            lb_fun_round(i),
            ub_fun_round(i)
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

        nadir <- "N/A"
        cell_nadir <- "N/A"
        risk_lvl_85 <- "N/A"
        cell_85 <- "N/A"

        select_rr <- data.table()
        for (i in exp_levels) {
          select_rr[[i]] <- "N/A"
        }

        # burden of proof risk function
        average_bprf <- round(data[, bprf], 2)
      }

      # average RR
      if (summary$risk_type == "dichotomous") {
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
        cause
      )

      # (average) change in risk
      if (!is.nan(summary$score) && summary$score > 0) {
        average_risk <- paste0(
          round(ifelse(harmful, (average_bprf - 1), (1 - average_bprf)) * 100, 2), "%"
        )
      } else {
        average_risk <- "N/A"
      }

      average_bprf <- format(average_bprf, nsmall = 2)

      # risk-outcome score
      alt_round <- function(x) {
        max(abs(round(x, 2)), abs(signif(x, 1))) * sign(x)
      }
      ros <- alt_round(summary$score)

      # star rating
      thresholds <- c(-Inf, log(1 + c(0.0, 0.15, 0.5, 0.85)), Inf)
      stars <- c(
        emo::ji("star"),
        paste(rep(emo::ji("star"), 2), collapse = ""),
        paste(rep(emo::ji("star"), 3), collapse = ""),
        paste(rep(emo::ji("star"), 4), collapse = ""),
        paste(rep(emo::ji("star"), 5), collapse = "")
      )
      star_rating <- cut(summary$score, thresholds,
        right = FALSE,
        labels = stars, include.lowest = TRUE
      )

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

      # compile results
      temp <- as.data.table(
        cbind(
          risk,
          health_outcome,
          risk_type = summary$risk_type,
          cell_mean_inner,
          cell_mean_outer,
          select_rr,
          nadir,
          cell_nadir,
          risk_lvl_85,
          cell_85,
          average_bprf,
          average_risk,
          ros,
          star_rating,
          pub_bias,
          n_studies
        )
      )
      results <- rbind(results, temp)
    }

    # order by ros
    results[, ros := as.character(ros)]
    results[is.na(ros), ros := "N/A"]

    # Define the desired order of health_outcome
    if (analysis_name == "manuscript") {
      desired_order <- c(
        "cvd_ihd_total",
        "cvd_ihd_incidence",
        "cvd_ihd_females_incidence",
        "cvd_ihd_males_incidence",
        "cvd_ihd_mortality",
        "cvd_ihd_females_mortality",
        "cvd_ihd_males_mortality",
        "cvd_ihd_case_control",
        "cvd_ihd_cohort",
        "cvd_ihd_mr_2sls",
        "cvd_ihd_mi",
        "cvd_ihd_mi_incidence",
        "cvd_ihd_mi_mortality"
      )

      # Change the order of rows in the dataframe
      results <- results[match(desired_order, results$health_outcome), ]

      # rename
      results$health_outcome[results$health_outcome == "cvd_ihd_total"] <- "Ischemic heart disease"
      results$health_outcome[results$health_outcome == "cvd_ihd_incidence"] <- "IHD Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_mortality"] <- "IHD Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi"] <- "Myocardial infarction"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi_incidence"] <- "MI Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi_mortality"] <- "MI Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_females_incidence"] <- "Females – Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_males_incidence"] <- "Males – Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_females_mortality"] <- "Females – Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_males_mortality"] <- "Males – Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_case_control"] <- "Case-control studies"
      results$health_outcome[results$health_outcome == "cvd_ihd_cohort"] <- "Cohort studies"
      results$health_outcome[results$health_outcome == "cvd_ihd_mr_2sls"] <- "Mendelian randomization studies"
    }

    if (analysis_name == "mr_sensitivity") {
      desired_order <- c(
        "cvd_ihd_mr_ivw",
        "cvd_ihd_mr_mvmr",
        "cvd_ihd_mr_nlmr",
        "cvd_ihd_mr_conventional",
        "cvd_ihd_cohort_location"
      )

      # Change the order of rows in the dataframe
      results <- results[match(desired_order, results$health_outcome), ]

      # rename
      results$health_outcome[results$health_outcome == "cvd_ihd_mr_ivw"] <- "Mendelian randomization – IVW"
      results$health_outcome[results$health_outcome == "cvd_ihd_mr_mvmr"] <- "Mendelian randomization – MVMR"
      results$health_outcome[results$health_outcome == "cvd_ihd_mr_nlmr"] <- "Mendelian randomization – non-linear MR"
      results$health_outcome[results$health_outcome == "cvd_ihd_mr_conventional"] <- "Mendelian randomization – Conventional estimates"
      results$health_outcome[results$health_outcome == "cvd_ihd_cohort_location"] <- "Cohort studies in the same location"
    }

    # Define the desired order of health_outcome
    if (analysis_name == "untrimmed") {
      desired_order <- c(
        "cvd_ihd_total",
        "cvd_ihd_incidence",
        "cvd_ihd_females_incidence",
        "cvd_ihd_males_incidence",
        "cvd_ihd_mortality",
        "cvd_ihd_females_mortality",
        "cvd_ihd_males_mortality",
        "cvd_ihd_case_control",
        "cvd_ihd_cohort",
        "cvd_ihd_mi",
        "cvd_ihd_mi_incidence",
        "cvd_ihd_mi_mortality"
      )

      # Change the order of rows in the dataframe
      results <- results[match(desired_order, results$health_outcome), ]

      # rename
      results$health_outcome[results$health_outcome == "cvd_ihd_total"] <- "Ischemic heart disease"
      results$health_outcome[results$health_outcome == "cvd_ihd_incidence"] <- "IHD Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_mortality"] <- "IHD Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi"] <- "Myocardial infarction"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi_incidence"] <- "MI Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_mi_mortality"] <- "MI Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_females_incidence"] <- "Females – Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_males_incidence"] <- "Males – Morbidity"
      results$health_outcome[results$health_outcome == "cvd_ihd_females_mortality"] <- "Females – Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_males_mortality"] <- "Males – Mortality"
      results$health_outcome[results$health_outcome == "cvd_ihd_case_control"] <- "Case-control studies"
      results$health_outcome[results$health_outcome == "cvd_ihd_cohort"] <- "Cohort studies"
    }

    if (unique(results$risk_type) == "continuous") {
      results$risk_type <- NULL
      results$cell_mean_inner <- NULL
      results$cell_mean_outer <- NULL

      # rename columns
      colnames(results) <- c(
        "Risk",
        "Health outcome",
        exp_levels,
        "Nadir exposure level",
        "RR (95% UI) at nadir",
        "85th percentile risk level",
        "RR (95% UI) at 85th percentile risk level",
        "Exposure-averaged BPRF",
        "Conservative interpretation of the average risk increase/decrease",
        "ROS",
        "Star rating",
        "Pub. bias",
        "No. of studies"
      )
    } else {
      # drop irrelevant columns
      results <- results[, c(
        "risk",
        "health_outcome",
        "cell_mean_inner",
        "cell_mean_outer",
        "average_bprf",
        "ros",
        "star_rating",
        "pub_bias",
        "n_studies"
      ),
      with = FALSE
      ]

      # rename columns
      colnames(results) <- c(
        "Risk",
        "Health outcome",
        "RR (95% UI without γ)",
        "RR (95% UI with γ)",
        "BPRF",
        "ROS",
        "Star rating",
        "Pub. bias",
        "No. of studies"
      )
    }

    # check if multiple exposures
    if (length(unique(results$`Risk`)) == 1) {
      results <- results[, -c("Risk")]
    }

    # create table
    flextable(results) %>%
      add_header_lines(values = "Table 1: Strength of the evidence for the relationship between alcohol consumption and ischemic heart disease") %>%
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
      save_as_docx(path = paste0(output, "summary_table_", analysis_name, ".docx"), align = "center")
  }

  # main text
  summary_table(
    input = "filepath",
    output = "filepath",
    exposure = "alcohol consumption",
    analysis_name = "manuscript",
    exp_levels = c(10, 30, 50)
  )

  # sensitivity analysis MR
  summary_table(
    input = "filepath",
    output = "filepath",
    exposure = "alcohol consumption",
    analysis_name = "mr_sensitivity",
    exp_levels = c(10, 30, 50)
  )

  # untrimmed
  summary_table(
    input = "filepath",
    output = "filepath",
    exposure = "alcohol consumption",
    analysis_name = "untrimmed",
    exp_levels = c(10, 30, 50)
  )
}

## bias covariate table
if (F) {
  rm(list = ls())

  ## create customer r library folder
  user <- Sys.getenv("USER")
  user_rlibs <- file.path("filepath")

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
  source("filepath")

  ## main function
  #' This function creates a table for selected bias covariates and gamma solution to be shown in a Burden of Proof study
  #'
  #' @param input is the path to the folder where the results of the BoP pipeline are stored
  #' @param output is the path to the folder where the table should be saved
  #' @param heading optional heading of the table; "Table S12. Bias covariates and estimated parameters" by default.
  #' @param footnote optional footnote of the table; vector(mode = "character") by default.
  #'
  #' @return
  #' @export
  #'
  #' @examples
  #' bc_gamma_table(input = "filepath",
  #'                output = "filepath",
  #'                heading = "Table SX. Selected bias covariates and estimated parameters")
  bc_gamma_table <- function(input, output, heading = "Table SX. Bias covariates and estimated parameters", footnote = vector(mode = "character")) {
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
      gamma <- paste0(round(summary$gamma[1], 2), " (", round(summary$gamma[2], 2), ")")

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

    desired_order <- c(
      "cvd_ihd_total",
      "cvd_ihd_incidence",
      "cvd_ihd_females_incidence",
      "cvd_ihd_males_incidence",
      "cvd_ihd_mortality",
      "cvd_ihd_females_mortality",
      "cvd_ihd_males_mortality",
      "cvd_ihd_case_control",
      "cvd_ihd_cohort",
      "cvd_ihd_mr_2sls",
      "cvd_ihd_mi",
      "cvd_ihd_mi_incidence",
      "cvd_ihd_mi_mortality"
    )

    # Change the order of rows in the dataframe
    results <- results[match(desired_order, results$health_outcome), ]

    # rename
    results$health_outcome[results$health_outcome == "cvd_ihd_total"] <- "Ischemic heart disease"
    results$health_outcome[results$health_outcome == "cvd_ihd_incidence"] <- "IHD Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_mortality"] <- "IHD Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi"] <- "Myocardial infarction"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_incidence"] <- "MI Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_mortality"] <- "MI Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_incidence"] <- "Females – Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_incidence"] <- "Males – Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_mortality"] <- "Females – Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_mortality"] <- "Males – Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_case_control"] <- "Case-control studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_cohort"] <- "Cohort studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_mr_2sls"] <- "Mendelian randomization studies"

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

  bc_gamma_table(
    input = "filepath",
    output = "filepath",
    heading = "Table S12. Bias covariates and estimated parameters"
  )
}

## burden of proof pipeline parameters
if (F) {
  rm(list = ls())

  ## create customer r library folder
  user <- Sys.getenv("USER")
  user_rlibs <- file.path("filepath")

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
  source("filepath")

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
  #' parameter_table(input = "filepath",
  #'                 output = "filepath",
  #'                 heading = "Table SX. MR-BRT spline and prior specifications")
  parameter_table <- function(input, output, heading = "Table SX. MR-BRT spline and prior specifications") {
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

      if (summary$risk_type == "continuous") {
        # Spline degree and interior knots
        spline_degree <- switch(settings$fit_signal_model$cov_model$spline_degree,
          "linear",
          "quadratic",
          "cubic",
          "N/A"
        )
        num_knots <- paste(
          length(settings$fit_signal_model$cov_model$spline_knots) - 2,
          "interior knots"
        )

        degree_knots <- paste(spline_degree, num_knots, sep = ", ")

        # Gaussian prior
        gauss_prior_mean <- settings$fit_signal_model$cov_model$prior_spline_maxder_gaussian[[1]]
        gauss_prior_sd <- settings$fit_signal_model$cov_model$prior_spline_maxder_gaussian[[2]]
        gauss_prior <- paste0(
          "Gaussian max derivative prior on the right tail (",
          gauss_prior_mean[length(gauss_prior_mean)], ", ",
          gauss_prior_sd[length(gauss_prior_sd)], ")"
        )

        # monotonicity constraint
        prior_spline_monotonicity <- settings$fit_signal_model$cov_model$prior_spline_monotonicity

        monotonicity <- ifelse(is.null(prior_spline_monotonicity), NA,
          ifelse(prior_spline_monotonicity == "decreasing",
            "monotonic decreasing", "monotonic increasing"
          )
        )

        # linear tails
        r_linear_tail <- settings$fit_signal_model$cov_model$spline_r_linear
        l_linear_tail <- settings$fit_signal_model$cov_model$spline_l_linear

        if (!r_linear_tail & !l_linear_tail) {
          linear_tail <- NA
        } else if (r_linear_tail & !l_linear_tail) {
          linear_tail <- "right linear tail"
        } else if (!r_linear_tail & l_linear_tail) {
          linear_tail <- "left linear tail"
        } else {
          linear_tail <- "left and right linear tail"
        }

        priors_constraints <- paste0(
          ifelse(is.na(monotonicity), "", paste0(monotonicity, ", ")),
          linear_tail, ", ", gauss_prior
        )
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

    desired_order <- c(
      "cvd_ihd_total",
      "cvd_ihd_incidence",
      "cvd_ihd_females_incidence",
      "cvd_ihd_males_incidence",
      "cvd_ihd_mortality",
      "cvd_ihd_females_mortality",
      "cvd_ihd_males_mortality",
      "cvd_ihd_case_control",
      "cvd_ihd_cohort",
      "cvd_ihd_mr_2sls",
      "cvd_ihd_mi",
      "cvd_ihd_mi_incidence",
      "cvd_ihd_mi_mortality"
    )

    # Change the order of rows in the dataframe
    results <- results[match(desired_order, results$health_outcome), ]

    # rename
    results$health_outcome[results$health_outcome == "cvd_ihd_total"] <- "Ischemic heart disease"
    results$health_outcome[results$health_outcome == "cvd_ihd_incidence"] <- "Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_mortality"] <- "Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi"] <- "Myocardial infarction"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_incidence"] <- "Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_mortality"] <- "Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_incidence"] <- "Females"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_incidence"] <- "Males"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_mortality"] <- "Females"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_mortality"] <- "Males"
    results$health_outcome[results$health_outcome == "cvd_ihd_case_control"] <- "Case-control studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_cohort"] <- "Cohort studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_mr_2sls"] <- "Mendelian randomization studies"

    # make first letters uppercase
    results <- data.table(
      sapply(
        results,
        function(x) paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep = "")
      )
    )


    # rename columns
    colnames(results) <- c("Risk", "Health outcome", "Spline degree, number of interior knots", "Priors & constraints")

    # check if multiple exposures
    if (length(unique(results$`Risk`)) == 1) {
      results <- results[, -c("Risk")]
    }

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

  parameter_table(
    input = "filepath",
    output = "filepath",
    heading = "Table S11. MR-BRT splines and prior specifications"
  )
}

## quantified bias covariates table

if (F) {
  rm(list = ls())
  bundle_id <- 9023
  decomp_step <- "iterative"
  gbd_round_id <- 7
  rei <- 102
  model <- "mr"

  ## load cc functions
  source("filepath")

  if (model == "total") {
    raw <- fread("filepath")
  } else {
    raw_2sls <- fread("filepath")
    raw_nlmr <- fread("filepath")
    raw <- merge(raw_2sls, raw_nlmr, all = T)
  }

  raw <- raw[, -c(
    "ln_rr", "ln_rr_se", "alt_risk_lower", "alt_risk_upper", "ref_risk_lower", "ref_risk_upper",
    "percent_male", "seq", "rei_id", "cause_id", "bundle_id", "bundle_version_id", "risk_type",
    "risk_unit"
  )]
  raw <- unique(raw)
  setnames(raw, "study_id", "nid")

  old_data <- as.data.table(get_bundle_data(bundle_id))
  old_data <- old_data[acause == "cvd_ihd"]
  new_data <- fread("filepath")
  old_data <- unique(old_data[, .(nid, field_citation_value, design)])
  new_data <- unique(new_data[, .(nid, field_citation_value, design)])
  data <- rbind(old_data, new_data)

  data <- merge(raw, data, by = "nid", all.x = TRUE)

  # data <- data[!(nid == 437427)]
  data[, name := sapply(strsplit(field_citation_value, " "), function(x) x[1])]

  cov_columns <- grep("^cov_", names(data), value = TRUE)

  data <- data[, c("name", "design", ..cov_columns)]

  data[name == "GÃ©mes", name := "Gémes"]
  data[name == "de", name := "de Labry"]
  data[name == "Iona", name := "Millwood"]
  data[name == "Joanna", name := "Lankester"]
  data[name == "JunYoung", name := "Chang"]
  data[name == "Philipp", name := "Reddiess"]
  data[name == "Setor", name := "Kuntsor"]
  data[name == "Sudhir", name := "Kurl"]
  data[name == "Rudolph", name := "Schutte"]
  data[name == "Kivela", name := "Kivelä"]
  data[name == "Makela", name := "Makelä"]
  data[name == "Romelsjo", name := "Romelsjö"]
  data[name == "Schroder", name := "Schröder"]

  if (model == "mr") {
    data[is.na(name), name := "Biddinger"]
    data[, design := "Mendelian randomization"]
    data$cov_exposure_selfreport <- NULL
    data$cov_sick_quitters <- "N/A"
    data$cov_design <- "N/A"
  }

  data <- unique(data)

  # order alphabetically
  data <- data[order(name), ]

  data[, design := sub("^(\\w)", "\\U\\1", design, perl = TRUE)]
  setnames(data, "design", "Study design")
  setnames(data, "name", "Author")
  setnames(data, "cov_adjusted_4", "cov_adjusted_3")

  # change column order for mr
  if (model == "mr") {
    data <- data[, c(
      "Author", "Study design", "cov_adjusted_0",
      "cov_adjusted_1", "cov_adjusted_2", "cov_adjusted_3",
      "cov_rep_prevalent_disease", "cov_rep_geography",
      "cov_exposure_study", "cov_outcome_selfreport",
      "cov_sick_quitters", "cov_non_drinker",
      "cov_older", "cov_incidence", "cov_mortality",
      "cov_bmi", "cov_blood_pressure", "cov_cholesterol",
      "cov_apolipoprotein", "cov_fibrinogen", "cov_adiponectin",
      "cov_ihd", "cov_mi", "cov_design"
    )]
  }

  # create table
  flextable(data) %>%
    add_header_lines(values = "Table S7. Quantified bias covariates") %>%
    style(part = "header", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "header", pr_p = fp_par(text.align = "left")) %>%
    style(part = "body", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "body", pr_p = fp_par(text.align = "left")) %>%
    width(width = 0.55) %>%
    width(j = 3:ncol(data), width = 0.45) %>%
    # width(j = "Name", width = 1.0) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 6, part = "all") %>%
    fontsize(size = 6, part = "header") %>%
    line_spacing(space = 1.15) %>%
    save_as_docx(path = ifelse(model == "total", "filepath",
      "filepath"
    ), align = "center")
}


## table with RR at each exposure level
if (F) {
  rm(list = ls())

  ## create customer r library folder
  user <- Sys.getenv("USER")
  user_rlibs <- file.path("filepath")

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

  select_rr_table <- function(input, output, exp_levels) {
    if (!dir.exists(input)) stop("Please provide a valid input path")
    if (!dir.exists(output)) stop("Please provide a valid output path")
    if (!(is.numeric(exp_levels) && length(exp_levels) > 0)) {
      stop("Please provide valid exposure levels")
    }

    results <- data.table()

    # loop through all risk-outcome pairs
    for (ro_pair in list.files(input)) {
      ro_path <- file.path(input, ro_pair)

      summary <- yaml.load_file(file.path(ro_path, "summary.yaml"))
      inner_ui <- fread(file.path(ro_path, "inner_quantiles.csv"), header = T)
      outer_ui <- fread(file.path(ro_path, "outer_quantiles.csv"), header = T)

      # mean RR and 95% UI with between-study heterogeneity
      data <- data.table(
        exposure = outer_ui$risk,
        rr = exp(outer_ui[["0.5"]]),
        lb_i = exp(inner_ui[["0.025"]]),
        ub_i = exp(inner_ui[["0.975"]]),
        lb_o = exp(outer_ui[["0.025"]]),
        ub_o = exp(outer_ui[["0.975"]])
      )

      # get mean RR and 95% UI
      create_cell <- function(rr, lb, ub) {
        paste0(rr, " (", lb, ", ", ub, ")")
      }

      rr_fun <- approxfun(data$exposure, data$rr)
      lb_o_fun <- approxfun(data$exposure, data$lb_o)
      ub_o_fun <- approxfun(data$exposure, data$ub_o)
      lb_i_fun <- approxfun(data$exposure, data$lb_i)
      ub_i_fun <- approxfun(data$exposure, data$ub_i)
      rr_fun_round <- function(x) format(round(rr_fun(x), 2), nsmall = 2)
      lb_o_fun_round <- function(x) format(round(lb_o_fun(x), 2), nsmall = 2)
      ub_o_fun_round <- function(x) format(round(ub_o_fun(x), 2), nsmall = 2)
      lb_i_fun_round <- function(x) format(round(lb_i_fun(x), 2), nsmall = 2)
      ub_i_fun_round <- function(x) format(round(ub_i_fun(x), 2), nsmall = 2)


      # at defined exposure values
      select_rr_o <- list()
      select_rr_i <- list()
      for (i in exp_levels) {
        select_rr_o[[as.character(i)]] <- create_cell(
          rr_fun_round(i),
          lb_o_fun_round(i),
          ub_o_fun_round(i)
        )
        select_rr_i[[as.character(i)]] <- create_cell(
          rr_fun_round(i),
          lb_i_fun_round(i),
          ub_i_fun_round(i)
        )
      }

      health_outcome <- stringr::str_split(ro_pair, "-")[[1]][2]

      # compile results
      temp <- as.data.table(
        cbind(
          health_outcome,
          exp_levels,
          select_rr_o,
          select_rr_i
        )
      )
      results <- rbind(results, temp)
    }

    # desired order
    desired_order <- c(
      "cvd_ihd_total",
      "cvd_ihd_incidence",
      "cvd_ihd_females_incidence",
      "cvd_ihd_males_incidence",
      "cvd_ihd_mortality",
      "cvd_ihd_females_mortality",
      "cvd_ihd_males_mortality",
      "cvd_ihd_case_control",
      "cvd_ihd_cohort",
      "cvd_ihd_mr_2sls",
      "cvd_ihd_mi",
      "cvd_ihd_mi_incidence",
      "cvd_ihd_mi_mortality"
    )

    results[, exp_levels := unlist(exp_levels)]

    results <- results[order(factor(health_outcome, levels = desired_order), exp_levels)]

    results[, select_rr_o := as.character(select_rr_o)]
    results[, select_rr_i := as.character(select_rr_i)]

    results[select_rr_o == "NA (NA, NA)", select_rr_o := "N/A"]
    results[select_rr_i == "NA (NA, NA)", select_rr_i := "N/A"]

    # rename
    results$health_outcome[results$health_outcome == "cvd_ihd_total"] <- "Ischemic heart disease"
    results$health_outcome[results$health_outcome == "cvd_ihd_incidence"] <- "Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_incidence"] <- "Females"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_incidence"] <- "Males"
    results$health_outcome[results$health_outcome == "cvd_ihd_mortality"] <- "Mortality"
    results$health_outcome[results$health_outcome == "cvd_ihd_females_mortality"] <- "Females"
    results$health_outcome[results$health_outcome == "cvd_ihd_males_mortality"] <- "Males"
    results$health_outcome[results$health_outcome == "cvd_ihd_case_control"] <- "Case-control studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_cohort"] <- "Cohort studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_mr_2sls"] <- "Mendelian randomization studies"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi"] <- "Myocardial infarction"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_incidence"] <- "Morbidity"
    results$health_outcome[results$health_outcome == "cvd_ihd_mi_mortality"] <- "Mortality"

    # rename columns
    colnames(results) <- c(
      "Type",
      "Alcohol consumption (g/day)",
      "RR (95% UI with between-study heterogeneity)",
      "RR (95% UI without between-study heterogeneity)"
    )

    # create table
    flextable(results) %>%
      add_header_lines(values = "Table S10. Relative risks across exposure range") %>%
      add_footer_lines("N/A = not available.") %>%
      style(part = "header", pr_p = fp_par(text.align = "center")) %>%
      style(j = 1, part = "header", pr_p = fp_par(text.align = "left")) %>%
      style(part = "body", pr_p = fp_par(text.align = "center")) %>%
      style(j = 1, part = "body", pr_p = fp_par(text.align = "left")) %>%
      width(width = 1.6) %>%
      font(fontname = "Times New Roman", part = "all") %>%
      fontsize(size = 8, part = "all") %>%
      fontsize(size = 10, part = "header") %>%
      line_spacing(space = 1.0) %>%
      save_as_docx(path = paste0(output, "select_rr_table.docx"), align = "center")
  }

  select_rr_table(
    input = "filepath",
    output = "filepath",
    exp_levels = seq(10, 100, by = 10)
  )
}

## table with effect sizes and exposure ranges
if (F) {
  rm(list = ls())
  bundle_id <- 9023

  ## load cc functions
  source("filepath")

  # load final datasets for modeling
  data_obs <- fread("filepath")
  data_obs <- data_obs[, .(
    study_id, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper,
    ln_rr, ln_rr_se
  )]
  data_obs[, design := ifelse(study_id %in% c(501650, 428498, 501522, 501158), "prospective cohort", NA)]

  data_mr_2sls <- fread("filepath")
  data_mr_2sls <- data_mr_2sls[, .(
    study_id, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper,
    ln_rr, ln_rr_se
  )]
  data_mr_2sls$design <- "Mendelian randomization"

  data_mr_ivw <- fread("filepath")
  data_mr_ivw <- data_mr_ivw[, .(
    study_id, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper,
    ln_rr, ln_rr_se
  )]
  data_mr_ivw$design <- "Mendelian randomization"

  data_mr_mvmr <- fread("filepath")
  data_mr_mvmr <- data_mr_mvmr[, .(
    study_id, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper,
    ln_rr, ln_rr_se
  )]
  data_mr_mvmr$design <- "Mendelian randomization"

  data_mr_nlmr <- fread("filepath")
  data_mr_nlmr <- data_mr_nlmr[, .(
    study_id, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper,
    ln_rr, ln_rr_se
  )]
  data_mr_nlmr$design <- "Mendelian randomization"

  data_mr <- merge(data_mr_2sls, data_mr_ivw, all = T)
  data_mr <- merge(data_mr, data_mr_mvmr, all = T)
  data_mr <- merge(data_mr, data_mr_nlmr, all = T)

  # merge datasets
  results <- rbind(data_obs, data_mr, fill = T)

  # load metadata with names
  old_data <- as.data.table(get_bundle_data(bundle_id))
  old_data <- old_data[acause == "cvd_ihd"]
  new_data <- fread("filepath")
  old_data <- unique(old_data[, .(nid, field_citation_value, design)])
  new_data <- unique(new_data[, .(nid, field_citation_value, design)])
  metadata <- rbind(old_data, new_data)

  # data <- merge(raw, data, by = "nid", all.x = TRUE)

  # data <- data[!(nid == 437427)]
  metadata[, name := sapply(strsplit(field_citation_value, " "), function(x) x[1])]

  metadata[name == "GÃ©mes", name := "Gémes"]
  metadata[name == "de", name := "de Labry"]
  metadata[name == "Iona", name := "Millwood"]
  metadata[name == "Joanna", name := "Lankester"]
  metadata[name == "JunYoung", name := "Chang"]
  metadata[name == "Philipp", name := "Reddiess"]
  metadata[name == "Setor", name := "Kuntsor"]
  metadata[name == "Sudhir", name := "Kurl"]
  metadata[name == "Rudolph", name := "Schutte"]
  metadata[name == "Kivela", name := "Kivelä"]
  metadata[name == "Makela", name := "Makelä"]
  metadata[name == "Romelsjo", name := "Romelsjö"]
  metadata[name == "Schroder", name := "Schröder"]

  setnames(metadata, "nid", "study_id")

  metadata$field_citation_value <- NULL

  # merge with names
  results <- merge(results, metadata, by = "study_id", all.x = T)

  # order alphabetically
  results <- results[order(name), .(name, design.y, ref_risk_lower, ref_risk_upper, alt_risk_lower, alt_risk_upper, ln_rr, ln_rr_se)]

  # round columns to one or two decimal places
  cols <- names(results)[3:6]
  results[, (cols) := round(.SD, 1), .SDcols = cols]

  cols <- names(results)[7:8]
  results[, (cols) := round(.SD, 2), .SDcols = cols]

  # remove .0s
  results[, ref_risk_lower := str_replace(ref_risk_lower, "//.0", "")]
  results[, ref_risk_upper := str_replace(ref_risk_upper, "//.0", "")]
  results[, alt_risk_lower := str_replace(alt_risk_lower, "//.0", "")]
  results[, alt_risk_upper := str_replace(alt_risk_upper, "//.0", "")]

  # combine reference and alternative exposure groups
  results[, reference := paste(ref_risk_lower, ref_risk_upper, sep = " - ")]
  results[, alternative := paste(alt_risk_lower, alt_risk_upper, sep = " - ")]

  # remove duplicate 0s in reference
  results[reference == "0-0", reference := 0]

  # make design uppercase
  results[, design.y := sub("^(\\w)", "\\U\\1", design.y, perl = TRUE)]


  # paste grams/day
  results[, reference := paste0(reference, " g/day")]
  results[, alternative := paste0(alternative, " g/day")]

  # subset to relevant variables
  results <- results[, .(name, design.y, reference, alternative, ln_rr, ln_rr_se)]

  # change column names
  setnames(results, "name", "Author")
  setnames(results, "design.y", "Study design")
  setnames(results, "alternative", "Alternative exposure group")
  setnames(results, "reference", "Reference exposure group")
  setnames(results, "ln_rr", "Log effect size")
  setnames(results, "ln_rr_se", "Standard error of effect size")

  # create table
  flextable(results) %>%
    add_header_lines(values = "Table S9. Results from input studies") %>%
    style(part = "header", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "header", pr_p = fp_par(text.align = "left")) %>%
    style(part = "body", pr_p = fp_par(text.align = "center")) %>%
    style(j = 1, part = "body", pr_p = fp_par(text.align = "left")) %>%
    width(width = 1.0) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 8, part = "all") %>%
    fontsize(size = 10, part = "header") %>%
    line_spacing(space = 1.0) %>%
    save_as_docx(path = "filepath", align = "center")
}