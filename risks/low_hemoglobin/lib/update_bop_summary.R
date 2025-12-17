# Function to update summary.yaml from BoP pipeline ----
#'
#' This function processes directories containing BOP analysis outputs,
#' specifically by applying various calculations and formatting to generate
#' summaries. It performs BoP output data extraction, risk score calculations,
#' relative risk (RR) estimations, and results formatting to facilitate
#' interpretation and reporting.
#'
#' @param parent_dir A directory path containing BoP results, which are used to
#'   write updated summary information to a `.yaml` file in that directory.
#' @param split_at A single numeric value or the string "tmrel" indicating the
#'   value at which the risk should be split. Default is "tmrel".   
#' @param split_method A character string to indicate desired method for
#'   splitting at TMREL. Default is "simple" and alternate option is "points".
#'
#' @export
update_summary <- function(parent_dir,
                           split_at = "tmrel",
                           split_method = "simple") {
  
  # Validate inputs
  checkmate::assert(
    checkmate::check_choice(split_at, "tmrel"), 
    checkmate::check_number(split_at, null.ok = TRUE)
  )
  checkmate::assert_subset(split_method, choices = c("simple", "points"))
  
  # For each trial directory, write an updated summary file
  trial_dirs <- list.dirs(parent_dir, full.names = TRUE, recursive = TRUE)
  #trial_dirs <- trial_dirs[grepl("low_", trial_dirs)]

  for (dir in trial_dirs) {

    tryCatch({
      # Load necessary objects and data
      study_data <- data.table::fread(list.files(dir, pattern = "^raw.*\\.csv$", full.names = TRUE)[1])

      output <- tryCatch(
        {
          # First attempt with split_at argument value
          if (split_at == "tmrel") {
            ref_risk <- NULL
          } else {
            ref_risk <- split_at
          }
          
          bophf::bop(
            directory = dir,
            split_at = split_at,
            reference_risk = ref_risk,
            nll = TRUE,
            split_method = split_method
          )
        },
        error = function(e) {
          # Fallback if the first attempt fails
          warning("Attempt to split at ", split_at, " failed with error: ",
                  e$message, "\nRetrying without split_at argument...")
          bophf::bop(
            directory = dir,
            reference_risk = NULL,
            nll = TRUE,
            split_method = split_method
          )
        }
      )

      all_split_at <- output[grep("^split_at$", names(output))]
      split_at_value_present <- any(
        sapply(all_split_at, function(x) split_at %in% x)
      )
      
      # For use of identifying "tmrel" in output$split_at
      if (is.null(split_at) | split_at == "tmrel") {
        split_at <- output$tmrel
      }

      # Set risk score bounds
      risk_15th_pct_overall <- output$overall$risk_score_bounds[1]
      risk_85th_pct_overall <- output$overall$risk_score_bounds[2]

      if (split_at_value_present) {
        risk_15th_pct_lower <- output$lower$risk_score_bounds[1]
        risk_85th_pct_lower <- output$lower$risk_score_bounds[2]
        risk_15th_pct_upper <- output$upper$risk_score_bounds[1]
        risk_85th_pct_upper <- output$upper$risk_score_bounds[2]
      } else {
        message(
          "Skipping setting the risk score bounds for lower & upper: \n",
          "bophf::bop was NOT run with split_at = ",
          split_at
        )
      }

      # Knot bounds
      output$knot_bounds <- get_best_spline_knots(dir)

      # RMSE
      output$rmse <- get_rmse(
        quantiles_name = "inner_quantiles",
        study_data = study_data,
        bop = output
      )

      # Derive mean BPRF
      ## Overall (15th-85th pct)
      outer_quantiles_15th_to_85th_pct_overall <- output$outer_quantiles[output$outer_quantiles$risk >= risk_15th_pct_overall & output$outer_quantiles$risk <= risk_85th_pct_overall, ]
      output$mean_bprf_overall <- exp(mean(outer_quantiles_15th_to_85th_pct_overall$lower_bprf, na.rm = TRUE))
      output$overall$mean_bprf <- output$mean_bprf_overall

      if (split_at_value_present) {
        ## Lower (15th-85th pct below reference risk value)
        outer_quantiles_15th_to_85th_pct_below_ref <- output$outer_quantiles[output$outer_quantiles$risk >= risk_15th_pct_lower & output$outer_quantiles$risk <= risk_85th_pct_lower, ]
        output$mean_bprf_lower <- exp(mean(outer_quantiles_15th_to_85th_pct_below_ref$lower_bprf, na.rm = TRUE))
        output$lower$mean_bprf <- output$mean_bprf_lower

        ## Upper (15th-85th pct above reference risk value)
        outer_quantiles_15th_to_85th_pct_above_ref <- output$outer_quantiles[output$outer_quantiles$risk >= risk_15th_pct_upper & output$outer_quantiles$risk <= risk_85th_pct_upper, ]
        output$mean_bprf_upper <- exp(mean(outer_quantiles_15th_to_85th_pct_above_ref$lower_bprf, na.rm = TRUE))
        output$upper$mean_bprf <- output$mean_bprf_upper

      } else {
        message(
          "Skipping BPRF calculation for lower & upper: \n",
          "bophf::bop was NOT run with split_at = ",
          split_at
        )
      }

      # mean RR & UI by overall (15th-85th pct overall), lower (15th-85th pct below reference risk value) & upper (15th-85th pct above reference risk value), fe
      overall_range <- c(risk_15th_pct_overall, risk_85th_pct_overall)
      output$mean_rr_overall_fe <- mean_rr_for(output$inner_quantiles, overall_range)

      if (split_at_value_present) {
        lower_range <- c(risk_15th_pct_lower, risk_85th_pct_lower)
        output$mean_rr_lower_fe <- mean_rr_for(output$inner_quantiles, lower_range)

        upper_range <- c(risk_15th_pct_upper, risk_85th_pct_upper)
        output$mean_rr_upper_fe <- mean_rr_for(output$inner_quantiles, upper_range)

      } else {
        message(
          "Skipping RR & UI calculation for lower & upper: \n",
          "bophf::bop was NOT run with split_at = ",
          split_at
        )
      }

      # Derive estimates
      ## For non-anemic, fe
      output$mean_rr_non_anemic_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 110)

      ## For non-anemic, re
      output$mean_rr_non_anemic_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 110)

      ## For mild anemia, fe
      output$mean_rr_mild_anemia_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 100, upper_risk_value = 110)

      ## For mild anemia, re
      output$mean_rr_mild_anemia_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 100, upper_risk_value = 110)

      ## For moderate anemia, fe
      output$mean_rr_moderate_anemia_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 70, upper_risk_value = 100)

      ## For moderate anemia, re
      output$mean_rr_moderate_anemia_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 70, upper_risk_value = 100)

      ## For severe anemia, fe
      output$mean_rr_severe_anemia_fe <- calculate_rr_output(output, "inner_quantiles", upper_risk_value = 70)

      ## For severe anemia, re
      output$mean_rr_severe_anemia_re <- calculate_rr_output(output, "outer_quantiles", upper_risk_value = 70)

      ## At anemia cutoff, fe
      output$rr_at_anemia_cutoff_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 110, upper_risk_value = 110)

      ## At anemia cutoff, re
      output$rr_at_anemia_cutoff_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 110, upper_risk_value = 110)

      ## 110 to 120 g/L, fe
      output$mean_rr_110_to_120_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 110, upper_risk_value = 120)

      ## 110 to 120 g/L, re
      output$mean_rr_110_to_120_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 110, upper_risk_value = 120)

      ## At 115 g/L, fe
      output$rr_at_115_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 115, upper_risk_value = 115)

      ## At 115 g/L, re
      output$rr_at_115_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 115, upper_risk_value = 115)

      ## At 120 g/L, fe
      output$rr_at_120_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 120, upper_risk_value = 120)

      ## At 120 g/L, re
      output$rr_at_120_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 120, upper_risk_value = 120)

      ## At 125 g/L, fe
      output$rr_at_125_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 125, upper_risk_value = 125)

      ## At 125 g/L, re
      output$rr_at_125_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 125, upper_risk_value = 125)

      ## 110 to reference risk value, fe
      output$mean_rr_110_to_ref_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 110, upper_risk_value = split_at)

      ## 110 to reference risk value, re
      output$mean_rr_110_to_ref_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 110, upper_risk_value = split_at)

      ## reference risk value to 140 g/L, fe
      output$mean_rr_ref_to_140_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = split_at, upper_risk_value = 140)

      ## reference risk value to 140 g/L, re
      output$mean_rr_ref_to_140_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = split_at, upper_risk_value = 140)

      ## reference risk value to 150 g/L, fe
      output$mean_rr_ref_to_150_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = split_at, upper_risk_value = 150)

      ## reference risk value to 150 g/L, re
      output$mean_rr_ref_to_150_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = split_at, upper_risk_value = 150)

      ## reference risk value to 160 g/L, fe
      output$mean_rr_ref_to_160_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = split_at, upper_risk_value = 160)

      ## reference risk value to 160 g/L, re
      output$mean_rr_ref_to_160_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = split_at, upper_risk_value = 160)

      ## reference risk value to 85th percentile, fe
      output$mean_rr_ref_to_85th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = split_at, upper_risk_value = output$overall$risk_score_bounds[2])

      ## reference risk value to 85th percentile, re
      output$mean_rr_ref_to_85th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = split_at, upper_risk_value = output$overall$risk_score_bounds[2])

      ## reference risk value to 100th percentile, fe
      output$mean_rr_ref_to_100th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = split_at, upper_risk_value = output$risk_bounds[2])

      ## reference risk value to 100th percentile, re
      output$mean_rr_ref_to_100th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = split_at, upper_risk_value = output$risk_bounds[2])

      ## >=140 g/L, fe
      output$mean_rr_over_140_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 140)

      ## >=140 g/L, re
      output$mean_rr_over_140_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 140)

      ## >=150 g/L, fe
      output$mean_rr_over_150_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 150)

      ## >=150 g/L, re
      output$mean_rr_over_150_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 150)

      ## >=160 g/L, fe
      output$mean_rr_over_160_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = 160)

      ## >=160 g/L, re
      output$mean_rr_over_160_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = 160)

      ## For 0th to 15th percentiles, fe
      output$mean_rr_0th_to_15th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = output$risk_bounds[1], upper_risk_value = output$overall$risk_score_bounds[1])

      ## For 0th to 15th percentiles, re
      output$mean_rr_0th_to_15th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = output$risk_bounds[1], upper_risk_value = output$overall$risk_score_bounds[1])

      ## For 15th percentile, fe
      output$rr_15th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = output$overall$risk_score_bounds[1], upper_risk_value = output$overall$risk_score_bounds[1])

      ## For 15th percentile, re
      output$rr_15th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = output$overall$risk_score_bounds[1], upper_risk_value = output$overall$risk_score_bounds[1])

      ## For 15th percentile to reference value, fe
      output$mean_rr_15th_pct_to_ref_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = output$overall$risk_score_bounds[1], upper_risk_value = split_at)

      ## For 15th percentile to reference value, re
      output$mean_rr_15th_pct_to_ref_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = output$overall$risk_score_bounds[1], upper_risk_value = split_at)

      ## For 85th percentile, fe
      output$rr_85th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = output$overall$risk_score_bounds[2], upper_risk_value = output$overall$risk_score_bounds[2])

      ## For 85th percentile, re
      output$rr_85th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = output$overall$risk_score_bounds[2], upper_risk_value = output$overall$risk_score_bounds[2])

      ## For 85th to 100th percentiles, fe
      output$mean_rr_85th_to_100th_pct_fe <- calculate_rr_output(output, "inner_quantiles", lower_risk_value = output$overall$risk_score_bounds[2], upper_risk_value = output$risk_bounds[2])

      ## For 85th to 100th percentiles, re
      output$mean_rr_85th_to_100th_pct_re <- calculate_rr_output(output, "outer_quantiles", lower_risk_value = output$overall$risk_score_bounds[2], upper_risk_value = output$risk_bounds[2])

      # Total unique NID & observation count
      output$unique_nids_total <- as.integer(length(unique(study_data$study_id)))
      output$obs_total <- as.integer(nrow(study_data))

      # Unique NID & observation count for low hb
      low_hb_rows <- which(
        study_data$alt_risk_upper <= study_data$ref_risk_lower &
          study_data$alt_risk_lower < study_data$ref_risk_lower
      )
      low_hb_study_data <- study_data[low_hb_rows, ]
      output$unique_nids_lowhb <- as.integer(length(unique(low_hb_study_data$study_id)))
      output$obs_lowhb <- as.integer(nrow(low_hb_study_data))

      # Unique NID & observation count for high hb
      high_hb_rows <- which(
        study_data$ref_risk_upper <= study_data$alt_risk_lower &
          study_data$ref_risk_upper < study_data$alt_risk_upper
      )
      high_hb_study_data <- study_data[high_hb_rows, ]
      output$unique_nids_highhb <- as.integer(length(unique(high_hb_study_data$study_id)))
      output$obs_highhb <- as.integer(nrow(high_hb_study_data))

      ### Reformat "mean_rr" element for reporting (fe)
      output$lower$mean_rr_fe <- reformat_mean_rr(output$mean_rr_lower_fe)
      output$upper$mean_rr_fe <- reformat_mean_rr(output$mean_rr_upper_fe)
      output$overall$mean_rr_fe <- reformat_mean_rr(output$mean_rr_overall_fe)

      ### Remove old mean_rr fields (fe)
      output$mean_rr_lower_fe <- NULL
      output$mean_rr_upper_fe <- NULL
      output$mean_rr_overall_fe <- NULL

      # Reformat "mean_rr" element for reporting (re)
      output$lower$mean_rr_re <- reformat_mean_rr(output$lower$mean_rr)
      output$upper$mean_rr_re <- reformat_mean_rr(output$upper$mean_rr)
      output$overall$mean_rr_re <- reformat_mean_rr(output$overall$mean_rr)

      # Remove old mean_rr fields (re)
      output$lower$mean_rr <- NULL
      output$upper$mean_rr <- NULL
      output$overall$mean_rr <- NULL

      # Calculate star ratings for lower (15th-reference value) and upper (reference value-85th)
      if(split_at_value_present) {
        output$lower$star_rating <- calculate_star_rating(output$lower$score)
        output$upper$star_rating <- calculate_star_rating(output$upper$score)

      } else {
        message(
          "Skipping Star rating calculation for lower and upper: \n",
          "bophf::bop was NOT run with split_at = ",
          split_at
        )
      }
      
      # Remove lower & upper elements where unexpected (i.e., bophf::bop() run w/o split at risk value)
      if(!split_at_value_present) {
        output <- output[!(names(output) %in% c("lower", "upper"))]
        
      } else {
        message(
          "Leaving lower & upper elements intact: \n",
          "bophf::bop was run with split_at = ",
          split_at
        )
      }

      # Rename stand-alone "mean_bprf" element to "mean_bprf_risk_bounds" to avoid confusion with mean_bprf calculations under various exposure ranges
      names(output)[names(output) == "mean_bprf"] <- "mean_bprf_risk_bounds"

      # Remove "split_at" with non-numeric value (i.e., 'split_at = "tmrel"')
      output <- output[!(names(output) == "split_at" & !sapply(output, is.numeric))]
      
      # Keep only the first "split_at" element if there are multiple
      if (sum(names(output) == "split_at") > 1) {
        first_split_index <- which(names(output) == "split_at")[1]
        output <- output[setdiff(seq_along(output), which(names(output) == "split_at")[-1])]
      }

      # Remove redundant mean_bprf metrics
      output <- output[!(names(output) %in% c("mean_bprf_overall", "mean_bprf_lower", "mean_bprf_upper"))]

      # Remove fields unnecessary for reporting to minimize file size
      output <- output[!(names(output) %in% c("data", "outer_quantiles", "inner_quantiles"))]

      # Write out .yaml file
      if (split_at == 120) {
        output_file <- glue::glue("{dir}/summary_updated_global_tmrel_120.yaml")
      } else {
        output_file <- file.path(dir, "summary_updated.yaml")
      }
      cat("Output file path:", output_file, "\n")
      yaml::write_yaml(output, output_file)

      cat("Output saved for directory:", dir, "\n")
    }, error = function(e) {
      cat("Error occurred in directory:", dir, "- skipping this directory.\n")
      cat("Error message:", e$message, "\n")
      traceback()
    })
  }
}

#' Retrieve the exposure value closest to the desired exposure value ----
#'
#' @param x A data object.
#' @param target The target value to which you want to find the closest value in
#'   data object `x`.
#'
#' @return The value in `x` that is closest to `target`.
#'
#' @export
closest_val <- function(x, target) {
  x[which.min(abs(x - target))]
}

#' Calculate and format the estimates & UI ----
#'
#' @param output A list containing quantiles data frames.
#' @param quantiles_name The name of the quantile within `output` to use.
#' @param lower_risk_value Optional lower bound to filter the quantiles data.
#' @param upper_risk_value Optional upper bound to filter the quantiles data.
#'
#' @return A formatted string representing the calculated mean, lower, and upper
#'   relative risks.
#'
#' @export
calculate_rr_output <- function(output, quantiles_name, lower_risk_value = NULL, upper_risk_value = NULL) {
  assertthat::assert_that(quantiles_name %in% names(output))
  quantiles_data <- output[[quantiles_name]]

  # Handle risk_value logic
  if (!is.null(lower_risk_value) && !is.null(upper_risk_value)) {
    closest_lower_risk_value <- closest_val(quantiles_data$risk, lower_risk_value)
    closest_upper_risk_value <- closest_val(quantiles_data$risk, upper_risk_value)
    quantiles_data <- quantiles_data |> dplyr::filter(risk >= closest_lower_risk_value & risk <= closest_upper_risk_value)
  } else if (!is.null(lower_risk_value)) {
    closest_lower_risk_value <- closest_val(quantiles_data$risk, lower_risk_value)
    quantiles_data <- quantiles_data |> dplyr::filter(risk >= closest_lower_risk_value)
  } else if (!is.null(upper_risk_value)) {
    closest_upper_risk_value <- closest_val(quantiles_data$risk, upper_risk_value)
    quantiles_data <- quantiles_data |> dplyr::filter(risk <= closest_upper_risk_value)
  }

  # Calculate mean, lower, and upper values, then exponentiate and round them
  mean_val <- purrr::chuck(quantiles_data, 'mean') |> mean() |> exp() |> round(2)
  lower_val <- purrr::chuck(quantiles_data, 'lower') |> mean() |> exp() |> round(2)
  upper_val <- purrr::chuck(quantiles_data, 'upper') |> mean() |> exp() |> round(2)

  # Format values to 2 decimal places
  mean_val <- sprintf("%.2f", mean_val)
  lower_val <- sprintf("%.2f", lower_val)
  upper_val <- sprintf("%.2f", upper_val)

  # Combine the values into a single string
  result <- paste0(mean_val, " (", lower_val, " - ", upper_val, ")")

  # Return the result with a label
  return(result)
}

#' Reformat mean_rr ----
#'
#' @param mean_rr A data frame containing columns `mean`, `lower`, and `upper`
#'   for relative risk.
#'
#' @return A formatted string representing the mean, lower, and upper relative
#'   risks, or a message indicating that no data is available.
#'
#' @export
reformat_mean_rr <- function(mean_rr) {
  if (!is.null(mean_rr) && all(c("mean", "lower", "upper") %in% colnames(mean_rr))) {
    mean_val <- sprintf("%.2f", mean_rr$mean)
    lower_val <- sprintf("%.2f", mean_rr$lower)
    upper_val <- sprintf("%.2f", mean_rr$upper)
    return(paste0(mean_val, " (", lower_val, " - ", upper_val, ")"))
  } else {
    return("No data available")
  }
}

#' Derive knot bounds based on the highest weight in the model ----
#'
#' @param dir A character string specifying the directory path where the pickle
#'   files are located.
#'
#' @return A list containing the bounds of the best spline knots.
#'
#' @export
get_best_spline_knots <- function(dir) {
  reticulate::use_condaenv("PATH TO ENVIRONMENT")
  signal_model <- reticulate::py_load_object(
    list.files(dir, pattern = "^signal.*\\.pkl$", full.names = TRUE)
  )
  best_knot <- which.max(signal_model$weights)
  knot_bounds <- as.list(signal_model$ensemble_knots[best_knot, ])
  return(knot_bounds)
}

#' Calculate star_rating based on score ----
#'
#' @param score A numeric value representing the score for which the star rating
#'   is calculated.
#'
#' @return An integer representing the star rating. Must be between 0 and 5.
#'
#' @export
calculate_star_rating <- function(score) {
  if (is.nan(score)) {
    return(as.integer(0))
  } else if (score < 0.000) {
    return(as.integer(1))
  } else if (score <= 0.1398) {
    return(as.integer(2))
  } else if (score <= 0.4055) {
    return(as.integer(3))
  } else if (score <= 0.6152) {
    return(as.integer(4))
  } else if (score > 0.6152) {
    return(as.integer(5))
  } else {
    stop("Unexpected score value provided, cannot calculate star rating")
  }
}

#' Calculate RMSE ----
#'
#' @param quantiles_name Character string specifying the name of the quantiles
#'   within the `bop` list.
#' @param study_data A data frame containing observed study data, including
#'   lower and upper risk bounds.
#' @param bop A list containing model quantiles, among other elements.
#'
#' @return A numeric value representing the RMSE.
#'
#' @export
get_rmse <- function(quantiles_name, study_data, bop) {
  quantiles_data <- bop[[quantiles_name]]
  sse_val <- 0

  # Calculate midpoint for the study_data
  study_data$alt_mid_point <- 0.5 * (study_data$alt_risk_lower + study_data$alt_risk_upper)

  # Loop through each row in study_data
  for (r in seq_len(nrow(study_data))) {
    temp_exp <- study_data$alt_mid_point[r]
    temp_study_rr <- study_data$ln_rr[r]

    # Find the closest value in quantiles_data$risk to temp_exp
    selector <- quantiles_data$risk == closest_val(quantiles_data$risk, temp_exp)
    est_rr <- quantiles_data[selector, "mean"]

    # Calculate the squared difference
    point_diff <- (temp_study_rr - est_rr)^2
    sse_val <- sse_val + point_diff
  }

  # Calculate RMSE
  rmse <- sqrt(sse_val / nrow(study_data))

  return(unlist(rmse))
}

#' Calculate mean estimates for given quantiles and risk range ----
#'
#' @param quantile_dat A data table containing risk values and their
#'   corresponding quantiles.
#' @param risk_range Numeric vector of length two specifying the lower and upper
#'   bounds of risk.
#'
#' @note Used for overall, lower & upper - FE
#' @return A summarized data table with mean estimates for overall, lower, and
#'   upper quantiles.
#'
#' @export
mean_rr_for <- function(quantile_dat, risk_range) {
  print(risk_range)
  checkmate::assert_data_table(quantile_dat)
  checkmate::assert_names(names(quantile_dat), must.include = "risk")
  checkmate::assert_numeric(risk_range,
                            lower = min(quantile_dat$risk),
                            upper = max(quantile_dat$risk), len = 2,
                            any.missing = FALSE,
                            sorted = TRUE)
  quantile_dat |>
    dplyr::summarize(
      dplyr::across(
        tidyselect::all_of(c("mean", "lower", "upper")),
        \(val) mean_val_for(.data$risk, val, risk_range) |> exp()
      )
    )
}

mean_val_for <- function(risk, val, risk_range) {
  mean(
    stats::approx(
      x = risk,
      y = val,
      xout = seq(risk_range[1], risk_range[2], length.out = 1000)
    )$y
  )
}
