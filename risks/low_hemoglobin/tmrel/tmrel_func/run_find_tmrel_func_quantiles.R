# Quick check of first and second derivative curves and potential empirical TMREL

run_find_tmrel_func_quantiles <- function(outputs_dir, analysis) {
  
  box::use(
    data.table,
    yaml,
    ggplot2,
    grDevices,
    gridExtra,
    gtable,
    stats,
  )
  
  # Validate outputs_dir
  checkmate::assertTRUE(dir.exists(outputs_dir))
  checkmate::assertTRUE(
    any(stringr::str_detect(outputs_dir, c(
      "vlbw", "elbw", "lbw", "lga", "antenatal_hem", "mat_hem", "mat_mort", 
      "mat_sepsis", "neo_mort", "neo_sepsis", "ppd", "preecl", "mptb", 
      "vptb", "eptb", "ptb", "sga", "stillbirth"
    )))
  )
  
  outcome_title <- dplyr::case_when(
    stringr::str_detect(outputs_dir, "eptb") ~ "Extremely preterm birth",
    stringr::str_detect(outputs_dir, "vptb") ~ "Very preterm birth",
    stringr::str_detect(outputs_dir, "mptb") ~ "Moderate preterm birth",
    stringr::str_detect(outputs_dir, "ptb") ~ "Preterm birth",
    stringr::str_detect(outputs_dir, "sga") ~ "Small-for-gestational-age birth",
    stringr::str_detect(outputs_dir, "lga") ~ "Large-for-gestational-age birth",
    stringr::str_detect(outputs_dir, "elbw") ~ "Extremely low birth weight",
    stringr::str_detect(outputs_dir, "vlbw") ~ "Very low birth weight",
    stringr::str_detect(outputs_dir, "lbw") ~ "Low birth weight",
    stringr::str_detect(outputs_dir, "antenatal_hem") ~ "Antepartum hemorrhage",
    stringr::str_detect(outputs_dir, "mat_hem") ~ "Postpartum hemorrhage",
    stringr::str_detect(outputs_dir, "mat_mort") ~ "All-cause maternal mortality",
    stringr::str_detect(outputs_dir, "mat_sepsis") ~ "Maternal sepsis",
    stringr::str_detect(outputs_dir, "neo_mort") ~ "All-cause neonatal mortality",
    stringr::str_detect(outputs_dir, "neo_sepsis") ~ "Neonatal sepsis",
    stringr::str_detect(outputs_dir, "ppd") ~ "Peripartum depression",
    stringr::str_detect(outputs_dir, "preecl") ~ "Preeclampsia",
    stringr::str_detect(outputs_dir, "stillbirth") ~ "Stillbirth"
  )
  
  # For each trial folder, calculate TMREL values and write .yaml to folder ----
  
  if (stringr::str_detect(outputs_dir, "trial")) {
    trial_dirs <- outputs_dir
  } else {
    trial_dirs <- list.dirs(outputs_dir, full.names = TRUE, recursive = TRUE)
  }
  
  for (dir in trial_dirs) {
    
    # Define file paths
    yaml_path <- file.path(dir, "summary_updated.yaml")
    csv_path  <- file.path(dir, "inner_quantiles.csv")
    
    # Skip if missing required files
    if (!file.exists(yaml_path) || !file.exists(csv_path)) {
      warning(paste("Missing file(s) for:", dir))
      next
    }
    
    message(paste("Processing:", dir))
    
    # Load data
    a_summary <- data.table$fread(csv_path)
    a_yaml <- yaml$read_yaml(yaml_path)
    
    # Normalize risk estimates to the minimum risk point
    tmrel_index <- which.min(a_summary$'0.5')
    a_summary$normalized_mean  <- a_summary$'0.5'   - a_summary$'0.5'[tmrel_index]
    a_summary$normalized_lower <- a_summary$'0.975' - a_summary$'0.975'[tmrel_index]
    a_summary$normalized_upper <- a_summary$'0.025' - a_summary$'0.025'[tmrel_index]
    
    # Fit smoothing spline
    slope_func <- stats$smooth.spline(x = a_summary$risk, y = a_summary$normalized_mean)
    
    # Derivatives
    a_summary$slope <- stats$predict(slope_func, deriv = 1)$y
    a_summary$inflection <- stats$predict(slope_func, deriv = 2)$y
    
    # Define range for functional TMREL search
    lower_hb_bound <- 110
    upper_hb_bound <- 140
    to_check_index_vec <- which(a_summary$risk >= lower_hb_bound & a_summary$risk <= upper_hb_bound)
    non_anemic_index <- min(to_check_index_vec)
    
    # Flat slope threshold
    flat_slope <- as.numeric(stats$quantile(abs(a_summary$slope[to_check_index_vec]), stats$pnorm(-2)))
    
    # Retrieve BOP TMREL and compute functional
    bop_tmrel <- a_yaml$tmrel
    tmrel_result <- find_tmrel(a_summary, to_check_index_vec, flat_slope, bop_tmrel, non_anemic_index)
    
    # Write TMRELs to trial output folder
    tmrel_results_list <- list(
      functional_tmrel = tmrel_result,
      bop_tmrel = bop_tmrel
    )
    yaml$write_yaml(tmrel_results_list, file.path(dir, "tmrel_result.yaml"))
    
    # Plot Derivatives + Risk Curve ----
    
    vline_data <- data.frame(
      label = c("BOP TMREL", "Functional TMREL"),
      xintercept = c(bop_tmrel, tmrel_result),
      col = c("red", "blue")
    )
    
    common_theme <- ggplot2$theme_minimal() +
      ggplot2$theme(
        legend.position = "none",
        axis.line = ggplot2$element_line(color = "grey70", linewidth = 0.4),
        axis.line.y = ggplot2$element_line(color = "grey70", linewidth = 0.4),
        axis.line.x = ggplot2$element_line(color = "grey70", linewidth = 0.4),
        plot.title = ggplot2::element_text(size = 12),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
    
    if (tmrel_result == bop_tmrel) {
      ftmrel_title <- "N/A"
    } else {
      ftmrel_title <- paste(round(tmrel_result, 0), "g/L")
    }
    
    # Adjust outcome title for special analyses
    outcome_title_modified <- outcome_title
    if (analysis == "reverse_causal_agnostic_incl_severe_preecl") {
      outcome_title_modified <- paste0(outcome_title, " & severe preeclampsia")
    } else if (analysis == "reverse_causal_agnostic_incl_severe_preecl_ecl") {
      outcome_title_modified <- paste0(outcome_title, ", severe preeclampsia & eclampsia")
    }
    
    # Determine descriptive label for analysis
    analysis_label <- dplyr::case_when(
      analysis == "trim1" ~ "(trimester 1)",
      analysis == "trim2" ~ "(trimester 2)",
      analysis == "trim3" ~ "(trimester 3)",
      analysis == "reverse_causal_agnostic" ~ "(agnostic to reverse causality)",
      analysis == "reverse_causal_agnostic_incl_severe_preecl" ~ "(agnostic to reverse causality)",
      analysis == "reverse_causal_agnostic_incl_severe_preecl_ecl" ~ "(agnostic to reverse causality)",
      analysis == "socio_conf_optimal" ~ "(restricted to adjustment for socioeconomic confounders)",
      analysis == "incl_undef" ~ "(including observations with undefined outcomes)",
      analysis == "overall" ~ "",
      TRUE ~ paste0("(", analysis, ")")  # fallback
    )
    
    make_title <- function(prefix, outcome, analysis_label) {
      if (analysis_label == "") {
        paste0(prefix, ": ", outcome)
      } else {
        paste0(prefix, ": ", outcome, "\n", analysis_label)
      }
    }
    
    title_text_1st_deriv  <- make_title("First derivative",  outcome_title_modified, analysis_label)
    title_text_2nd_deriv  <- make_title("Second derivative", outcome_title_modified, analysis_label)
    title_text_mean_risk  <- make_title("Mean risk",         outcome_title_modified, analysis_label)
    
    # First Derivative
    p0 <- ggplot2$ggplot(a_summary, ggplot2$aes(x = risk, y = slope)) +
      ggplot2$geom_line() +
      ggplot2$geom_hline(yintercept = 0, color = "black", linewidth = 1) +  # Thick y=0 line
      ggplot2$geom_vline(data = vline_data, ggplot2$aes(xintercept = xintercept, color = label), linetype = "dashed") +
      ggplot2$scale_color_manual(name = "Legend", values = c("BOP TMREL" = "red", "Functional TMREL" = "blue")) +
      ggplot2$labs(
        title = title_text_1st_deriv,
        subtitle = paste0("BoP TMREL = ", round(bop_tmrel, 0), " g/L; Functional TMREL = ", ftmrel_title), 
        x = "Hemoglobin (g/L)",
        y = "ln (RR)"
      ) +
      common_theme +
      ggplot2$theme(legend.position = "bottom")
    
    # Second Derivative
    p1 <- ggplot2$ggplot(a_summary, ggplot2$aes(x = risk, y = inflection)) +
      ggplot2$geom_line() +
      ggplot2$geom_hline(yintercept = 0, color = "black", linewidth = 1) +  # Thick y=0 line
      ggplot2$geom_vline(data = vline_data, ggplot2$aes(xintercept = xintercept, color = label), linetype = "dashed") +
      ggplot2$scale_color_manual(values = c("red", "blue")) +
      ggplot2$labs(
        title = title_text_2nd_deriv,
        subtitle = paste0("BoP TMREL = ", round(bop_tmrel, 0), " g/L; Functional TMREL = ", ftmrel_title), 
        x = "Hemoglobin (g/L)",
        y = "ln (RR)"
      ) +
      common_theme
    
    # Mean Risk Plot
    p2 <- ggplot2$ggplot(a_summary, ggplot2$aes(x = risk, y = normalized_mean)) +
      ggplot2$geom_line() +
      ggplot2$geom_hline(yintercept = 0, color = "black", linewidth = 1) +  # Thick y=0 line
      ggplot2$geom_vline(data = vline_data, ggplot2$aes(xintercept = xintercept, color = label), linetype = "dashed") +
      ggplot2$scale_color_manual(values = c("red", "blue")) +
      ggplot2$labs(
        title = title_text_mean_risk,
        subtitle = paste0("BoP TMREL = ", round(bop_tmrel, 0), " g/L; Functional TMREL = ", ftmrel_title), 
        x = "Hemoglobin (g/L)",
        y = "ln (RR)"
      ) +
      common_theme
    
    # Save to PDF
    legend <- gtable$gtable_filter(ggplot2$ggplotGrob(p0), "guide-box")
    grDevices$pdf(file.path(dir, "derivative_plot.pdf"), width = 16, height = 7.5)
    
    gridExtra$grid.arrange(
      gridExtra$arrangeGrob(p0 + ggplot2$theme(legend.position = "none"), p1, p2, ncol = 3),
      legend,
      nrow = 2,
      heights = c(9, 1.5)
    )
    grDevices$dev.off()
  }
}


#' Calculate functional TMREL
find_tmrel <- function(a_summary, to_check_index_vec, flat_slope, bop_tmrel, non_anemic_index) {
  found_inflection <- FALSE
  found_good_point <- FALSE
  tmrel <- NA  
  root_second_derivative <- sign(a_summary$inflection[(non_anemic_index - 1)])
  
  for (i in to_check_index_vec) {
    current_second_derivative <- sign(a_summary$inflection[i])
    
    if (!found_inflection && root_second_derivative != current_second_derivative) {
      found_inflection <- TRUE
    }
    
    current_slope <- abs(a_summary$slope[i])
    
    if (found_inflection && current_slope < flat_slope) {
      tmrel <- a_summary$risk[i]
      found_good_point <- TRUE
      message(paste("Functional TMREL Found:", tmrel))
      break
    }
  }
  
  if (!found_good_point) {
    tmrel <- bop_tmrel
    message("No functional TMREL found, using BOP TMREL.")
  }
  
  return(tmrel)
}
