# -------------------------------
# Load Required Libraries
# -------------------------------
library(data.table)
library(yaml)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)

# -------------------------------
# Load Outcome Metadata
# -------------------------------

# List of potential file paths
paths <- c(
  ## Overall
  "PATH TO SUMMARY TABLE",
  "PATH TO SUMMARY TABLE",
  
  ## Trimester 1
  "PATH TO SUMMARY TABLE",
  "PATH TO SUMMARY TABLE",
  
  ## Trimester 2
  "PATH TO SUMMARY TABLE",
  "PATH TO SUMMARY TABLE",
  
  ## Trimester 3
  "PATH TO SUMMARY TABLE",
  "PATH TO SUMMARY TABLE",
  
  ## LBW PTB severities
  "PATH TO SUMMARY TABLE",
  
  ## cov_socio_adj
  "PATH TO SUMMARY TABLE",
  
  ## incl_reverse_causal
  "PATH TO SUMMARY TABLE",
  
  ## incl_undefined_outcomes
  "PATH TO SUMMARY TABLE"
)

# Filter to only paths where the file exists
valid_paths <- paths[file.exists(paths)]

# Read and combine all valid CSVs into one data.table
if (length(valid_paths) > 0) {
  outcome_map <- rbindlist(
    lapply(valid_paths, fread, stringsAsFactors = FALSE),
    use.names = TRUE,
    fill = TRUE
  )
} else {
  stop("None of the specified paths contain the file.")
}


# Create parent_dir by removing the filename from the filepath
outcome_map[, parent_dir := gsub("/summary_updated.yaml$", "", filepath)]

# -------------------------------
# Helper: Extract ggplot legend
# -------------------------------
g_legend <- function(a.gplot) {
  tmp <- ggplotGrob(a.gplot)
  gtable::gtable_filter(tmp, "guide-box")
}

# -------------------------------
# Helper: Functional TMREL Function
# -------------------------------
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

# -------------------------------
# TMREL Evaluation Loop
# -------------------------------
tmrel_results_list <- list()

for (i in 1:nrow(outcome_map)) {
  parent_dir <- outcome_map$parent_dir[i]
  f <- outcome_map$outcome[i]
  
  # Define file paths
  yaml_path <- file.path(parent_dir, "summary_updated.yaml")
  csv_path  <- file.path(parent_dir, "inner_quantiles.csv")
  
  # Skip if missing required files
  if (!file.exists(yaml_path) || !file.exists(csv_path)) {
    warning(paste("Missing file(s) for:", f, "in", parent_dir))
    next
  }
  
  message(paste("Processing:", f))
  
  # Load data
  a_summary <- fread(csv_path)
  a_yaml <- yaml::read_yaml(yaml_path)
  
  # Normalize risk estimates to the minimum risk point
  tmrel_index <- which.min(a_summary$'0.5')
  a_summary$normalized_mean  <- a_summary$'0.5'   - a_summary$'0.5'[tmrel_index]
  a_summary$normalized_lower <- a_summary$'0.975' - a_summary$'0.975'[tmrel_index]
  a_summary$normalized_upper <- a_summary$'0.025' - a_summary$'0.025'[tmrel_index]
  
  # Fit smoothing spline
  slope_func <- smooth.spline(x = a_summary$risk, y = a_summary$normalized_mean)
  
  # Derivatives
  a_summary$slope     <- predict(slope_func, deriv = 1)$y
  a_summary$inflection <- predict(slope_func, deriv = 2)$y
  
  # Define range for functional TMREL search
  lower_hb_bound <- 110
  upper_hb_bound <- 140
  to_check_index_vec <- which(a_summary$risk >= lower_hb_bound & a_summary$risk <= upper_hb_bound)
  non_anemic_index <- min(to_check_index_vec)
  
  # Flat slope threshold
  flat_slope <- as.numeric(quantile(abs(a_summary$slope[to_check_index_vec]), pnorm(-2)))
  
  # Retrieve BOP TMREL and compute functional
  bop_tmrel <- a_yaml$tmrel
  tmrel_result <- find_tmrel(a_summary, to_check_index_vec, flat_slope, bop_tmrel, non_anemic_index)
  
  # Save TMRELs
  tmrel_results_list[[f]] <- list(
    functional_tmrel = tmrel_result,
    bop_tmrel = bop_tmrel
  )
  yaml::write_yaml(tmrel_results_list[[f]], file.path(parent_dir, "tmrel_result.yaml"))
  
  # -------------------------------
  # Plot Derivatives + Risk Curve
  # -------------------------------
  vline_data <- data.frame(
    label = c("BOP TMREL", "Functional TMREL"),
    xintercept = c(bop_tmrel, tmrel_result),
    col = c("red", "blue")
  )
  
  common_theme <- theme_minimal() +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "grey70", size = 0.4),
      axis.line.y = element_line(color = "grey70", size = 0.4),
      axis.line.x = element_line(color = "grey70", size = 0.4)
    )
  
  # First Derivative
  p0 <- ggplot(a_summary, aes(x = risk, y = slope)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "black", size = 1) +  # Thick y=0 line
    geom_vline(data = vline_data, aes(xintercept = xintercept, color = label), linetype = "dashed") +
    scale_color_manual(name = "Legend", values = c("BOP TMREL" = "red", "Functional TMREL" = "blue")) +
    labs(title = paste("First Derivative -", f), x = "Hb (g/L)", y = "ln (RR)") +
    common_theme +
    theme(legend.position = "bottom")
  
  # Second Derivative
  p1 <- ggplot(a_summary, aes(x = risk, y = inflection)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "black", size = 1) +  # Thick y=0 line
    geom_vline(data = vline_data, aes(xintercept = xintercept, color = label), linetype = "dashed") +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = paste("Second Derivative -", f), x = "Hb (g/L)", y = "ln (RR)") +
    common_theme
  
  # Mean Risk Plot
  p2 <- ggplot(a_summary, aes(x = risk, y = normalized_mean)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "black", size = 1) +  # Thick y=0 line
    geom_vline(data = vline_data, aes(xintercept = xintercept, color = label), linetype = "dashed") +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = paste("Mean Risk -", f), x = "Hb (g/L)", y = "ln (RR)") +
    common_theme
  
  # Save to PDF
  legend <- g_legend(p0)
  pdf(file.path(parent_dir, "derivative_plot.pdf"), width = 15, height = 7.5)
  
  grid.arrange(
    arrangeGrob(p0 + theme(legend.position = "none"), p1, p2, ncol = 3),
    legend,
    nrow = 2,
    heights = c(9, 1.5)
  )
  dev.off()
}

