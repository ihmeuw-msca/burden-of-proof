#####################
##
## BoP General Edits Script (w/ options for further formatting based on R-O pair)
##
## NCH Team
##
#####################

library(dplyr)
library(tidyr)
library(data.table)
set.seed(123)

# source all necessary IHME and BoP helper functions ----------------------

withr::local_libpaths("PATH TO PACKAGES", action = "prefix")
source("PATH TO format_bop_bundle.R")

# set global variables ----------------------------------------------------
pipeline_type <- "continuous"
cause_id <- 567
acause <- "mental_unipolar"
rei_id <- 376
risk_name <- "low_hemoglobin"
bundle_id <- 999
bundle_version_id <- 999
rr_or_name_vec <- c("Relative risk (RR)", "Odds ratio (OR)")

# load in the data from a bundle ------------------------------------------
df <- fread("PATH TO MAIN DATASET") # exclcrx
# df <- fread("PATH TO POSTATAL DEPRESSION DATASET") # exclcrx_postnataldep

df$seq <- 1:nrow(df)

df <- format_bop_bundle(
  input_df = df,
  release_id = 16
)

df$acause <- acause
df$cause_id <- cause_id

df$risk_name <- risk_name
df$rei_id <- rei_id

df$bundle_id <- bundle_id
df$bundle_version_id <- bundle_version_id

df$cv_pregnant <- 1

df <- df[cause_id == cause_id & rei_id == rei_id & Parameter.Measure %in% rr_or_name_vec]

# impute missing reference bounds -----------------------------------------
df_list <- bophf::impute_missing_hb_main(
  input_df = df,
  biomarker_column_name = "biomarker_type_1",
  biomarker_unit_column_name = "unit_of_measure_1",
  ref_risk_upper_column_name = "ref_upper_limit_b1",
  ref_risk_lower_column_name = "ref_lower_limit_b1",
  alt_risk_upper_column_name = "alt_upper_limit_1",
  alt_risk_lower_column_name = "alt_lower_limit_1",
  apply_elevation_adjustment = TRUE,
  conf_elevation_adjust_column_name = "covariate_altitude_adjusted"
)

df <- copy(df_list$cleaned_df)
no_cutoff_df <- copy(df_list$no_cutoff_df)
high_hb_df <- copy(df_list$high_hb)
bad_unit_df <- copy(df_list$bad_conversion_df)
bad_bounds_df <- copy(df_list$bad_bounds_df)

# set the names of the mean and standard error columns --------------------
df_list <- bophf::transform_mean_se_bop(
  input_df = df,
  mean_column_name = "mean",
  upper_ui_column_name = "upper",
  lower_ui_column_name = "lower",
  se_column_name = "Standard.deviation"
)

# remove invalid rows -----------------------------------------------------
bad_transform_df <- copy(df_list$na_df)
df <- copy(df_list$df)

# format covariate columns ------------------------------------------------
df <- bophf::check_for_linear_dependent_covs(df)
df$cov_con_demo <- NULL # for exccrxl analysis

# df <- df %>% # for prenataldep analysis (3 observations only; need to drop all bias covs for BoP to run)
#   select(-starts_with("cov_"))
# df$cov_con_demo <- NULL # for postnataldep analysis
# df$cov_out_method <- NULL # for postnataldep analysis

# finalize BoP data set ---------------------------------------------------
df$out_lower <- NA
df$out_upper <- NA

df <- df |>
  dplyr::group_by(nid) |>
  dplyr::mutate(
    inflate_se = dplyr::case_when(
      dplyr::n_distinct(exp_timepoint) > 1 ~ dplyr::n_distinct(exp_timepoint),
      TRUE ~ 1
    ),
    ln_rr_se = ln_rr_se * sqrt(inflate_se)
  ) |>
  dplyr::ungroup()

bop_df <- bophf::get_final_bop_df(
  input_df = df,
  pipeline_type = pipeline_type,
  other_cols_to_add = c("exp_timepoint","hb_risk_type")
)

ro_pair <- paste(risk_name, acause, sep = "-")
write.csv(bop_df,paste0("bop_model_sandbox/data/",ro_pair,".csv"),row.names = F)

# plot diagnostics before running BoP -------------------------------------
source("diagnostic_code/mrbrt_analysis.R")

meta_analysis_result <- run_mrbrt(input_df = bop_df)

mr_brt_df <- copy(meta_analysis_result$mr_brt_df)

bophf::bop_forest_plot_w_mrbrt(
  input_df = mr_brt_df,
  plot_linear_space = T,
  num_decimal_points = 2,
  main_title = "MR-BRT Meta Analysis - Low Hgb as a Risk Factor for Perinatal Depression",
  output_directory = "OUTPUT PATH"
)

bophf::plot_rr_vs_alt(df = bop_df)

knot_df <- bophf::plot_spline_knots(
  df = bop_df
)

knot_df_frequency <- bophf::plot_spline_knots(
  df = bop_df
)

plotly::ggplotly(knot_df_frequency$knot_plot)

domain_spline_knots <- bophf::get_spline_knots(
  min_exposure = min(bop_df$alt_risk_lower),
  max_exposure = max(bop_df$ref_risk_upper),
  knot_value_vec = c(70,100,110)
)

knot_df_domain <- bophf::plot_spline_knots(
  df = bop_df,
  knot_type = "domain",
  knot_vec = domain_spline_knots
)

plotly::ggplotly(knot_df_domain$knot_plot)

# save copies of output files post model run ------------------------------
notes <- "Frequency knots, left linear = T, right linear = T, 10% outliered, inflated SE in Maeda article for multiple timepoints- sqrt(3) x ln_rr_se."

bophf::save_bop_outputs(
  output_directory = "OUTPUT PATH",
  risk_outcome_pair = ro_pair,
  notes = notes
)

# compare different model outputs -----------------------------------------
compare_df <- bophf::compare_summary_outputs(
  bop_snapshot_directory = file.path("OUTPUT PATH"),
  output_directory = file.path("OUTPUT PATH")
)

# Vetting 1
library(ggplot2)

df$alt_risk_midpoint <- (df$alt_risk_lower + df$alt_risk_upper) / 2
df$ref_risk_midpoint <- (df$ref_risk_lower + df$ref_risk_upper) / 2

df$hb_risk_type <- factor(df$hb_risk_type, levels = c("low_hb", "high_hb"))

df_long <- df %>%
  pivot_longer(
    cols = c(alt_risk_midpoint, ref_risk_midpoint),
    names_to = "risk_type",
    values_to = "midpoint"
  )

ggplot(df_long, aes(x = midpoint, fill = risk_type)) +
  geom_histogram(
    position = "stack", 
    bins = 20, 
    alpha = 0.7
  ) +
  facet_wrap(~hb_risk_type) +
  labs(
    title = "Distribution of Midpoints by Risk Type",
    x = "Hb Midpoint",
    y = "Count",
    fill = "Risk Type"
  ) +
  theme_minimal()

# Vetting 2
meta_analysis_result <- run_mrbrt(input_df = bop_df[hb_risk_type == "low_hb"])

mr_brt_df <- copy(meta_analysis_result$mr_brt_df)

bophf::bop_forest_plot_w_mrbrt(
  input_df = mr_brt_df,
  plot_linear_space = T,
  num_decimal_points = 2,
  main_title = "MR-BRT Meta Analysis - Low Hgb as a Risk Factor for Perinatal Depression (Low-Hb)"
)

meta_analysis_result <- run_mrbrt(input_df = bop_df[hb_risk_type == "high_hb"])

mr_brt_df <- copy(meta_analysis_result$mr_brt_df)

bophf::bop_forest_plot_w_mrbrt(
  input_df = mr_brt_df,
  plot_linear_space = T,
  num_decimal_points = 2,
  main_title = "MR-BRT Meta Analysis - Low Hgb as a Risk Factor for Perinatal Depression (High-Hb)"
)