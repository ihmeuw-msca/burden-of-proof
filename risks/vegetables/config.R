# Configuration of pipeline

# Directory settings
# ------------------------------------------------------------------------------
OUT_DIR <- FILEPATH
INPUT_DATA_DIR <- FILEPATH

# Output directory for each stage
SUB_DIRS <- c(
  paste0(OUT_DIR, "00_prepped_data"),
  paste0(OUT_DIR, "01_template_pkl_files"),
  paste0(OUT_DIR, "01_template_models"),
  paste0(OUT_DIR, "02_loglinear_models"),
  paste0(OUT_DIR, "02_loglinear_pkl_files"),
  paste0(OUT_DIR, "03_covariate_selection_models"),
  paste0(OUT_DIR, "04_mixed_effects_pkl_files"),
  paste0(OUT_DIR, "05_evidence_score"),
  paste0(OUT_DIR, "05_all_plots"),
  paste0(OUT_DIR, "05_all_csvs"),
  paste0(OUT_DIR, "05_pub_bias"),
  paste0(OUT_DIR, "05_draw_csvs")
) 

# data settings
# ------------------------------------------------------------------------------
ALL_RO_PAIRS <- gsub(".csv", "", list.files(INPUT_DATA_DIR))
EXCLUDED_RO_PAIRS <- c("sugar_cvd", "sugar_obesity", "fruit_oral", "fruit_larynx", 
                       ALL_RO_PAIRS[grepl("original", ALL_RO_PAIRS)], 
                       ALL_RO_PAIRS[grepl("_stroke", ALL_RO_PAIRS)], 
                       ALL_RO_PAIRS[grepl("sugar", ALL_RO_PAIRS)])
RO_PAIRS <- ALL_RO_PAIRS[!(ALL_RO_PAIRS %in% EXCLUDED_RO_PAIRS)]
RO_PAIRS <- ALL_RO_PAIRS[grepl("veg", ALL_RO_PAIRS)] # Alternative Option 1: select a subset of RO pairs to run


OBS_VAR <- "ln_effect"
OBS_SE_VAR <- "ln_se"
STUDY_ID_VAR <- "nid"

ALT_EXPOSURE_COLS <- c("b_0", "b_1")
REF_EXPOSURE_COLS <- c("a_0", "a_1")

USE_GLOBAL_DIST_PREDICT <- F # use the data to predict = F; use the exposure model to predict = T


# model settings
# ------------------------------------------------------------------------------
BIAS_COVARIATES_AS_INTX <- TRUE

# For diet
DIRECTION = list(
    veg = "decreasing",
  )

BETA_PRIOR_MULTIPLIER = 0.1 #used in covfinder and final model on covs


colnames(PRE_SELECTED_COVS) <- c("ro_pair", "cov")

COV_FINDER_CONFIG = list(
  #pre_selected_covs = list("signal"), 
  num_samples = 1000L,
  power_range = list(-4, 4), 
  power_step_size = 0.05,
  laplace_threshold = 1e-5,
  inlier_pct = 1, #since we trim in stage 1
  bias_zero = TRUE
)

INLIER_PCT <- 0.9 # 0.9 standard trimming


N_I_KNOTS <- 2L
PRIOR_VAR_RSLOPE = 1e-6 #originally 1e-6
PRIOR_VAR_MAXDER <- 1e-4

# Monotonic risks will have monotonicity constraint included
CONFIG = list(
  use_spline = TRUE,
  use_re = FALSE,
  spline_degree = 2L, 
  spline_knots_type = 'domain',
  spline_r_linear = TRUE,
  spline_l_linear = FALSE,
  prior_spline_funval_uniform = array(c(-1 + 1e-6, 19)),     
  prior_spline_num_constraint_points = 150L,
  spline_knots = array(seq(0, 1, length.out = N_I_KNOTS + 2)),
  prior_spline_maxder_gaussian = cbind(rbind(rep(0, N_I_KNOTS),
                                             rep(Inf, N_I_KNOTS)),
                                       c(0, sqrt(PRIOR_VAR_RSLOPE)))
)


J_N_I_KNOTS <- 3L
# 
J_SHAPED_CONFIG = list(
  use_spline = TRUE,
  use_re = FALSE,
  spline_degree = 2L, 
  spline_knots_type = 'domain',
  spline_r_linear = TRUE,
  spline_l_linear = TRUE,
  prior_spline_funval_uniform = array(c(-1 + 1e-6, 19)),     
  prior_spline_num_constraint_points = 150L,
  spline_knots = array(seq(0, 1, length.out = J_N_I_KNOTS + 2)),
  prior_spline_maxder_gaussian = cbind(c(0, sqrt(PRIOR_VAR_RSLOPE)), 
                                       rbind(rep(0, J_N_I_KNOTS-1),
                                             rep(Inf, J_N_I_KNOTS-1)),
                                       c(0, sqrt(PRIOR_VAR_RSLOPE)))
)

