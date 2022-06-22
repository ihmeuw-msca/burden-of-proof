# Configuration of pipeline

# User settings
# ------------------------------------------------------------------------------
USER <- Sys.getenv("USER")
WORK_DIR <- "/ihme/homes/jiaweihe/msca/mrbrt/evidence_score_pipeline"
CODE_PATH <- paste0(WORK_DIR, "/src/")
ARCHIVE <- "/mnt/team/msca/pub/archive/evidence-score/gbd2020"

# Cluster settings
# ------------------------------------------------------------------------------
PROJ <- "proj_mscm"
SINGULARITY_IMG <- "/ihme/singularity-images/rstudio/ihme_rstudio_3631.img"

# Version settings
# ------------------------------------------------------------------------------
VERSION_ID <- "prod"

# Directory settings
# ------------------------------------------------------------------------------
OUT_DIR <- paste0("/ihme/scratch/users/", USER, "/evidence_score_pipeline/", VERSION_ID, "/")
INPUT_DATA_DIR = "/home/j/temp/hkl1/mr_brt/03_evidence_score/input_data/for_ryan/version15"

# Output directory for each stage
SUB_DIRS <- c(
  paste0(OUT_DIR, "00_prepped_data"),
  paste0(OUT_DIR, "01_template_pkl_files"),
  paste0(OUT_DIR, "01_template_models"),
  paste0(OUT_DIR, "02_loglinear_models"),
  paste0(OUT_DIR, "02_loglinear_pkl_files"),
  paste0(OUT_DIR, "03_covariate_selection_models"),
  paste0(OUT_DIR, "03_covariate_selection_pkl_files"),
  paste0(OUT_DIR, "04_mixed_effects_models"),
  paste0(OUT_DIR, "04_mixed_effects_pkl_files"),
  paste0(OUT_DIR, "05_evidence_score")
) 

# data settings
# ------------------------------------------------------------------------------
ALL_RO_PAIRS <- gsub(".csv", "", list.files(INPUT_DATA_DIR))
EXCLUDED_RO_PAIRS <- c("dairy_stroke", "fruit_oral", "fruit_larynx")
RO_PAIRS <- ALL_RO_PAIRS[!(ALL_RO_PAIRS %in% EXCLUDED_RO_PAIRS)]
RO_PAIRS <- c("alcohol")
RO_PAIRS <- c("redmeat_colorectal")
RO_PAIRS <- c("lpa_ihd", "bmi_diabetes", "bmi_uter_canc", 
              "bmi_leuk", "fruit_ihd", "redmeat_colorectal", 
              "nuts_ihd", "wholegrain_ihd", "fiber_stroke", 
              "alcohol_lri", "alcohol_tb", "lung_cancer", 
              "copd", "ihd_19", "diabetes", "peptic_ulcer")

OBS_VAR <- "ln_effect"
OBS_SE_VAR <- "ln_se"
STUDY_ID_VAR <- "nid"

ALT_EXPOSURE_COLS <- c("b_0", "b_1")
REF_EXPOSURE_COLS <- c("a_0", "a_1")

# Sarah's
# RO_PAIRS <- c("air_pmhap_neo_lung", "air_pmhap_lri", "air_pmhap_t2_dm", 
#               "air_pmhap_resp_copd", "air_pmhap_cvd_stroke_60")
# ALT_EXPOSURE_COLS <- c("conc")
# REF_EXPOSURE_COLS <- c("conc_den")


# model settings
# ------------------------------------------------------------------------------
BIAS_COVARIATES_AS_INTX <- TRUE

# For diet
# DIRECTION = list(
#   calcium = "decreasing",
#   cheese = "decreasing",
#   dairy = "decreasing",
#   fiber = "decreasing",
#   fish = "decreasing",
#   fruit = "decreasing",
#   legumes = "decreasing",
#   milk = "decreasing",
#   nuts = "decreasing",
#   omega3 = "decreasing",
#   veg = "decreasing",
#   wholegrain = "decreasing",
#   pufa = "decreasing",
#   yogurt = "decreasing",
#   procmeat = "increasing",
#   redmeat = "increasing",
#   sodium = "increasing",
#   ssb = "increasing",
#   sugar = "increasing",
#   transfat = "increasing"
# )

# Get monotonicity direction.
# tmp <- read.csv("/ihme/homes/jiaweihe/msca/mrbrt/evidence_score_pipeline/all_pairs.csv")
# tmp$mono <- ifelse(tmp$type=='protective', 'decreasing', 'increasing')
# DIRECTION <- setNames(as.character(tmp$mono), tmp$ro_pair)
DIRECTION = list(
  lpa_ihd = "increasing"
)


BETA_PRIOR_MULTIPLIER = 0.1
COV_FINDER_CONFIG = list(
    pre_selected_covs = list("exposure_linear"), 
    num_samples = 1000L,
    power_range = list(-4, 4), 
    power_step_size = 0.05,
    laplace_threshold = 1e-5,
    inlier_pct = 1.0,
    bias_zero = TRUE
)

# Not used by new pipeline

# N_I_KNOTS <- 3
# PRIOR_VAR_RSLOPE = 1e-6
# PRIOR_VAR_MAXDER <- 1e-4
# MONOSPLINE_SLOPE_MULTIPLIER <- 2

# MONOSPLINE_CONFIG = list(
#     use_re = TRUE,
#     use_spline = TRUE,
#     spline_degree = 3L,
#     spline_knots_type = 'domain',
#     spline_r_linear = TRUE,
#     prior_spline_funval_uniform = array(c(-1 + 1e-6, 19)),
#     prior_spline_num_constraint_points = 150L,
#     spline_knots = array(seq(0, 1, length.out = N_I_KNOTS + 2)),
#     prior_spline_maxder_gaussian = cbind(rbind(rep(0, N_I_KNOTS), 
#       rep(sqrt(PRIOR_VAR_MAXDER), N_I_KNOTS)), c(0, sqrt(PRIOR_VAR_RSLOPE))),
#     prior_spline_der2val_gaussian = NULL,
#     prior_spline_der2val_gaussian_domain = array(c(0.0, 1.0)),
#     name = "exposure"
# )

# MONOSPLINE_BIAS_CONFIG = list(
#   spline_degree = 3L
# )

# LOGLINEAR_BIAS_CONFIG = list(
#   spline_degree = 3L
# )
