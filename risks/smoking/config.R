# Configuration of pipeline

# User settings
# ------------------------------------------------------------------------------
USER <- Sys.getenv("USER")
WORK_DIR <- "[working directory]"
CODE_PATH <- paste0(WORK_DIR, "/src/")

# Cluster settings
# ------------------------------------------------------------------------------
PROJ <- "[project name]"
SINGULARITY_IMG <- "[R image directory]"

# Version settings
# ------------------------------------------------------------------------------
VERSION_ID <- "prod"

# Directory settings
# ------------------------------------------------------------------------------
OUT_DIR <- "[output directory]"
INPUT_DATA_DIR = "[input data directory]"

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

# create directories
for (direc in SUB_DIRS){
  dir.create(direc, showWarnings = F)
}

# data settings
# ------------------------------------------------------------------------------
ALL_RO_PAIRS <- gsub(".csv", "", list.files(INPUT_DATA_DIR))
EXCLUDED_RO_PAIRS <- c("dairy_stroke", "fruit_oral", "fruit_larynx")
RO_PAIRS <- ALL_RO_PAIRS[!(ALL_RO_PAIRS %in% EXCLUDED_RO_PAIRS)]

OBS_VAR <- "ln_effect"
OBS_SE_VAR <- "ln_se"
STUDY_ID_VAR <- "nid"

ALT_EXPOSURE_COLS <- c("b_0", "b_1")
REF_EXPOSURE_COLS <- c("a_0", "a_1")

# model settings
# ------------------------------------------------------------------------------
BIAS_COVARIATES_AS_INTX <- TRUE

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

PRIOR_VAR_RSLOPE = 1e-6
PRIOR_VAR_MAXDER <- 1e-4
MONOSPLINE_SLOPE_MULTIPLIER <- 2

MONOSPLINE_BIAS_CONFIG = list(
  spline_degree = 3L
)

LOGLINEAR_BIAS_CONFIG = list(
  spline_degree = 3L
)
