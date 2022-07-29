
prep_diet_data <- function(
    ro_pair, obs_var, obs_se_var, ref_vars, alt_vars, allow_ref_gt_alt = FALSE,
    study_id_var = NA,
    drop_x_covs = NA, keep_x_covs = NA, drop_z_covs = NA, keep_z_covs = NA,
    diet_dir = NA,
    verbose = TRUE) {
  
  require(dplyr)
  require(rlang)
  
  if (verbose) cat(ro_pair, "\n")
  if (!is.na(drop_x_covs) & !is.na(keep_x_covs)) stop("Cannot specify both drop and keep X-covs")
  if (!is.na(drop_z_covs) & !is.na(keep_z_covs)) stop("Cannot specify both drop and keep Z-covs")
  
  df <- read.csv(paste0(diet_dir, "/", ro_pair, ".csv")) %>%
    filter(complete.cases(.[, c(ref_vars, alt_vars)]))
  
  if (nrow(df) == 0) stop("No observations with non-missing exposure columns")
  
  # # convert non-binary covariates into dummy variables
  # # and create list of bias covariates for the analysis
  create_dummy_vars <- function(dat, varname, reference_level) {
    dev <- FALSE
    if (dev) {
      dat <- data.frame(x1 = sample(c("a", "b", "c"), 30, TRUE))
      varname <- "x1"
      reference_level <- "a"
    }
    vec <- as.data.frame(dat)[, varname]
    lvls <- unique(vec)[!unique(vec) == reference_level]
    dat2 <- as.data.frame(do.call("cbind", lapply(lvls, function(x) as.integer(vec == x))))
    if(!is_empty(dat2)){
      names(dat2) <- paste0(varname, "_", lvls)
    }
    return(dat2)
  }
  
  confounders <- names(df)[grepl('confounders_', names(df))]
  cvs <- names(df)[grepl('cv_', names(df))]
  
  data_cols <- c(cvs)
  
  bias_covs <- c()
  for (cov in data_cols[data_cols %in% names(df)]) {
    
    dev <- FALSE
    if (dev) {
      cov <- "follow_up"
    }
    
    if (any(is.na(df[, cov]))) next 
    
    if (all(df[, cov] == round(df[, cov]))) {
      df[, cov] <- as.integer(df[, cov])
    } else {
      stop(paste0("Bias covariate '", cov, "' is not of type integer"))
    }
    bias_covs <- c(bias_covs, cov)
  }
  
  # use SVD to prevent adding collinear variables
  bias_covs_tmp <- c()
  
  # sort the bias_covs to make sure the cv_adj is always included
  bias_covs <- sort(bias_covs)
  
  for (bias_cov in bias_covs) {
    dev <- FALSE
    if (dev) {
      bias_cov <- "exposure_3"
    }
    d <- svd(cbind(df[, bias_covs_tmp], df[, bias_cov]))$d
    if (d[length(d)] > 1e-10) bias_covs_tmp <- c(bias_covs_tmp, bias_cov)
  }
  
  bias_covs <- bias_covs_tmp
  
  # warn if cv_adj is not selected
  if(!'cv_adj' %in% bias_covs) message("Warning: cv_adj is not selected")
  
  # dataset
  # NOTE: these covs cannot have missingness!
  df <- df[, c("nid", "ln_effect", "ln_se", ALT_EXPOSURE_COLS, REF_EXPOSURE_COLS, 'percent_male', 'age_start', 'age_end', 'age_ref', bias_covs)] %>%
    filter(complete.cases(.)) %>%
    arrange(nid)
  
  ##cov inclusion/exclusion
  # -- X
  if (!is.na(keep_x_covs)) {
    if (!all(keep_x_covs %in% bias_covs)) {
      stop("One or more provided X-covs not allowed.")
    } else {
      x_covs <- keep_x_covs
    }
  } else if (!is.na(drop_x_covs)) {
    x_covs <- bias_covs[!bias_covs %in% drop_x_covs]
  } else {
    x_covs <- bias_covs
  }
  
  #-- Z
  if (!is.na(keep_z_covs)) {
    if (!all(keep_z_covs %in% bias_covs)) {
      stop("One or more provided Z-covs not allowed.")
    } else {
      z_covs <- keep_z_covs
    }
  } else if (!is.na(drop_z_covs)) {
    z_covs <- bias_covs[!bias_covs %in% drop_z_covs]
  } else {
    z_covs <- bias_covs
  }
  
  out <- list(
    df=df, ro_pair=ro_pair, x_covs=x_covs, z_covs=z_covs,
    obs_var=obs_var, obs_se_var=obs_se_var,
    ref_vars=ref_vars, alt_vars=alt_vars,
    study_id_var=study_id_var,
    allow_ref_gt_alt=allow_ref_gt_alt
  )
  return(out)
}

