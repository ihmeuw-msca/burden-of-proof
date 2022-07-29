#----------------------------------------------------------------------------------
# Prep bias covariates for diet MR-BRT
#----------------------------------------------------------------------------------



prep_diet_data <- function(
  ro_pair, obs_var, obs_se_var, ref_vars, alt_vars, allow_ref_gt_alt = FALSE,
  study_id_var = NA,
  drop_x_covs = NA, keep_x_covs = NA, drop_z_covs = NA, keep_z_covs = NA,
  diet_dir = "FILEPATH",
  verbose = TRUE) {
  
  
  require(dplyr)
  
  if (verbose) cat(ro_pair, "\n")
  if (!is.na(drop_x_covs) & !is.na(keep_x_covs)) stop("Cannot specify both drop and keep X-covs")
  if (!is.na(drop_z_covs) & !is.na(keep_z_covs)) stop("Cannot specify both drop and keep Z-covs")
  
  df <- read.csv(paste0(diet_dir, "/", ro_pair, ".csv")) %>%
    filter(complete.cases(.[, c(ref_vars, alt_vars)]))
  
  if (nrow(df) == 0) stop("No observations with non-missing exposure columns")
  
  # identify and fix rows where reference (a) is higher than alternative (b)
  # -- doesn't apply to J-shaped curves where ref isn't necessarily lower than alt
  
  if (!allow_ref_gt_alt) {
    rev_idx <- apply(df[, ref_vars], 1, mean) > apply(df[, alt_vars], 1, mean)
    
    df[rev_idx, paste0(alt_vars, "_tmp")] <- df[rev_idx, ref_vars]
    df[rev_idx, ref_vars] <- df[rev_idx, alt_vars]
    df[rev_idx, alt_vars] <- df[rev_idx, paste0(alt_vars, "_tmp")]
    df[rev_idx, obs_var] <- -1 * df[rev_idx, obs_var]
  }
  
  
  # set incidence AND mortality to 0 for both of those variables
  if ("incidence" %in% names(df) & "mortality" %in% names(df)) {
    if (!any(is.na(c(df$incidence, df$mortality)))) {
      inc_mort_idx = (df$incidence + df$mortality) == 2
      df[inc_mort_idx, "incidence"] <- 0
      df[inc_mort_idx, "mortality"] <- 0
    }
  }
  
  
  # turn follow-up into dichotomous
  if ("follow_up" %in% names(df)) {
    df <- df %>%
      mutate(
       value_of_duration_fup = follow_up,
        follow_up = ifelse(follow_up < 10, 1, 0) )
  }
  
  # convert non-binary covariates into dummy variables
  # and create list of bias covariates for the analysis
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
    names(dat2) <- paste0(varname, "_", lvls)
    return(dat2)
  }
  
  data_cols = c(
    'exposure_1', 'exposure_2', 'outcome_1', 'outcome_2',
    'exposure_3', 'confounder_1', 'confounder_2',
    'selection_bias', 'reverse_causation', 
    'incidence', 'mortality','follow_up'
  )
  
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
    has_multiple_categories <- !all(unique(df[, cov]) %in% c(0, 1, NA))
    
    if (has_multiple_categories) {
      df_newvars <- create_dummy_vars(
        dat = df, 
        varname = cov, 
        reference_level = min(df[, cov], na.rm = TRUE)
      )
      df <- as.data.frame(cbind(df, df_newvars))
      bias_covs <- c(bias_covs, names(df_newvars))
    } else {
      bias_covs <- c(bias_covs, cov)
    }
  }
  
  # if no incidence + mortality data is present, don't use those variables
  if (all(c("incidence", "mortality") %in% names(df))) {
    if (max(df$incidence + df$mortality, na.rm = TRUE) - min(df$incidence + df$mortality, na.rm = TRUE) == 0) {
      bias_covs <- bias_covs[!bias_covs %in% c("incidence", "mortality")]
    }
  } else {
    # remove 'incidence' and 'mortality' if 0 or 1 of them are in the dataset
    bias_covs <- bias_covs[!bias_covs %in% c("incidence", "mortality")] 
  }
  
  
  # make sure bias covariates have at least two studies in both
  # reference and alternative (assuming study count > 4?)
  for (bias_cov in bias_covs) {
    dev <- FALSE
    if (dev) {
      bias_cov <- "exposure_3"
    }

    source_counts <- df %>%
      group_by(!! rlang::sym(bias_cov)) %>%
      summarize(n_studies = length(unique(nid)) ) %>%
      filter(n_studies >= 2) %>%
      as.data.frame(.)

    if (!all(0:1 %in% source_counts[,1]) & length(unique(df$nid)) > 4) {
      bias_covs <- bias_covs[!bias_covs == bias_cov]
    }
  }

  
  # use SVD to prevent adding collinear variables
  bias_covs_tmp <- c()
  
  for (bias_cov in bias_covs) {
    dev <- FALSE
    if (dev) {
      bias_cov <- "exposure_3"
    }
    d <- svd(cbind(df[, bias_covs_tmp], df[, bias_cov]))$d
    if (d[length(d)] > 1e-10) bias_covs_tmp <- c(bias_covs_tmp, bias_cov)
  }
  
  bias_covs <- bias_covs_tmp
  
  
  #make sure each cov has multiple values
  for (bias_cov in bias_covs){
    
    tmp_df <- df[, c(bias_cov)]
    
    if(length(unique(tmp_df)) < 2 ){
      bias_covs <- bias_covs[!bias_covs == bias_cov]
    }
  }
  
  
  # dataset
  df <- df[, c("nid", "ln_effect", "ln_se", ALT_EXPOSURE_COLS, REF_EXPOSURE_COLS, bias_covs)] %>%
    filter(complete.cases(.)) %>%
    arrange(nid)
  
  
  # cov inclusion/exclusion
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
