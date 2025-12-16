# source necessary libraries for mrbrt ------------------------------------

library(reticulate, lib.loc = "PATH TO SOURCE")
Sys.setenv("RETICULATE_PYTHON" = "PATH TO SOURCE")
library(mrbrt003, lib.loc = "PATH TO PACKAGE")
withr::local_libpaths("PATH TO PACKAGE", action = "prefix")
library(data.table)
library(checkmate)
set.seed(123)

run_mrbrt_publication <- function(input_df){
  
  # check for necessary columns in bop data.table ---------------------------
  
  assert_data_frame(input_df)
  assert_subset(
    x = c("ln_rr", "ln_rr_se", "study_id"),
    choices = colnames(input_df),
    empty.ok = FALSE
  )
  
  # prep mrbrt model --------------------------------------------------------
  
  dt <- copy(input_df)
  
  dat1 <- MRData()
  dat1$load_df(
    data = dt,
    col_obs = "ln_rr",
    col_obs_se = "ln_rr_se",
    col_study_id = "study_id"
  )
  
  #add intercept
  mrbrt_covs <- list()
  mrbrt_covs <- append(LinearCovModel("intercept", use_re = T), mrbrt_covs)
  
  # run mr-brt model --------------------------------------------------------
  
  fit1 <- MRBRT(
    data = dat1,
    cov_models = mrbrt_covs
  )
  
  fit1$fit_model(inner_print_level = 5L, inner_max_iter = 50000L)
  
  # predict values ----------------------------------------------------------
  
  dat_pred <- MRData()
  dat_pred$load_df(
    data = dt,
    col_study_id = "study_id"
  )
  dt$pred_ln_rr <- fit1$predict(
    data = dat_pred,
    predict_for_study = F, # this gives us overall ln(RR), not study specific
    sort_by_data_id = F
  )
  
  # get beta and gamma values -----------------------------------------------
  
  coefs_dt <- as.data.table(fit1$summary())
  
  # get UI from FE ----------------------------------------------------------
  
  n_samples1 <- 1000L
  samples1 <- fit1$sample_soln(
    sample_size = n_samples1
  )
  draws1 <- fit1$create_draws(
    data = dat_pred,
    beta_samples = samples1[[1]],
    gamma_samples = matrix(rep(0, n_samples1), ncol = 1),
    random_study = FALSE
  )
  
  dt$pred1_lo <- apply(draws1, 1, function(x) quantile(x, 0.025))
  dt$pred1_up <- apply(draws1, 1, function(x) quantile(x, 0.975))
  dt$pred1_ln_se <- (dt$pred1_up - dt$pred1_lo) / (2 * qnorm(0.975))
  
  # get UI from RE and FE ---------------------------------------------------
  
  n_samples2 <- 1000L
  # each row is the same; point estimates of gamma rather than samples
  gamma_vals2 <- matrix(
    data = fit1$gamma_soln,
    nrow = n_samples2,
    ncol = length(fit1$gamma_soln),
    byrow = TRUE # 'byrow = TRUE' is important to include
  )
  samples2 <- fit1$sample_soln(sample_size = n_samples2)
  draws2 <- fit1$create_draws(
    data = dat_pred,
    beta_samples = samples2[[1]],
    gamma_samples = gamma_vals2,
    random_study = TRUE
  )
  
  dt$pred2_lo <- apply(draws2, 1, function(x) quantile(x, 0.025))
  dt$pred2_up <- apply(draws2, 1, function(x) quantile(x, 0.975))
  dt$pred2_ln_se <- (dt$pred2_up - dt$pred2_lo) / (2 * qnorm(0.975))
  
  # calculate linear space rr and se ----------------------------------------
  
  dt$original_rr <- exp(dt$ln_rr)
  dt$original_se <- dt$original_rr * dt$ln_rr_se
  dt$pred_rr <- exp(dt$pred_ln_rr)
  dt$pred1_se <- dt$pred_rr * dt$pred1_ln_se
  dt$pred2_se <- dt$pred_rr * dt$pred2_ln_se
  
  # add input weights from each model ---------------------------------------
  
  dt$input_weights <- 1 / dt$ln_rr_se ^ 2
  
  # return mr-brt results ---------------------------------------------------
  
  return(
    list(
      mr_brt_df = dt,
      beta_gamma_df = coefs_dt,
      mr_brt_model = fit1
    )
  )
  
}
