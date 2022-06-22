# unpack information from results and create draws
rm(list = ls())
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")

# define functions
# =================================================================================================
get_cov_names <- function(signal_model) {
  cov_model <- signal_model$sub_models[[1]]$cov_models[[1]]
  list(alt_covs = cov_model$alt_cov,
       ref_covs = cov_model$ref_cov)
}

get_risk_limits <- function(signal_model) {
  cov_names <- get_cov_names(signal_model)
  risk_data <- signal_model$data$get_covs(unlist(cov_names))
  c(min(risk_data), max(risk_data))
}

get_signal <- function(signal_model, risk) {
  cov_names <- get_cov_names(signal_model)
  risk_limits <- get_risk_limits(signal_model)
  df_covs <- data.frame(
    c(sapply(cov_names$ref_covs, function(x) rep(risk_limits[1], length.out = length(risk)),
             simplify = FALSE, USE.NAMES = TRUE),
      sapply(cov_names$alt_covs, function(x) risk,
             simplify = FALSE, USE.NAMES = TRUE))
  )
  data <- MRData()
  data$load_df(df_covs, col_covs=unlist(cov_names))
  signal_model$predict(data)
}

get_beta <- function(linear_model) {
  beta <- linear_model$beta_soln
  names(beta) <- linear_model$cov_names
  specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(linear_model)
  beta_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(specs)
  beta_sd <- 1/sqrt(diag(beta_hessian))
  names(beta_sd) <- linear_model$cov_names
  c(beta["signal"], beta_sd["signal"])
}

get_gamma <- function(linear_model) {
  gamma <- linear_model$gamma_soln[[1]]
  gamma_fisher <- linear_model$lt$get_gamma_fisher(linear_model$gamma_soln)
  gamma_sd <- 1/sqrt(diag(gamma_fisher)[[1]])
  c(gamma, gamma_sd)
}

get_soln <- function(linear_model) {
  list(
    beta_soln = get_beta(linear_model),
    gamma_soln = get_gamma(linear_model)
  )
}

get_ln_rr_draws <- function(signal_model,
                            linear_model,
                            risk,
                            num_draws = 1000L,
                            normalize_to_tmrel = FALSE,
                            include_re = TRUE) {
  # set seed inside function
  set.seed(1234)
  
  signal <- get_signal(signal_model, risk)
  re_signal <- signal
  soln <- get_soln(linear_model)
  
  fe_samples <- rnorm(num_draws, mean=soln$beta[1], sd=soln$beta[2])
  re_samples <- rnorm(num_draws, mean=0, sd=sqrt(soln$gamma[1] + 2*soln$gamma[2]))
  
  draws <- outer(signal, fe_samples)
  if (include_re) {
    draws <- draws + outer(re_signal, re_samples)
  }
  
  if (normalize_to_tmrel) {
    tmrel_index <- which.min(signal)
    draws <- apply(draws, 2, function(x) x - x[tmrel_index])
  }
  
  df <- as.data.frame(cbind(risk, draws))
  names(df) <- c("risk", sapply(1:num_draws, function(i) paste0("draw_", i)))
  return(df)
}

# process results
# =================================================================================================
# load models
signal_model_path <- "/mnt/team/msca/pub/archive/evidence-score/gbd2020-process/nuts_ihd/signal_model.pkl"
linear_model_path <- "/mnt/team/msca/pub/archive/evidence-score/gbd2020-process/nuts_ihd/new_linear_model.pkl"

signal_model <- py_load_object(filename = signal_model_path, pickle = "dill")
linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")

# specify risk, you need to input the exposures that you want to predict
risk <- 0:100

# get_draws
df <- get_ln_rr_draws(signal_model,
                      linear_model,
                      risk,
                      num_draws = 1000L,
                      normalize_to_tmrel = FALSE)

# visual check draws
draws <- df[, 2:ncol(df)]
draw_mean <- apply(draws, 1, function(x) mean(x))
draw_lower <- apply(draws, 1, function(x) quantile(x, probs=.05))
draw_upper <- apply(draws, 1, function(x) quantile(x, probs=.95))

lines(risk, draw_mean)
lines(risk, draw_lower)
lines(risk, draw_upper)
