###
# Post analysis related functions
###
library(scales)
source("./src/utils/egger_functions.R")

extract_data_info <- function(model,
                              ref_covs,
                              alt_covs,
                              exp_quantiles = c(0.15, 0.85),
                              exp_bounds = NULL,
                              pred_exp_bounds = NULL,
                              num_points = 100L) {
  # extract dataframe
  df <- model$data$to_df()
  
  # get exposure information
  if (is.null(ref_covs)) {
    ref_exp <- matrix(0, nrow = nrow(df), ncol = 1)
  } else {
    ref_exp <- df[ref_covs]
  }
  alt_exp <- df[alt_covs]
  df$ref_mid <- rowMeans(ref_exp)
  df$alt_mid <- rowMeans(alt_exp)
  
  exp_limits <- c(min(min(ref_exp), min(alt_exp)), max(max(ref_exp), max(alt_exp)))
  if (is.null(pred_exp_bounds)) {
    pred_exp <- seq(exp_limits[1], exp_limits[2], length = num_points)  
  } else {
    pred_exp <- seq(pred_exp_bounds[1], pred_exp_bounds[2], length = num_points)
  }
  
  if (is.null(exp_bounds)) {
    exp_bounds = c(quantile(df$ref_mid, exp_quantiles[1]), quantile(df$alt_mid, exp_quantiles[2]))
  }
  
  # create prediction of log relative risk
  pred_lrr <- pred_exp*model$beta_soln[[1]]
  
  # add prediction column
  df$pred <- (df$alt_mid - df$ref_mid)*model$beta_soln[[1]]
  
  # add outlier column
  df$outlier <- model$get_w_soln() < 0.1
  
  # create prediction at reference and alternative mid points
  df$ref_pred <- (df$ref_mid - exp_limits[1])*model$beta_soln[[1]]
  df$alt_pred <- df$ref_pred + df$obs
  df$residual <- df$obs - df$pred
  df$residual_se <- sqrt(df$obs_se^2 + df$pred^2*model$gamma_soln[[1]])
  
  list(
    df = df,
    ref_covs = ref_covs,
    alt_covs = alt_covs,
    exp_limits = exp_limits,
    exp_bounds = exp_bounds,
    exp_quantiles = exp_quantiles,
    pred_exp = pred_exp,
    pred_lrr = pred_lrr
  )
}


get_df_fill <- function(df) {
  # get number of filled data points
  absr_rank <- rep(0L, length=nrow(df))
  absr_rank[order(abs(df$residual))] <- seq(1, nrow(df))
  sort_index <- order(df$residual, decreasing = mean(df$residual/df$residual_se) > 0)
  num_fill <- nrow(df) - absr_rank[tail(sort_index, n=1)]
  
  # get the fill-dataframe
  df_fill <- df[sort_index[1:num_fill],]
  df_fill$study_id <- paste0('fill_', df_fill$study_id)
  df_fill$residual <- - df_fill$residual
  df_fill$obs <- df_fill$pred + df_fill$residual
  df_fill$alt_pred <- df_fill$obs + df_fill$ref_pred
  df_fill
}

plot_residual <- function(df, title){
  residual <- df$residual
  obs_se <- df$residual_se
  max_obs_se <- quantile(obs_se, 0.99)
  fill_index <- rownames(df)[grepl('fill', df$study_id, fixed=TRUE)]
  # funnel plot
  plot(residual, obs_se, pch=19, col=alpha('gray', 0.6), 
       ylim=c(max_obs_se,0), xlim=c(-2*max_obs_se, 2*max_obs_se),
       yaxs='i', xlab="residual", ylab="ln_rr_se", main=title)
  if (length(fill_index) > 0){
    points(df[fill_index,'residual'], df[fill_index,'residual_se'], 
           col=alpha('#008080',0.6), pch=19,)
  }
  
  if (sum(df$outlier) > 0) {
    outlier_index <- rownames(df)[df$outlier==TRUE]
    points(df[outlier_index, 'residual'], df[outlier_index, 'residual_se'],
           col=alpha('red', 0.6), pch=4)
  }
  x <- c(0, -1.96*max_obs_se, 1.96*max_obs_se)
  y <- c(0, max_obs_se, max_obs_se)
  polygon(x,y,col=alpha('#B0E0E6', 0.4))
  lines(c(0, 0), c(0, max_obs_se), lty='dashed')
  lines(c(0, -1.96*max_obs_se), c(0, max_obs_se), col='#87CEFA')
  lines(c(0, 1.96*max_obs_se), c(0, max_obs_se), col='#87CEFA')
}

get_gamma_sd <- function(model){
  gamma <- model$gamma_soln
  gamma_fisher <- model$lt$get_gamma_fisher(gamma)
  return(1/sqrt(gamma_fisher[1,1]))
}

get_beta_sd <- function(model){
  model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
  return(1/sqrt(beta_hessian[1,1]))
}

get_uncertainty_info <- function(data_info, model) {
  # extract info
  beta <- model$beta_soln[[1]]
  gamma <- model$gamma_soln[[1]]
  beta_sd <- get_beta_sd(model)
  gamma_sd <- get_gamma_sd(model)
  
  # compute uncertainty of fixed effects
  inner_fe_sd <- sqrt(beta_sd^2 + gamma)
  outer_fe_sd <- sqrt(beta_sd^2 + gamma + 2*gamma_sd)
  inner_betas <- c(qnorm(0.05, mean=beta, sd=inner_fe_sd),
                   qnorm(0.95, mean=beta, sd=inner_fe_sd))
  outer_betas <- c(qnorm(0.05, mean=beta, sd=outer_fe_sd),
                   qnorm(0.95, mean=beta, sd=outer_fe_sd))
  
  df <- data_info$df
  
  # compute the lower and upper draws
  pred_exp <- data_info$pred_exp
  pred_lrr <- data_info$pred_lrr
  inner_draws <- rbind(inner_betas[1]*pred_exp, inner_betas[2]*pred_exp) 
  outer_draws <- rbind(outer_betas[1]*pred_exp, outer_betas[2]*pred_exp)
  # compute and score
  pred_exp <- data_info$pred_exp
  exp_bounds <- data_info$exp_bounds
  index <- (pred_exp >= exp_bounds[1]) & (pred_exp <= exp_bounds[2])
  sign_score <- sign(mean(pred_lrr[index]))
  score <- min(rowMeans(outer_draws[,index])*sign_score)
  list(
    score = score,
    inner_draws = inner_draws,
    outer_draws = outer_draws
  )
}

plot_model <- function(data_info,
                       uncertainty_info,
                       model,
                       uncertainty_info_fill = NULL,
                       model_fill = NULL) {
  df <- data_info$df
  main <- paste0(ro_pair, ": score=", round(uncertainty_info$score,3))
  pred_exp <- data_info$pred_exp
  pred_lrr <- data_info$pred_lrr
  inner_draws <- uncertainty_info$inner_draws
  outer_draws <- uncertainty_info$outer_draws
  xrange <- max(pred_exp) - min(pred_exp)
  
  # plot data
  plot(df$alt_mid, df$alt_pred,
       col = alpha('gray', 0.8), cex = 0.05/df$obs_se,
       xlim = c(min(pred_exp) - xrange/20, max(pred_exp) + xrange/20),
       xlab = "", ylab = "", pch = 19)
  df_outlier <- df[df$outlier,]
  points(df_outlier$alt_mid, df_outlier$alt_pred,
         cex = 0.05/df_outlier$obs_se, col = alpha('red', 0.9), pch = 4)

  # plot prediction
  lines(pred_exp, pred_lrr, col="#008080", lwd=1)
  
  # plot uncertainties
  polygon(c(pred_exp, rev(pred_exp)), c(inner_draws[1,], rev(inner_draws[2,])),
          col=alpha("#69b3a2", 0.2), border=FALSE)
  polygon(c(pred_exp, rev(pred_exp)), c(outer_draws[1,], rev(outer_draws[2,])),
          col=alpha("#69b3a2", 0.3), border=FALSE)
  
  # plot filled model
  if (!is.null(model_fill)){
    lines(pred_exp, uncertainty_info_fill$outer_draws[1,], lty='dashed', col=alpha('gray', 0.9))
    lines(pred_exp, uncertainty_info_fill$outer_draws[2,], lty='dashed', col=alpha('gray', 0.9))
    df_fill <- df[grepl('fill', df$study_id, fixed=TRUE),]
    points(df_fill$alt_mid, df_fill$alt_pred,
           cex = 0.05/df_fill$obs_se, pch=1, col = alpha("#008080", 0.8))
    main <- paste0(main, ", score_fill=", round(uncertainty_info_fill$score, 3))
  }
  
  # plot bounds
  abline(v = data_info$exp_bounds[1], lty = 'dashed', lwd = 1, col = 'black')
  abline(v = data_info$exp_bounds[2], lty = 'dashed', lwd = 1, col = 'black')
  # plot 0 line
  abline(h = 0.0, lwd = 1, col = 'black')
  title(main = main) 
}

summarize_model <- function(data_info,
                            uncertainty_info,
                            model,
                            egger_model,
                            egger_model_all,
                            uncertainty_info_fill = NULL,
                            model_fill = NULL) {
  num_fill <- sum(grepl('fill', data_info$df$study_id, fixed=TRUE))
  summary <- list(
    ro_pair = data_info$ro_pair,
    num_fill = num_fill, 
    gamma = model$gamma_soln[[1]], 
    gamma_sd = get_gamma_sd(model),
    score = uncertainty_info$score
  )
  
  if (!is.null(model_fill)){
    summary$gamma_adjusted <- model_fill$gamma_soln[[1]]
    summary$gamma_sd_adjusted <- get_gamma_sd(model_fill)
    summary$score_adjusted <- uncertainty_info_fill$score
  }
  
  for (name in names(egger_model)){
    summary[[paste0("egger_",name)]] <- egger_model[[name]]
  }
  
  for (name in names(egger_model_all)){
    summary[[paste0("egger_all_",name)]] <- egger_model_all[[name]]
  }
  
  data.frame(summary)
}

get_draws <- function(data_info, model, num_draws = 1000L) {
  # extract info
  beta <- model$beta_soln[[1]]
  gamma <- model$gamma_soln[[1]]
  beta_sd <- get_beta_sd(model)
  gamma_sd <- get_gamma_sd(model)
  
  outer_fe_sd <- sqrt(beta_sd^2 + gamma + 2*gamma_sd)
  beta_samples <- rnorm(num_draws, mean = beta, sd = outer_fe_sd)
  draws <- as.data.frame(
    cbind(data_info$pred_exp, outer(data_info$pred_exp, beta_samples))
  )
  names(draws) <- c("exposure", sapply(1:num_draws, function(i) paste0("draw_", i)))
  draws
}