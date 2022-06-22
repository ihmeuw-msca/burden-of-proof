###
# Post analysis related functions
###
library(scales)
source("./src/utils/egger_functions.R")

extract_data_info <- function(model,
                              cont_cov,
                              cont_quantiles = c(0.15, 0.85),
                              cont_bounds = NULL,
                              num_points = 100L) {
  # get observation information
  df <- model$data$to_df()
  df$pred <- model$predict(model$data)
  df$residual <- df$obs - df$pred
  df$residual_se <- sqrt(df$obs_se^2 + model$gamma_soln[[1]])
  df$outlier <- model$get_w_soln() < 0.1
  
  # get continuous covariate information and prediction
  cov <- df[[cont_cov]]
  if (is.null(cont_bounds)) {
    cont_bounds <- quantile(cov, cont_quantiles)
  }
  cont_limits <- c(min(cov), max(cov))
  pred_cov <- seq(cont_limits[1], cont_limits[2], length = num_points)
  
  df_pred <- data.frame(pred_cov)
  names(df_pred) <- cont_cov
  for (name in model$cov_names) {
    if (name == "intercept") {
      df_pred[name] <- 1
    } else if (name != cont_cov) {
      df_pred[name] <- 0
    }
  }
  data <- MRData()
  data$load_df(df_pred, col_covs = as.list(model$cov_names))
  pred_lrr <- model$predict(data)
  
  list(
    df = df,
    cont_cov = cont_cov,
    cont_limits = cont_limits,
    cont_bounds = cont_bounds,
    pred_cov = pred_cov,
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
  df_fill
}

get_gamma_sd <- function(model){
  gamma <- model$gamma_soln
  gamma_fisher <- model$lt$get_gamma_fisher(gamma)
  return(1/sqrt(gamma_fisher[1, 1]))
}

get_beta_sd <- function(model){
  model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
  beta_sd <- sqrt(diag(solve(beta_hessian)))
  names(beta_sd) <- model$cov_names
  return(beta_sd)
}


get_uncertainty_info <- function(data_info, model){
  beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
  gamma_sd <- get_gamma_sd(model)
  gamma_inner_samples <- matrix(rep(model$gamma_soln[1], 1000L), nrow = 1000L)
  gamma_outer_samples <- matrix(rep(model$gamma_soln[1] + 2*gamma_sd, 1000L), nrow = 1000L)
  
  df_pred <- data.frame(data_info$pred_cov)
  names(df_pred) <- data_info$cont_cov
  for (name in model$cov_names) {
    if (name == "intercept") {
      df_pred[name] <- 1
    } else if (name != data_info$cont_cov) {
      df_pred[name] <- 0
    }
  }
  data <- MRData()
  data$load_df(df_pred, col_covs = as.list(model$cov_names))
  
  inner_draws <- model$create_draws(data,
                                    beta_samples = beta_samples,
                                    gamma_samples = gamma_inner_samples)
  outer_draws <- model$create_draws(data,
                                    beta_samples = beta_samples,
                                    gamma_samples = gamma_outer_samples)
  inner_draws <- apply(inner_draws, 1, function(x) quantile(x, c(0.05, 0.95)))
  outer_draws <- apply(outer_draws, 1, function(x) quantile(x, c(0.05, 0.95)))
  
  index <- (data_info$pred_cov >= data_info$cont_bounds[1]) & (data_info$pred_cov <= data_info$cont_bounds[2])
  sign_score <- sign(mean(data_info$pred_lrr[index]))
  score <- 0.5*min(rowMeans(outer_draws[,index])*sign_score)
  list(
    score = score,
    inner_draws = inner_draws,
    outer_draws = outer_draws
  )
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

plot_model <- function(data_info,
                       uncertainty_info,
                       model,
                       uncertainty_info_fill = NULL,
                       model_fill = NULL) {
  main <- paste0(ro_pair, ": score=", round(uncertainty_info$score,3))
  df <- data_info$df
  pred_cov <- data_info$pred_cov
  pred_lrr <- data_info$pred_lrr
  inner_draws <- uncertainty_info$inner_draws
  outer_draws <- uncertainty_info$outer_draws
  xrange <- max(pred_cov) - min(pred_cov)
  
  # plot data
  plot(df[[data_info$cont_cov]], df$obs,
       col = alpha('gray', 0.8), cex = 0.05/df$obs_se,
       xlim = c(min(pred_cov) - xrange/20, max(pred_cov) + xrange/20),
       xlab = "", ylab = "", pch = 19)
  df_outlier <- df[df$outlier,]
  points(df_outlier[[data_info$cont_cov]], df_outlier$obs,
         cex = 0.05/df_outlier$obs_se, col = alpha('red', 0.9), pch = 4)
  
  # plot prediction
  lines(pred_cov, pred_lrr, col="#008080", lwd=1)
  
  # plot uncertainties
  polygon(c(pred_cov, rev(pred_cov)), c(inner_draws[1,], rev(inner_draws[2,])),
          col=alpha("#69b3a2", 0.2), border=FALSE)
  polygon(c(pred_cov, rev(pred_cov)), c(outer_draws[1,], rev(outer_draws[2,])),
          col=alpha("#69b3a2", 0.3), border=FALSE)
  
  # plot filled model
  if (!is.null(model_fill)){
    lines(pred_cov, uncertainty_info_fill$outer_draws[1,], lty='dashed', col=alpha('gray', 0.9))
    lines(pred_cov, uncertainty_info_fill$outer_draws[2,], lty='dashed', col=alpha('gray', 0.9))
    df_fill <- df[grepl('fill', df$study_id, fixed=TRUE),]
    points(df_fill[[data_info$cont_cov]], df_fill$obs,
           cex = 0.05/df_fill$obs_se, pch=1, col = alpha("#008080", 0.8))
    main <- paste0(main, ", score_fill=", round(uncertainty_info_fill$score, 3))
  }
  
  # plot bounds
  abline(v = data_info$cont_bounds[1], lty = 'dashed', lwd = 1, col = 'black')
  abline(v = data_info$cont_bounds[2], lty = 'dashed', lwd = 1, col = 'black')
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
  beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(num_draws, model)
  gamma_sd <- get_gamma_sd(model)
  gamma_outer_samples <- matrix(rep(model$gamma_soln[1] + 2*gamma_sd, num_draws), nrow = num_draws)
  
  df_pred <- data.frame(data_info$pred_cov)
  names(df_pred) <- data_info$cont_cov
  for (name in model$cov_names) {
    if (name == "intercept") {
      df_pred[name] <- 1
    } else if (name != data_info$cont_cov) {
      df_pred[name] <- 0
    }
  }
  data <- MRData()
  data$load_df(df_pred, col_covs = as.list(model$cov_names))
  outer_draws <- model$create_draws(data,
                                    beta_samples = beta_samples,
                                    gamma_samples = gamma_outer_samples)
  
  draws <- as.data.frame(
    cbind(data_info$pred_cov, outer_draws)
  )
  names(draws) <- c(data_info$cont_cov, sapply(1:num_draws, function(i) paste0("draw_", i)))
  draws
}