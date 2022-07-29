
library(dplyr)
library(data.table)
library(ggplot2)
library(reticulate)

summarize_draws <- function(data){
  
  df <- as.data.table(copy(data))
  draw_cols <- colnames(df)[grepl("draw_", colnames(df))]
  
  df[, mean := apply(.SD, 1, mean), .SDcols = draw_cols]
  df[, upper := apply(.SD, 1, quantile, 0.975), .SDcols = draw_cols]
  df[, lower := apply(.SD, 1, quantile, 0.025), .SDcols = draw_cols]
  
  df[, (draw_cols) := NULL]
  return(df)
  
}

exponentiate_draws <- function(df){
  # Function to exponentiate draws of ln(RR) so that each draw is in normal space
  # df - dataframe of wide draws where draws are named draw_0,...,draw_999 and some identifying columns
  
  id_vars <- setdiff(names(df), paste0("draw_", 0:999))
  df_long <- melt(df, id.vars = id_vars)
  
  df_long[, value := exp(value)]
  
  out <- as.data.table(dcast(df_long, ... ~ variable, value.var = "value")) # cast wide again
  
  return(out)
  
}

model_summary_table <- function(out_dir, ro_pair){
  
  WORK_DIR <- "/ihme/code/dbd/hkl1/diet_exposure/code/GBD2020/relative_risk/mrbrt_evidence_score/"
  source(paste0(out_dir, "/config.R"))
  
  # Load final model
  final_model <- py_load_object(filename=paste0(out_dir, "04_mixed_effects_pkl_files/", ro_pair, ".pkl"), 
                                pickle = "dill")
  
  input_data <- readRDS(paste0(out_dir, "/00_prepped_data/", ro_pair, ".RDS"))
  df <- data.table(input_data$df)
  data_used <- df[,.(num_observations = .N), by = c("nid", final_model$cov_names[!final_model$cov_names=="signal"])]
  studies_used <- data_used[,.(num_studies = .N, num_observations = sum(num_observations)),  by = c(final_model$cov_names[!final_model$cov_names=="signal"])]
  studies_used <- rbind(studies_used, data.table(model = "total in model", num_studies = length(unique(data_used$nid)), num_observations = nrow(df)), fill = T)
  
  # get model summary
  beta_result <- as.data.table(t(py_to_r(final_model$summary()[[1]])), keep.rownames = "parameter")
  setnames(beta_result, "0", "beta")
  gamma_result <- as.data.table(t(py_to_r(final_model$summary()[[2]])),keep.rownames = "parameter")
  setnames(gamma_result, "0", "gamma")
  
  sampling <- import("mrtool.core.other_sampling")
  num_samples <- 1000L
  beta_samples <- as.data.table(sampling$sample_simple_lme_beta(num_samples, final_model))
  setnames(beta_samples, beta_result$parameter)
  beta_samples_w <- melt(beta_samples, variable.name = "parameter")
  
  beta_ui <- beta_samples_w[, .(beta_lower = quantile(value, 0.025),
                                beta_upper = quantile(value, 0.975)), by = "parameter"]
  beta_result <- merge(beta_result, beta_ui, by = "parameter")
  
  parameter_results <- merge(beta_result, gamma_result, by = "parameter", all = T)

  return(list(parameter_results, final_model$data, data_used, studies_used))
}

generate_draws_nogamma <- function(signal_model, final_model, exposure){
  
  # Set up prediction data frame
  min_cov <- rep(min(exposure), length(exposure))
  df_signal_pred <- data.frame(a_0=min_cov, a_1=min_cov, b_0=exposure, b_1=exposure)
  
  # Predict using signal model and gridded exposure
  data_signal_pred <- MRData()
  data_signal_pred$load_df(
    df_signal_pred,
    col_covs = as.list(c("a_0","a_1", "b_0", "b_1"))
  )
  df_signal_pred$signal <- signal_model$predict(data_signal_pred)
  
  df_final_pred <- as.data.table(df_signal_pred)
  
  # For any other bias covariates selected, we assume we should predict out at 0
  bias_covs <- final_model$cov_names[final_model$cov_names != "signal"]
  for(var in bias_covs){
    df_final_pred[, paste0(var):= 0]
  }
  
  data_final_pred <- MRData()
  data_final_pred$load_df(
    df_final_pred,
    col_covs = as.list(c("signal", bias_covs))
  )
  
  # Create draws and prediction
  sampling <- import("mrtool.core.other_sampling")
  num_samples <- 1000L
  beta_samples <- sampling$sample_simple_lme_beta(num_samples, final_model)
  gamma_samples <- rep(final_model$gamma_soln, num_samples) * matrix(1, num_samples)
  
  y_draws_fe = as.data.table(final_model$create_draws(data_final_pred, beta_samples, gamma_samples, random_study=F))
  setnames(y_draws_fe, colnames(y_draws_fe), paste0("draw_",0:999))
  
  # Resample if NA draws AND drop draws that are greater than 6 MADs
  y_draws_fe <- cbind(data.table("exposure" = exposure), as.data.table(y_draws_fe))
  y_draws_fe <- outlier_draws(y_draws_fe, exposure_col = "exposure")
  return(y_draws_fe)
}

min_exp_draws <- function(df, exposure_col = "b_0"){
  # Function to normalize draws of ln(RR) so that each draw is normalized by the mean exposure value with the lowest value
  # df - dataframe of wide draws with identifying columns cause_id, b_0 (or other exposure column name specified in exposure_col argument), 
  #      draws are named "draw_0"..."draw_999"
  #      each draw value in ln(RR)
  # return_log_space - Boolean to specify whether draws should be returned in log space or normal space
  # exposure_col - name of column that identifies the name of the alternative exposure column name
  # returns wide dataframe by draw
  id_vars <- setdiff(names(df), paste0("draw_", 0:999))
  df_long <- melt(df, id.vars = id_vars)
  df_long[, mean := mean(value), by = exposure_col]
  min_exp <- filter(df_long, mean == min(df_long$mean)) %>%
    select(exposure_col) %>%
    unique() %>%
    as.numeric()
  
  # min_exp is the TMREL
  
  temp <- df[get(exposure_col) == min_exp]
  temp <- melt(temp, id.vars = id_vars) %>% setnames(., "value", "min_val")
  
  df_long <- merge(df_long, temp[,.(variable, min_val)], by = "variable")
  rm(temp)
  
  if(return_log_space){
    df_long[, value_rescaled := value - min_val] # Take the min of each draw
    df_long[, c("min_val", "value", "mean") := NULL] # drop old columns
  } else {
    df_long[, value_rescaled := exp(value - min_val)] # Take the min of each draw
    df_long[, c("min_val", "value", "mean") := NULL] # drop old columns
  }
  
  
  out <- dcast(df_long, ... ~ variable, value.var = "value_rescaled") # cast wide again
  
  return(list(out, min_exp))
}

normalize_draws <- function(df, return_log_space = TRUE, exposure_col = "b_0"){
  # Function to normalize draws of ln(RR) so that each draw is normalized by the mean exposure value with the lowest value
  # df - dataframe of wide draws with identifying columns cause_id, b_0 (or other exposure column name specified in exposure_col argument), 
  #      draws are named "draw_0"..."draw_999"
  #      each draw value in ln(RR)
  # return_log_space - Boolean to specify whether draws should be returned in log space or normal space
  # exposure_col - name of column that identifies the name of the alternative exposure column name
  # returns wide dataframe by draw
  id_vars <- setdiff(names(df), paste0("draw_", 0:999))
  df_long <- melt(df, id.vars = id_vars)
  df_long[, mean := mean(value), by = exposure_col]
  min_exp <- filter(df_long, mean == min(df_long$mean)) %>%
    select(exposure_col) %>%
    unique() %>%
    as.numeric()
  
  # min_exp is the TMREL
  
  temp <- df[get(exposure_col) == min_exp]
  temp <- melt(temp, id.vars = id_vars) %>% setnames(., "value", "min_val")
  
  df_long <- merge(df_long, temp[,.(variable, min_val)], by = "variable")
  rm(temp)
  
  if(return_log_space){
    df_long[, value_rescaled := value - min_val] # Take the min of each draw
    df_long[, c("min_val", "value", "mean") := NULL] # drop old columns
  } else {
    df_long[, value_rescaled := exp(value - min_val)] # Take the min of each draw
    df_long[, c("min_val", "value", "mean") := NULL] # drop old columns
  }
  
  
  out <- dcast(df_long, ... ~ variable, value.var = "value_rescaled") # cast wide again
  
  return(list(out, min_exp))
}

outlier_tmrel <- function(vec, offset=0){
  
  median_vec <- median(vec$min_exp, na.rm = T)
  mad_vec <- mad(vec$min_exp, na.rm = T)+offset
  if(mad_vec ==0){
    warning("MAD of tmrel is 0, cannot proceed")
    return(list(vec, c(), c()))
  }
  
  drop <- as.numeric(abs(vec$min_exp - median_vec)/mad_vec > 6)
  
  if(sum(drop) >0 ){
    print(paste0("Dropping ", sum(drop), " draws that are greater than 6 MADs."))
  }
  
  # drop if draw
  drop_draws <- vec$variable[drop==1]
  vec <- vec[!(variable %in% drop_draws)]
  
  add_draws <- sample(vec$variable, sum(drop), replace = T)
  
  new_draws <- rbindlist(lapply(add_draws, function(x) return(vec[variable==x])))
  
  out <- rbind(vec, new_draws)
  
  return(list(out, add_draws, drop_draws))
  
}

outlier_draws <- function(dt, exposure_col = "b_0"){
  # Function to replace NA draws AND drop draws that are greater than 6 MADs
  # df - dataframe of wide draws with identifying columns cause_id, b_0 (or other exposure column name specified in exposure_col argument), 
  #      draws are named "draw_1"..."draw_1000"
  #      each draw value in ln(RR)
  # exposure_col - name of column that identifies the name of the alternative exposure column name
  # returns wide dataframe by draw in log space
  
  temp <- copy(dt)
  # Replacing NA columns with adjacent draw
  df <- temp %>% 
    dplyr::select_if(~ !all(is.na(.)))
  
  missing_draws <- setdiff(names(temp), names(df)) # gets names of missing draws
  if(length(missing_draws) >0 ){
    print(paste0("There are ", length(missing_draws), " missing draws in the dataframe."))
  }
  
  
  if(length(missing_draws) != 0){
    for(draw in missing_draws){
      df[, paste0(draw) := get(paste0('draw_', as.integer(gsub('draw_', '', draw)) + 1))]
    }
  }
  
  all_draws <- grep("draw", names(df), value = T)
  id_vars <- setdiff(names(df), all_draws)
  df_long <- melt(df, id.vars = id_vars)
  df_long[, median := lapply(.SD, function(x) median(x, na.rm = T)), .SDcols = "value", by = exposure_col]
  df_long[, mad := lapply(.SD, function(x) mad(x, na.rm = T)), .SDcols = "value", by = exposure_col]
  df_long[, drop := ifelse(abs(value - median)/mad > 6, 1, 0)]
  
  drop_draws <- as.character(unique(df_long[drop == 1]$variable))
  if(length(drop_draws) >0 ){
    print(paste0("Dropping ", length(drop_draws), " draws that are greater than 6 MADs."))
  }
  if(length(drop_draws) != 0){
    df[, c(drop_draws) := NULL]
    
    draw_options <- grep("draw", names(df), value = T)
    new_draws <- sample(draw_options, length(drop_draws), replace = T)
    
    add_draws <- df[, new_draws, with = F]
    colnames(add_draws) <- drop_draws
    
    out <- cbind(df, add_draws)
    
    return(out)
  } else {
    return(df)
  }
}

# Function to make effect modifier column from bias or dummy covariate
make.eff.mod <- function(data, columns, exp){
  # data is your dataframe (must be data table)
  # columns is a character vector of the names of columns you want to convert to effect modifiers
  # exp is a character object of the name of column that represents the dose
  
  for(c in 1:length(columns)){
    effect_mod <- data.table(var = data[, get(columns[c])] * data[, get(exp)]) %>% setnames(., "var", paste0("e_m_", columns[c]))
    data <- cbind(data, effect_mod)
  }
  return(data)
}

#####
plot_signal_model <- function(model, ro_pair){
  
  # Get dataframe
  data_df <- as.data.table(model$data$to_df())
  
  # Extract weights
  w <- t(do.call(rbind, 
                 lapply(1:length(model$sub_models), 
                        function(i){model$sub_models[[i]]$w_soln}))
  ) %*% model$weights
  
  data_df$weight <- w
  
  # Auxilary columns
  data_df[, `:=` (b_mid = (b_0 + b_1)/2, a_mid = (a_0 + a_1)/2)]
  
  # Predict out for full range of exposure
  min_val <- head(model$sub_models[[1]]$get_cov_model(model$ensemble_cov_model_name)$spline$knots, 1)
  max_val <- tail(model$sub_models[[1]]$get_cov_model(model$ensemble_cov_model_name)$spline$knots, 1)
  exposure <- seq(min_val, max_val, length.out = 1000)
  
  
  pred_df <- data.table(b_0 = exposure,
                        b_1 = exposure,
                        a_0 = min_val,
                        a_1 = min_val)
  
  # Make sure each a_mid value in the data has a prediction
  pred_df <- rbind(pred_df, data.table(b_0 = unique(data_df$a_mid), b_1 = unique(data_df$a_mid), a_0 = min_val, a_1 = min_val))
  
  pred_dat <- MRData()
  pred_dat$load_df(
    data = pred_df,
    col_covs = as.list(names(pred_df))
  )
  
  pred_df$prediction <- model$predict(pred_dat)
  pred_df$a_mid <- pred_df$b_0
  
  # Merge prediction onto dataframe
  data_df <- merge(data_df, pred_df[,.(prediction, a_mid)], by = "a_mid")
  data_df[, outlier := ifelse(weight >= 0.1, 0, 1)]
  
  # Plot it up
  gg <- ggplot(data = pred_df, aes(x = b_0, y = prediction)) +
    geom_point(data = data_df, aes(x = b_mid, y = prediction + obs, color = as.factor(outlier), shape = as.factor(outlier), size = 1/obs_se), alpha = 0.7)+
    geom_line(size = 1.2, color = "black") +
    scale_color_manual(values = c("grey55", "red")) +
    scale_shape_manual(values = c(16, 4))+
    theme_bw() + 
    ggtitle(paste0("Signal model for ", ro_pair))+
    xlab("Exposure")+
    ylab("ln(RR)")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(legend.position = "none", plot.title = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.2)), panel.border = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(gg)
}



# Function to rescale predictions
rescale_predictions <- function(df, mean, uncertainty = T, log_scale = F){
  # df is the prediction dataframe/matrix 
  # mean is the name of the column you want to rescale (will also rescale the lower and upper bounds)
  pred_df <- copy(df)
  min_rr <- min(pred_df[, c(mean), with = F])
  
  if(!log_scale){
    
    if(uncertainty){
      pred_df <- eval(parse(text=paste0("pred_df[, `:=` (", 
                                        mean, "_scaled = ", mean, "/", min_rr, ", ",
                                        mean, "_lo_scaled = ", mean, "_lo/", min_rr, ",",
                                        mean, "_hi_scaled = ", mean, "_hi/", min_rr,
                                        ")]")))
    } else {
      
      pred_df <- eval(parse(text=paste0("pred_df[, `:=` (", 
                                        mean, "_scaled = ", mean, "/", min_rr,
                                        ")]")))
    }
    
  } else {
    
    if(uncertainty){
      pred_df <- eval(parse(text=paste0("pred_df[, `:=` (", 
                                        mean, "_scaled = ", mean, "-", min_rr, ", ",
                                        mean, "_lo_scaled = ", mean, "_lo-", min_rr, ",",
                                        mean, "_hi_scaled = ", mean, "_hi-", min_rr,
                                        ")]")))
    } else {
      
      pred_df <- eval(parse(text=paste0("pred_df[, `:=` (", 
                                        mean, "_scaled = ", mean, "-", min_rr,
                                        ")]")))
    }
    
  }
  
  
  return(pred_df)
}


#####
plot_agetrend_model <- function(model, ro_pair){
  
  # Get dataframe
  data_df <- as.data.table(model$data$to_df())
  
  # Extract weights
  data_df$weight <- model$w_soln
  pred_df <- data.table("intercept"=c(1), "age_start" = seq(25, 99, by = 5), "age_end" = seq(30, 100, by = 5))
  
  pred_dat <- MRData()
  pred_dat$load_df(
    data = pred_df,
    col_covs = as.list(names(pred_df))
  )
  
  pred_df$prediction <- model$predict(pred_dat)
  pred_df[, age_midpoint := age_start + (age_end - age_start)/2 ]
  
  data_df[, outlier := ifelse(weight >= 0.1, 0, 1)]
  data_df[, age_midpoint := age_start + (age_end - age_start)/2 ]
  
  # Plot it up
  gg <- ggplot(data = pred_df, aes(x = age_midpoint, y = prediction)) +
    geom_point(data = data_df, aes(x = age_midpoint, y = obs, color = as.factor(study_id), shape = as.factor(outlier), size = 1/obs_se), alpha = 0.7)+
    geom_line(size = 1.2, color = "black") +
    scale_shape_manual(values = c(16, 4), guide = F)+
    scale_size_continuous(guide = F)+
    theme_bw() + 
    labs(title = paste0("Model fit for ", ro_pair), x = "Age", y = "Excess RR", color = "")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(legend.position = "bottom", plot.title = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.2)), panel.border = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(gg)
  
  gg <- ggplot(data = pred_df, aes(x = age_midpoint, y = exp(prediction))) +
    geom_point(data = data_df, aes(x = age_midpoint, y = exp(obs), color = as.factor(study_id), shape = as.factor(outlier), size = 1/obs_se), alpha = 0.7)+
    geom_line(size = 1.2, color = "black") +
    scale_shape_manual(values = c(16, 4), guide = F)+
    scale_size_continuous(guide = F)+
    theme_bw() + 
    labs(title = paste0("Model fit for ", ro_pair), x = "Age", y = "RR", color = "")+
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(legend.position = "bottom", plot.title = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.2)), panel.border = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
}