#--------------------------------------------------------------------------------------
# Purpose: useful functions for age specific RR of smoking on CVDs
# Author: Xiaochen Dai, adapted from Haley's codes
# Data: 07/25/2022
#--------------------------------------------------------------------------------------

# Setup ---------------------------
source("get_age_metadata.R")
source("get_ids.R")
library(data.table)
library(ggplot2)
rei_table <- get_ids("rei")
cause_table <- get_ids("cause")

apply_age_pattern <- function(ro_pair, risk_curve_draws_df, age_pattern_draws_df, 
                              plot = T, 
                              plot_path){
  # Will apply the age pattern to an 'all-age' exposure-dependent risk curve
  #
  #  @causeid = integer cause_id for the RO pair of interest
  #  @reid = integer rei_id for the RO pair of interest
  #  @risk_curve_draws_df = dataframe of MR-BRt risk curve draws with columns 'exposure' and 'draw_0'...'draw_999'. 
  #                     Should not be age specific
  #  @age_pattern_draws_df = dataframe of age pattern attenuation factor draws. Either as returned from the 'weighted_average'
  #                     or 'simple_average' function, or read in directly from the draw csv
  #  @draws_in_log = boolean; T if risk curves draws are in log space, F if in normal space; default F
  #  @return_draws_log = boolean; T if you want to return and plot risk curve in log space, F if in normal space; default F
  
  if(sum(!c("exposure", paste0("draw_",0:999)) %in% colnames(risk_curve_draws_df))){
    stop("Risk curve draws should have 'exposure' and 'draw_1'...'draw_1000' as column names, apparently it does not")
  }
  if( sum(colnames(risk_curve_draws_df) %like% "age")){
    stop(paste0("Risk curve draws should not be age specific, remove column(s): ", paste0(
      colnames(risk_curve_draws_df)[colnames(risk_curve_draws_df) %like% "age"], collapse = ","
    )))
  }
  
  risk_byvars <- colnames(risk_curve_draws_df)[! colnames(risk_curve_draws_df) %like% "draw"]
  
  # reshape both df long
  risk_curve <- melt(risk_curve_draws_df, id.vars = risk_byvars, value.name = "rr")
  age_pct <- melt(age_pattern_draws_df, id.vars = "age_group_id", measure.vars = paste0("draw_", 0:999), value.name = "pct_change")
  
  # merge together and apply attenuation factor (the attenuation factor is on the RR scale)
  risk_curve <- merge(risk_curve, age_pct, by = "variable", allow.cartesian = T)
  risk_curve[, rr:= exp(rr)]
  
  risk_curve[, age_specific_rr := (rr-1)*pct_change+1]
  risk_curve[, `:=` (rr = NULL, pct_change = NULL)]
  risk_curve<- as.data.table(dcast(risk_curve, ... ~ variable, value.var = "age_specific_rr"))
  
  if(plot){
    ylab <- "RR"
    path <- plot_path
    visualize_age_specific_curve(risk_curve, ylab, plot_name = paste0("smoking-", ro_pair,", draws of AFs"), plot_path=path)
  }
  
  return(risk_curve)
}

apply_age_pattern_mean_af <- function(ro_pair, risk_curve_draws_df, age_pattern_mean_df, age_pattern_mean_log_df, log_af=T, plot = T, plot_path){
  # Will apply the age pattern to an 'all-age' exposure-dependent risk curve
  #
  #  @causeid = integer cause_id for the RO pair of interest
  #  @reid = integer rei_id for the RO pair of interest
  #  @risk_curve_draws_df = dataframe of MR-BRt risk curve draws with columns 'exposure' and 'draw_0'...'draw_999'. 
  #                     Should not be age specific
  #  @age_pattern_draws_df = dataframe of age pattern attenuation factor draws. Either as returned from the 'weighted_average'
  #                     or 'simple_average' function, or read in directly from the draw csv
  #  @log_af = boolean; T if the AFs are based on log_rr, F if the AFs are based on rr ; default T
  
  if(sum(!c("exposure", paste0("draw_",0:999)) %in% colnames(risk_curve_draws_df))){
    stop("Risk curve draws should have 'exposure' and 'draw_1'...'draw_1000' as column names, apparently it does not")
  }
  if( sum(colnames(risk_curve_draws_df) %like% "age")){
    stop(paste0("Risk curve draws should not be age specific, remove column(s): ", paste0(
      colnames(risk_curve_draws_df)[colnames(risk_curve_draws_df) %like% "age"], collapse = ","
    )))
  }
  
  risk_byvars <- colnames(risk_curve_draws_df)[! colnames(risk_curve_draws_df) %like% "draw"]
  
  # reshape risk draws and merge with template
  risk_curve <- melt(risk_curve_draws_df, id.vars = risk_byvars, value.name = "rr")
  template <-  expand.grid(exposure = 0:100, variable = paste0("draw_", 0:999), age_group_id = c(9:20, 30:32, 235)) %>% data.table
  risk_curve <- merge(template, risk_curve, by=c("exposure", "variable")) # draws of log_rr
  risk_curve[, log_rr:= rr]
  risk_curve[, rr:= exp(log_rr)]
  
  # merge mean AF and apply attenuation factor
  if(log_af==T){
    risk_curve <- merge(risk_curve, age_pattern_mean_log_df, by="age_group_id")
    # apply the log_af on log_rr
    risk_curve[, age_specific_rr := exp(log_rr*pct_mean)]
  } else {
    risk_curve <- merge(risk_curve, age_pattern_mean_df, by="age_group_id")
    # apply the original scale AF on the log_rr
    risk_curve[, age_specific_rr := rr*pct_mean-(pct_mean-1)] 
  }
  
  risk_curve<- as.data.table(dcast(risk_curve, exposure + age_group_id + age_start + age_end + grp ~ variable, value.var = "age_specific_rr"))
  
  if(plot){
    ages <- get_age_metadata(19)
    risk_curve <- merge(risk_curve, ages[,.(age_group_id, age_group_name)], by = "age_group_id")
    risk_curve[, rrmean := apply(.SD, 1, median), .SDcols = paste0("draw_", 0:999)]
    risk_curve[, rrupper := apply(.SD, 1, quantile, 0.975), .SDcols = paste0("draw_", 0:999)]
    risk_curve[, rrlower := apply(.SD, 1, quantile, 0.025), .SDcols = paste0("draw_", 0:999)]
    
    ylab <- "RR"
    
    pdf(plot_path)
    p <- ggplot(risk_curve, aes(x = exposure, y = rrmean, ymin=rrlower, ymax = rrupper))+
      geom_line()+geom_ribbon(alpha = 0.2)+theme_bw()+ geom_abline(slope = 0, intercept = 1, color = "blue", alpha = 1/5)+
      facet_wrap(~age_group_name)+labs(y= ylab, title = paste0("smoking-", ro_pair,", mean AF"))
    print(p) 
    dev.off()
    
    print(p)
  }
  
  return(risk_curve)
}

visualize_age_specific_curve <- function(risk_curve, ylab, plot_name, plot_path){
  # For use in plotting age specific curves
  #
  #  @risk_curve = dataframe of age specific curves. Must have 'draw_0'...'draw_999' and 'age_group_id'
  #  @ylab = character to be used as y-axis label in plot. ex: RR or log(RR)
  #  @plot_name = character to be used as plot title in plot. 
  
  ages <- get_age_metadata(19)
  risk_curve <- merge(risk_curve, ages, by = "age_group_id")
  risk_curve[, rrmean := apply(.SD, 1, median), .SDcols = paste0("draw_", 0:999)]
  risk_curve[, rrupper := apply(.SD, 1, quantile, 0.975), .SDcols = paste0("draw_", 0:999)]
  risk_curve[, rrlower := apply(.SD, 1, quantile, 0.025), .SDcols = paste0("draw_", 0:999)]
  
  pdf(plot_path)
  plot <- ggplot(risk_curve, aes(x = exposure, y = rrmean, ymin=rrlower, ymax = rrupper))+
    geom_line()+geom_ribbon(alpha = 0.2)+theme_bw()+ geom_abline(slope = 0, intercept = 1, color = "blue", alpha = 1/5)+
    facet_wrap(~age_group_name)+labs( y= ylab, title = plot_name)
  print(plot) 
  dev.off()
  print(plot)
}