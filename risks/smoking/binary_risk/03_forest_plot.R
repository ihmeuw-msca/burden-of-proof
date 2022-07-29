####################################################################################################################################################################
# 
# Author: Xiaochen Dai
# Purpose: Plot mr-brt results
#
####################################################################################################################################################################

rm(list=ls())

library(data.table)
library(dplyr)
library(openxlsx)
library(ggplot2)
# library(crosswalk, lib.loc = "/ihme/code/mscm/R/packages/")
library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")

args <- commandArgs(trailingOnly = TRUE)

## NEED TO CHANGE THE RO_PAIR HERE

if(interactive()){
  # NOTE: the ro_pair for this script does not include age-specific info
  ro_pair <- "fractures" # only works for fractures now
  cov_setting <- "cov_finder_no_sex" # option: ['cov_finder', 'cov_finder_no_sex', 'no_cov','percent_male_only','self_selected'(no percent_male)]
  trim <- 0.9
  #out_dir <- "/ihme/homes/xdai88/gbd_tobacco/gbd2019_alcohol/evidence_score/testing/test_run1_2020_09_05/"
  out_dir <- "/mnt/team/team/pub/sub_risks/tobacco/code/xdai88/gbd2020_smoking/relative_risk_curves/binary_risk/fracture_binary/"
} else {
  ro_pair <- args[1]
  cov_setting <- args[2]
  trim <- args[3]
  out_dir <- args[4]
}

# 1. plotting the results ---------------------------------------------------------------------------------------------------------------------------------------------------

#obs_data <- fread(paste0('/ihme/homes/xdai88/gbd_tobacco/gbd2020_smoking/test_run_2020_10_25/fracture_binary/mrbrt_output_', cov_setting, '_', trim, '.csv'))
#mod1 <- py_load_object(filename = paste0('/ihme/homes/xdai88/gbd_tobacco/gbd2020_smoking/test_run_2020_10_25/fracture_binary/', ro_pair, '_', cov_setting,'_', trim, '.pkl'), pickle = "dill")

obs_data <- fread(paste0(out_dir, ro_pair,"_",cov_setting, '_', trim, '.csv'))
mod1 <- py_load_object(filename = paste0(out_dir, ro_pair, '_', cov_setting,'_', trim, '.pkl'), pickle = "dill")

cov_names <- mod1$cov_names[!mod1$cov_names=="intercept"]
header <- paste0(as.character(mod1$data),"\ncovariates: ",paste0(cov_names, collapse=", ") ,"\nexp(beta): ", round(exp(mod1$beta_soln[1]), digits = 3),"     gamma: ", mod1$gamma_soln)

# 
test <- obs_data[!is.na(se)]
test[, val := exp(val)]; test[, lower := exp(lower)]; test[, upper := exp(upper)]
results <- obs_data[is.na(se)]

#lin_effect <- as.data.table(delta_transform(mean=test$val, sd=test$se, transformation='log_to_linear'))
#data <- cbind(lin_effect, test)
#data[, upper:=mean_linear+1.96*sd_linear]
#data[, lower:=mean_linear-1.96*sd_linear]

plot_data <- rbind(test, results, fill=T)
plot_data[, mean_linear:=val]
plot_data[study=='2019 Result', data:=3]

plot_data[included==0, data:=4]

# forest plot of the data
color_vals <- c("black", "blue", "darksalmon", "red")
names(color_vals) <- c(1,2,3,4)

if (cov_setting %in% c('cov_finder', 'percent_male_only')) {
  
  plot_data[sample_sex==0, sex:='female']
  plot_data[sample_sex==1, sex:='male']
  plot_data[!sample_sex%in%c(0,1) & !is.na(sample_sex), sex:='both']
  
  
  alpha_vals <- c(1, 0.75, 0.75, 1)
  names(alpha_vals) <- c('NA', 'male', 'female', 'both')
  
  
  pdf(paste0(out_dir, ro_pair, '_simple_forest_plot_', cov_setting, '_', trim, '.pdf'),
      height = 12)  
  
  p <- ggplot(data=plot_data,
              aes(x = row,y = mean_linear, ymin = lower, ymax = upper, color = as.factor(data)))+
    geom_pointrange(aes(shape = as.factor(data)))+
    scale_color_manual("", values = color_vals, guide = F) + 
    # scale_alpha_manual(values=c('NA'=1, 'male'=0.3, 'female'=0.3, 'both'=1)) + 
    geom_hline(yintercept =1, linetype=2, color = "black")+
    xlab('study id')+ ylab(paste0('Relative Risk', " (95% Confidence Interval)"))+
    scale_x_continuous(label = obs_data$study, breaks = obs_data$row)+
    geom_errorbar(aes(ymin=lower, ymax=upper),width=0.5,cex=1)+ 
    labs(subtitle = header, color = "")+
    scale_shape_discrete(guide = F)+
    theme_bw()+
    coord_flip()
  
  print(p)
  
  dev.off()
} else {
  pdf(paste0(out_dir, ro_pair, '_simple_forest_plot_', cov_setting, '_', trim, '.pdf'),
      height=15)  
  
  p <- ggplot(data=plot_data,
              aes(x = row,y = mean_linear, ymin = lower, ymax = upper, color = as.factor(data)))+
    geom_pointrange(aes(shape = as.factor(data)))+
    scale_color_manual("", values = color_vals, guide = F) + 
    geom_hline(yintercept =1, linetype=2, color = "black")+
    xlab('study id')+ ylab(paste0('Relative Risk', " (95% Confidence Interval)"))+
    scale_x_continuous(label = obs_data$study, breaks = obs_data$row)+
    geom_errorbar(aes(ymin=lower, ymax=upper),width=0.5,cex=1)+ 
    labs(subtitle = header, color = "")+
    scale_shape_discrete(guide = F)+
    theme_bw()+
    coord_flip()
  
  print(p)
  
  dev.off()
}
