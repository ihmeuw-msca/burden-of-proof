#----------------------------------------------------------------------------------
# calculate across-cause risk curve for redmeat
#----------------------------------------------------------------------------------

rm(list = ls())

# System info
os <- Sys.info()[1]
user <- Sys.info()[7]

# Drives
j <- if (os == "Linux") "/home/j/" else if (os == "Windows") "J:/"
h <- if (os == "Linux") paste0("/homes/", user, "/") else if (os == "Windows") "H:/"

out_dir <- "FILEPATH" 
WORK_DIR <- "FILEPATH"
central_model_output_folder <- "FILEPATH" 

library(dplyr)
library(ggplot2)
library(data.table)
library(reticulate)
library(mrbrt002, lib.loc = "FILEPATH") 
source("/FILEPATH/get_ids.R")
source("/FILEPATH/get_draws.R") 
source("FILEPATH/model_functions.R")
ids <- get_ids("measure")
set.seed(123)

# Set up arguments---------------------------------------
random_effects <- T
meas_id <- 1 # 1: deaths, 2:DALYs
mrbrt_version <- "DATE_redmeat" 
risk_grep <- "redmeat" 
mrbrt_dir <- paste0(out_dir, mrbrt_version,"/")
save_dir <- paste0(out_dir, mrbrt_version, "/") 
exclude_ro_pair <- c()
note <- ""
sex_id <- 3
MAD_outlier_outcomes <- F
MAD_outlier_tmrel <- F

#----------------------------------------------------------
# Step 1: set up directory and cause list
#----------------------------------------------------------
if(meas_id == 1){
  source = "codcorrect"
  version_id = "version" 
   }else if(meas_id == 2){
  source = "dalynator"
  version_id = "version" 
}

meas_name <- gsub(" .*", "", ids[measure_id==meas_id, measure_name])

save_dir <- paste0(save_dir, "/06_tmrel/")
if(!dir.exists(save_dir)){dir.create(save_dir)}
save_dir <- paste0(save_dir, "/", risk_grep, "/")
dir.create(save_dir)

source(paste0(mrbrt_dir, "/config.R"))
#pdf(paste0(save_dir, "/", risk_grep, "all_cause_mortality_curve_", ifelse(random_effects, "re_", "fe_"), meas_name,note, ".pdf"), width = 11, height = 8)

pairs <- gsub(".pkl","",list.files(paste0(mrbrt_dir, "/04_mixed_effects_pkl_files/"), pattern = ".pkl")) 
diet_ro_pair_map <- fread("/FILEPATH/diet_ro_map.csv")
diet_ro_pair_map <- diet_ro_pair_map[include==1]

use_pairs <- pairs[grepl(risk_grep, pairs)]
causes <-  diet_ro_pair_map[risk_cause %in% use_pairs, .(risk_cause, cause_id)]
causes <- unique(causes[risk_cause %like% "hemstroke", cause_id:=494])
causes <- causes[!(risk_cause %in% exclude_ro_pair)]

print(causes)


#----------------------------------------------------------
# Step 2: Pull burden draws from GBD 2019
#----------------------------------------------------------
drop_cols <- c("location_id", "age_group_id", "year_id", "metric_id", "measure_id", "sex_id")
# First all but hem stroke since it has two subcauses
mort <- get_draws(gbd_id_type = "cause_id", 
                  gbd_id = causes$cause_id[!causes$cause_id %in% c(494, 496, 497)],
                  year_id = 2019,
                  measure_id = meas_id,
                  metric_id = 1,
                  location_id = 1,
                  age_group_id = 22,
                  sex_id = sex_id,
                  source = source,
                  release_id = 7,
                  version_id = version_id) %>%
  .[, c(drop_cols) := NULL] %>%
  melt(., id.vars = "cause_id")

if(494 %in% causes$cause_id){ 
  
  # Pulling draws for children of hemorrhagic stroke
  # Intracerebral hemorrhage
  intr_mort <- get_draws(gbd_id_type = "cause_id",
                         gbd_id = 496,
                         year_id = 2019,
                         measure_id = meas_id,
                         metric_id = 1,
                         location_id = 1,
                         age_group_id = 22,
                         sex_id = sex_id,
                         release_id = 7,
                         source = source,
                         version_id = version_id) %>%
    .[, c(drop_cols) := NULL] %>%
    melt(., id.vars = "cause_id") %>%
    setnames(., "value", "intr_value")
  
  # Subarachnoid hemorrhage
  sub_mort <- get_draws(gbd_id_type = "cause_id",
                        gbd_id = 497,
                        year_id = 2019,
                        measure_id = meas_id,
                        metric_id = 1,
                        location_id = 1,
                        age_group_id = 22,
                        sex_id = sex_id,
                        source = source,
                        release_id = 7,
                        version_id = version_id) %>%
    .[, c(drop_cols) := NULL] %>%
    melt(., id.vars = "cause_id") %>%
    setnames(., "value", "sub_value")
  
  # Merge to get mortality for hemorrhagic stroke
  hem_mort <- merge(intr_mort, sub_mort, by = "variable") %>%
    .[, value := intr_value + sub_value] %>%
    .[, c("sub_value", "intr_value", "cause_id.x", "cause_id.y") := NULL] %>%
    .[, cause_id := 494]
  
}else{
  
  hem_mort <- data.table()
}
# Append other mortality estimates
mort <- rbind(mort, hem_mort) %>% 
  setnames(., "value", "deaths")

mort <- mort[cause_id %in% causes$cause_id]

# mean of mortality draws
mort[, avg_deaths := mean(deaths), by = "cause_id"]
mort <- unique(mort[,.(cause_id, avg_deaths)])
mort[, weight := avg_deaths/sum(avg_deaths)]


#----------------------------------------------------------
# Step 3: Pull data and risk curves for each outcome
#----------------------------------------------------------

# calculate exposure range for the risk
rr_data <- rbindlist(lapply(causes$pair, function(p){
  
  input_data <- readRDS(paste0(mrbrt_dir, "/00_prepped_data/", p, ".RDS"))
  df <- as.data.table(input_data$df)
  df[, b_midpoint := b_0 + (b_1-b_0)/2]
  df[, a_midpoint := a_0 + (a_1-a_0)/2]
  
  df[, pair := p]
  return(df)
  
}), fill = T)


rr_draws <- rbindlist(lapply(causes$pair, function(p){
  
  is_j_shaped <- F 
  
  signal_model_path <- paste0(central_model_output_folder, "/", p, "/signal_model.pkl")
  linear_model_path <- paste0(central_model_output_folder, "/", p, "/new_linear_model.pkl")
  signal_model <- py_load_object(filename = signal_model_path, pickle = "dill")
  linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")
  
  # get_draws
  y_draws <- get_ln_rr_draws(signal_model,
                              linear_model,
                              risk = seq(0, 200, length.out = 1000),  #use 200 as upper bound instead of 238 (data max)
                              num_draws = 1000L,
                              normalize_to_tmrel = is_j_shaped,
                              fe_only = !random_effects)
  y_draws <- as.data.table(y_draws)
  
  setnames(y_draws, colnames(y_draws), c("exposure", paste0("draw_",0:999)))
  
  df_final_pred <- summarize_draws(y_draws)
  
  y_draws[, ro := p]
  y_draws[, cause_id := causes[pair == p, cause_id]]

  return(y_draws)
  
}), use.names = TRUE)


#----------------------------------------------------------
# Step 4: Calculate all-cause mortality curve
#----------------------------------------------------------

# Melt RR draws long
rr_long <- melt(rr_draws, id.vars = c("ro","exposure", "cause_id"))
rr_long$variable <- as.character(rr_long$variable)

# Merge on deaths
rr_long <- merge(rr_long, mort, by = c("cause_id"))

rr_long[, rr := sum(value*weight), by = c("exposure", "variable")] # Take weighted mean of RRs

# save separate dataset for plotting purposes
rr_long[, `:=` (mean_rr = mean(rr), mean_value = mean(value)), by = c("cause_id", "exposure")]
out_rr <- unique(rr_long[,.(cause_id, ro, exposure, weight, mean_value)])
all_rr <- unique(rr_long[,.(exposure, mean_rr)])
all_rr[, `:=` (ro = "all-cause", weight = 1)]
setnames(out_rr, "mean_value", "mean_rr")
means_for_plot <- rbind(all_rr, out_rr, fill = T)
means_for_plot[, ro:=gsub("_fixed", "", ro)]

rr_long_adj <- unique(rr_long[,.(variable, exposure, rr)])
write.csv(rr_long_adj, paste0(save_dir, "/",risk_grep,"_allcausedraws_",meas_name,note,ifelse(random_effects, "_re", "_fe"),".csv"), row.names = F)


# what exposure value minimizes each draw ?
rr_long_adj[, min_rr := min(rr), by = "variable"]

min_exposure_details <- rr_long_adj[rr==min_rr, .(variable, exposure)]
setnames(min_exposure_details, "exposure", "min_exp")
min_exposure_vec <- min_exposure_details$min_exp


# summarize and move to normal space
rr_long_adj[, all_cause_rr := exp(mean(rr)), by = "exposure"] # Median behaves much better than mean in draws with random effects
rr_long_adj[, all_cause_rr_lo := exp(quantile(rr, probs = 0.025)), by = "exposure"]
rr_long_adj[, all_cause_rr_hi := exp(quantile(rr, probs = 0.975)), by = "exposure"]

all_cause_rr <- unique(rr_long_adj[,.(exposure, all_cause_rr, all_cause_rr_lo, all_cause_rr_hi)])

# now scale each draw to it's minimum
rr_long_adj[, rr_at_min := min(rr), by = "variable"]
rr_long_adj[, rr_shifted := rr - rr_at_min]

#----------------------------------------------------------
# Step 5: Plot, calculate TMREL, etc
#----------------------------------------------------------
tmrel <- mean(min_exposure_vec)
tmrel_lower <- quantile(min_exposure_vec, probs = 0.025) 
tmrel_upper <- quantile(min_exposure_vec, probs = 0.975)

# Density plot of TMREL draws
density <- ggplot(data.table(tmrel = min_exposure_vec), aes(x = tmrel)) + 
  geom_histogram() + 
  theme_bw() 
print(density)

# Plot of all-cause risk curve
plot <- ggplot(all_cause_rr, aes(x = exposure, y = all_cause_rr))+geom_line()+
  geom_ribbon(aes(ymin = all_cause_rr_lo, ymax = all_cause_rr_hi),alpha =0.4)+
  theme_bw()+
  labs(y = "Relative Risk", x = "Exposure", title = "All cause")
print(plot)


# Plot of mean all-cause risk curve components
plot <- ggplot(means_for_plot, aes(x = exposure, y = mean_rr*weight, color = ro))+geom_line()+
  theme_bw()+
  labs(y = "Relative Risk * Weight", x = "Exposure", title = "All cause as additive of components")
print(plot)
