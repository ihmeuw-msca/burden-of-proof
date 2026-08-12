## SET-UP #################################################################
rm(list = ls())

source("filepath")
library(ggpubr)
library(data.table)
library(ggplot2)
library(officer)
library(flextable)
library(purrr)
library(stringr)
library(gridExtra)

## FILEPATHS ######################################################################################
results_fp <- "filepath"
viz_fp <- "filepath"
tables_fp <- "filepath"
reference_fp <- "filepath"

# FIGURING OUT WHAT FOLDERS WE ARE PULLING FROM
options <- c("notrim", "trim")
alts <- c("outcome_sens", "perp_specific", "anyperp_specific", "nopregnancyrecall", "onlypregnancy", "maleonly", "femaleonly", "noadjustment")
full <- c()
for(i in options){
  temp <- paste0(alts, "_", i)
  full <- c(temp, full)
}
options <- c("main_models", "notrim", full)
age_of_interest <- "childhood"

# Setting parameters
for_paper <- T

## SPECIAL FUNCTIONS ##############################################################################
generate_stars <- function(n) {
  if(is.na(n)){
    return(paste0(""))
  } else {
    return(paste(rep("\U2605", n), collapse = ""))
  }
}
str_replace_nth <- function(x, pattern, replacement, n) {
  g <- gregexpr(pattern, x)[[1]][n]
  s <- scan(text = gsub("[()]", "", pattern),
            sep = "|",
            what = "")
  substr(x, g, g) <- replacement[match(substr(x, g, g), s)]
  x
}
aesth <- theme_bw() + theme(axis.title.x = element_text(size=7,face="plain", family='Sans'),
                            axis.title.y = element_text(size=7,face="plain", family='Sans'),
                            axis.text.x = element_text(size=7,face="plain", family='Sans'),
                            axis.text.y = element_text(size=7,face="plain", family='Sans'),
                            strip.text = element_text(size=7,face="bold"),
                            plot.background = element_blank(),
                            legend.text=element_text(size=7, family='Sans'),
                            legend.title=element_text(size=7, face='bold', family='Sans'),
                            plot.title=element_text(size=7, face='bold', family='Sans'),
                            plot.caption=element_text(size=7, family='Sans'))
# Order of causes
cause_order <- rev(c("Maternal abortion and miscarriage", "Maternal hemorrhage","Sexually transmitted infections excluding HIV", "HIV/AIDS", 
                     "Anxiety disorders", "Major depressive disorder", "Bipolar disorder","Alcohol use disorders", "Drug use disorders",
                     "Schizophrenia", "Eating disorders", "Anorexia nervosa", "Bulimia nervosa", "Conduct disorder", 
                     "Diabetes mellitus type 2", "Asthma", "Ischemic heart disease", "Stroke","Gynecological diseases", "Migraine",
                     "Self-harm"))

## PULLING IN EXTRA FILES ##########################################################################
# This file has the legible names for the outcomes that we want to use
outcome_maps <- fread(paste0(reference_fp, "outcome_name_mapping.csv"))
# This file has the legible descriptions of the different outcome-specific sensitivity analyses
sens_outcome_maps <- fread(paste0(reference_fp, "outcome_sens_mapping.csv"))
# This file has Cory's final results that will be used for comparison
corys_results <- fread(paste0(reference_fp, "corys_csa_ipv_results.csv")) %>% .[!is.na(mean_rr)]
# Extra cause metadata
cause <- get_cause_metadata(cause_set_id = 2, release_id = 16)

## PULLIN IN RESULTS ###############################################################################
for(versions in options){
  temporary <- data.table()
  temp_data <- data.table()
  if(versions == "main_models"){
    fp_of_interest <- paste0(results_fp)
  } else {
    fp_of_interest <- paste0(results_fp, "sens_analyses/")
  }
  filepaths <- list.dirs(paste0(fp_of_interest, versions, "/"), full.names = F)[-1]
  
  filepaths <- filepaths[!(filepaths %like% "archive")]
  if(for_paper){
    filepaths <- filepaths[!(filepaths %like% "adulthood_ipv")]
    filepaths <- filepaths[!(filepaths %like% "childhood_sexual")]
  }

  for(i in seq_along(filepaths)){
    data <- fread(paste0(fp_of_interest, versions, "/", filepaths[i], "/", filepaths[i], ".csv"))
    inner_quantiles <- fread(paste0(fp_of_interest, versions, "/", filepaths[i], "/inner_quantiles.csv")) # RR without gamma
    outer_quantiles <- fread(paste0(fp_of_interest, versions, "/", filepaths[i], "/outer_quantiles.csv")) # RR with gamma
    selected_covs <- yaml::read_yaml(paste0(fp_of_interest, versions, "/", filepaths[i], '/cov_finder_result.yaml'))$selected_covs # Selected covariates
    model_out <- yaml::read_yaml(paste0(fp_of_interest, versions, "/", filepaths[i], '/summary.yaml')) # Model results
    
    if(nrow(data) <= 10 & versions %like% "_trim"){
      next
    } else if((nrow(data) > 10 & (versions %like% "_notrim" & versions != "notrim"))){
      next
    } else if(versions %like% "_notrim" | (nrow(data) <= 10 & versions == "main_models")){
      trimming <- "0%"
    } else if(versions %like% "_trim" | (versions == "main_models" & nrow(data) > 10)){
      trimming <- "10%"
    }
    
    mean_rr <- sprintf("%.2f", exp(inner_quantiles$V3[2]))
    rr_without_gamma <- paste0(sprintf("%.2f", exp(inner_quantiles$V3[2])), " (", sprintf("%.2f", exp(inner_quantiles$V1[2])), "-", sprintf("%.2f", exp(inner_quantiles$V5[2])), ")")
    ui_without_gamma <- paste0(sprintf("%.2f", exp(inner_quantiles$V1[2])), "-", sprintf("%.2f", exp(inner_quantiles$V5[2])))
    
    rr_with_gamma <- paste0(sprintf("%.2f", exp(outer_quantiles$V3[2])), " (", sprintf("%.2f", exp(outer_quantiles$V1[2])), "-", sprintf("%.2f", exp(outer_quantiles$V5[2])), ")")
    ui_with_gamma <- paste0(sprintf("%.2f", exp(outer_quantiles$V1[2])), "-", sprintf("%.2f", exp(outer_quantiles$V5[2])))
    
    upper_with_gamma <- round(exp(outer_quantiles$V5[2]), 2)
    lower_with_gamma <- round(exp(outer_quantiles$V1[2]), 2)
    upper_without_gamma <- round(exp(inner_quantiles$V5[2]), 2)
    lower_without_gamma <- round(exp(inner_quantiles$V1[2]), 2)
    
    bprf <- sprintf("%.2f", exp(outer_quantiles$V2[2]))
    # risk-outcome score
    alt_round <- function(x) {
      max(abs(round(x, 2)), abs(signif(x, 2))) * sign(x)
    }
    gamma <- paste0(signif(model_out$gamma[1], 3), " (", signif(model_out$gamma[2], 3), ")")
    publication_bias <- ifelse(model_out$pub_bias == 1, "Yes", "No")
    cause_id <- model_out$cause_id
    ros <- alt_round(model_out$score)
    # ros <- round(model_out$score, 2)
    selected_covariates <- ifelse(length(selected_covs) > 1, paste0(selected_covs, collapse = ", "), 
                                  ifelse(length(selected_covs == 1), selected_covs, "None"))
    num_studies <- paste0(length(unique(data$study_id)), " (", nrow(data), ")")
    thresholds <- c(log(1 + c(0.15, 0.5, 0.85)), Inf)
    stars <- c("1", "2", "3", "4", "5")
    star_rating <- cut(model_out$score, c(-Inf, 0, thresholds), right = FALSE,
                       labels = stars, include.lowest = TRUE)
    
    temp <- as.data.table(
      cbind(
        age_of_exposure = sub("_(.*)", "", filepaths[i]),
        violence_type = sub(".*?_(.*?)-.*", "\\1", filepaths[i]),
        cause_name = sub(".*-", "", filepaths[i]),
        mean_rr, 
        ui_without_gamma, 
        ui_with_gamma,
        upper_with_gamma, 
        lower_with_gamma, 
        upper_without_gamma, 
        lower_without_gamma,
        rr_without_gamma,
        rr_with_gamma,
        bprf,
        ros, 
        selected_covariates, 
        num_studies,
        star_rating,
        gamma, 
        publication_bias, 
        trimming)
    )
    temporary <- rbindlist(list(temp, temporary), fill = T)

    data$model_name <- filepaths[i]
    data$cause_name <- sub(".*-", "", filepaths[i])
    data$violence_type <-  sub(".*?_(.*?)-.*", "\\1", filepaths[i])
    data$age_of_exposure <- sub("_(.*)", "", filepaths[i])
    
    data[, study_label := stringr::str_to_title(Study_ID)]
    
    temp_data <- rbindlist(list(data, temp_data), fill = T)
  }
  if(nrow(temp_data) > 0 & !(versions %like% "outcome_sens")){
    temporary <- merge(temporary, outcome_maps, by = "cause_name", all.x = T)
    temporary[, study_version := "Present analysis"]
    temporary[, violence_labels := ifelse(violence_type == "neglect", "Neglect", paste0(stringr::str_to_title(violence_type), "\nabuse"))]
    temp_data <- merge(temp_data, temporary, by = c("violence_type", "cause_name", "age_of_exposure"), all = T)
    
    assign(paste0(versions, "_results"), temporary)
    assign(paste0(versions, "_data"), temp_data)
  } else if(nrow(temp_data) > 0) {
    temporary[, study_version := "Present analysis"]
    temporary[, violence_labels := ifelse(violence_type == "neglect", "Neglect", paste0(stringr::str_to_title(violence_type), "\nabuse"))]
    temp_data <- merge(temp_data, temporary, by = c("violence_type", "cause_name", "age_of_exposure"), all = T)
    temporary <- merge(temporary, sens_outcome_maps, by = "cause_name", all.x = T)
    temp_data <- merge(temp_data, sens_outcome_maps, by = "cause_name", all.x = T)
    
    assign(paste0(versions, "_results"), temporary)
    assign(paste0(versions, "_data"), temp_data)
  }
}

if(for_paper){
  corys_results[, `:=` (study_version = "Previously published\nanalysis", violence_labels = ifelse(violence_type == "ipv", "Intimate partner\nsexual and/or\nphysical violence", "Sexual\nabuse"))]
  corys_results <- merge(corys_results, cause[, .(cause_id, cause_name)], by = "cause_id")
  setnames(corys_results, "cause_name", "full_cause_names")
  corys_results[, star_rating := as.numeric(star_rating)]
  corys_results[, ros := ifelse(ros %like% "N/A", NA, ros)]
  corys_results[, ros := as.numeric(ros)]
  corys_results[, star_rating := ifelse(star_rating == 0, NA, star_rating)]
  
  main_models_results <- rbindlist(list(main_models_results, corys_results), fill = T)
}

for(versions in c("main_models", "notrim", alts)){
  if(versions == "main_models" | versions == "notrim"){
    temporary <- get(paste0(versions, "_results"))
    temp_data <- get(paste0(versions, "_data"))
  } else {
    if(exists(paste0(versions, "_trim_results"))){
      temporary_trim <- get(paste0(versions, "_trim_results"))
      temp_data_trim <- get(paste0(versions, "_trim_data"))
      if(!exists(paste0(versions, "_notrim_results"))){
        temporary <- copy(temporary_trim)
        temp_data <- copy(temp_data_trim)
      }
    }
    if(exists(paste0(versions, "_notrim_results"))){
      temporary_notrim <- get(paste0(versions, "_notrim_results"))
      temp_data_notrim <- get(paste0(versions, "_notrim_data"))
      
      if(!exists(paste0(versions, "_trim_results"))){
        temporary <- copy(temporary_notrim)
        temp_data <- copy(temp_data_notrim)
        
      }
    }
    if(exists(paste0(versions, "_notrim_results")) & exists(paste0(versions, "_trim_results"))){
      temporary <- rbindlist(list(temporary_notrim, temporary_trim))
      temp_data <- rbindlist(list(temp_data_notrim, temp_data_trim), fill = T)
    }
    
  }

  if(age_of_interest == "adulthood"){
    temporary <- temporary[age_of_exposure == "adulthood"]
    temporary[, violence_labels := gsub("abuse", "violence", violence_labels)]
    temporary[, violence_labels := ifelse(violence_labels == "Ipv violence", "Intimate partner violence", violence_labels)]
    temporary[, violence_labels := ifelse(violence_labels %like% "Intimate partner", "Intimate partner\nsexual and/or\nphysical violence*", violence_labels)]
    temporary[, violence_labels := factor(violence_labels, levels = c("Neglect", "Psychological\nviolence", "Physical\nviolence", "Sexual\nviolence", "Intimate partner\nsexual and/or\nphysical violence*"))]
    
    temp_data <- temp_data[age_of_exposure == "adulthood"]
    temp_data[, violence_labels := gsub("abuse", "violence", violence_labels)]
    temp_data[, violence_labels := ifelse(violence_labels %like% "Intimate partner", "Intimate partner\nsexual and/or\nphysical violence*", violence_labels)]
    temp_data[, violence_labels := factor(violence_labels, levels = c("Neglect", "Psychological\nviolence", "Physical\nviolence", "Sexual\nviolence", "Intimate partner\nsexual and/or\nphysical violence*"))]
    
    if(for_paper){
      temporary <- temporary[!violence_type == "ipv"]
    }
  } else if(age_of_interest == "childhood"){
    temporary <- temporary[age_of_exposure == "childhood"]
    temporary[, violence_labels := factor(violence_labels, levels = c("Physical\nabuse", "Psychological\nabuse", "Neglect", "Sexual\nabuse", "Intimate partner\nsexual and/or\nphysical violence"))]
    
    temp_data <- temp_data[age_of_exposure == "childhood"]
    temp_data[, violence_labels := factor(violence_labels, levels = c("Physical\nabuse", "Psychological\nabuse", "Neglect", "Sexual\nabuse", "Intimate partner\nsexual and/or\nphysical violence"))]
    
  } else {
    temporary[age_of_exposure == "adulthood", violence_labels := gsub("abuse", "GBV", violence_labels)]
    temporary[, violence_labels := ifelse(violence_labels == "Ipv\nGBV", "Intimate partner\nsexual and/or\nphysical violence", violence_labels)]
    
    temporary[, violence_labels := factor(violence_labels, levels = c("Neglect", "Psychological\nabuse", "Physical\nabuse", "Sexual\nabuse", "Psychological\nGBV", "Physical\nGBV", "Sexual\nGBV", "Intimate partner\nsexual and/or\nphysical violence"))]
  }
  
  
    temporary[, star_rating := as.numeric(star_rating)]
    temporary$stars <- sapply(temporary$star_rating, generate_stars)
    
    temporary[, mean_rr_label := ifelse(star_rating == 0 | is.na(star_rating), "", mean_rr)]
    temporary[, full_cause_names := factor(full_cause_names, levels = cause_order)]
    temp_data[, study_label := ifelse(study_label %in% c("Llosamartínez 2019", "Llosamartinez 2019"), "Llosa Martínez 2019", 
                                      ifelse(study_label == "Xavierhall 2021", "Xavier Hall 2021", 
                                             ifelse(study_label == "Vanroode 2009", "van Roode 2009", 
                                                    ifelse(study_label == "Kiselydmedres 2022", "Kisely 2022", 
                                                           ifelse(study_label == "Tenhave 2019", "Ten Have 2019", study_label)))))]
    if(for_paper == F){
      temporary[, age_of_exposure := ifelse(age_of_exposure == "childhood", "Violence Against Children", ifelse(age_of_exposure == "adulthood", "Adulthood Gender-Based Violence", "Sexual Violence Regardless of Age"))]
      temporary[, age_of_exposure := factor(age_of_exposure, levels = c("Violence Against Children", "Adulthood Gender-Based Violence", "Sexual Violence Regardless of Age"))]
      temporary <- temporary[star_rating != 0]
    }
    
    if(versions == "outcome_sens"){
      setorder(temporary, full_cause_names)
      setnames(temporary, c("full_cause_names", "sens_analysis"), c("cause_names_full", "full_cause_names"))
      setnames(temp_data, c("full_cause_names", "sens_analysis"), c("cause_names_full", "full_cause_names"))
      temporary[, full_cause_names := str_replace_nth(full_cause_names, " ", "\n", round(str_count(full_cause_names, " ")/2,0))]
      temp_data$full_cause_names <- NULL
      temp_data <- merge(temp_data, unique(temporary[, .(cause_name, full_cause_names, violence_type, age_of_exposure)]), by = c("violence_type", "cause_name", "age_of_exposure"), all.x = T)
    }
    
    fill_options <- c("star_rating")
    for(filling_value in fill_options){
      temporary[, star_rating := ifelse(star_rating == 0, NA, star_rating)]
      for(i in unique(temporary$age_of_exposure)){
        if(filling_value == "ros"){
          fill_label <- "Risk Outcome Score"
        } else if(filling_value == "star_rating"){
          fill_label <- "Star rating"
          temporary[, star_rating := ifelse(is.na(star_rating), 0, star_rating)]
        } else if(filling_value == "mean_rr"){
          fill_label <- "Mean Relative Risk"
        }
        if(for_paper){
          p <- ggplot(temporary[age_of_exposure == i], aes(x = violence_labels, y = full_cause_names, fill = as.factor(get(filling_value)))) + 
            geom_tile()
          if(filling_value == "ros"){
            p <- p + geom_tile(data = temporary[age_of_exposure == i & is.na(star_rating)], aes(x = violence_labels, y = full_cause_names), fill = "black") +
              geom_text(aes(label = paste0(mean_rr)), size=6, colour = "white", vjust = -0.75) + 
              geom_text(aes(label = paste0(stars)), size=5, colour = "white", vjust = 0.75)
          } else if(filling_value == "star_rating"){
            p <- p + geom_text(aes(label = paste0(mean_rr_label)), size=1.5, colour = "white")
          } else if(filling_value == "mean_rr"){
            p <- p + 
              geom_text(aes(label = paste0(mean_rr, " (", ui_with_gamma, ")")), size=6, colour = "white", vjust = -0.75) + 
              geom_text(aes(label = paste0(stars)), size = 10, colour = "white", vjust = 0.75)
          }
          if(length(unique(temporary[age_of_exposure == i]$star_rating)) <= 4 & min(unique(as.numeric(temporary[age_of_exposure == i]$star_rating))) == 0){
            p <- p +
              scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(1,2,4,6,8)]),
                                guide = guide_legend(reverse = TRUE))
          } else if(length(unique(temporary[age_of_exposure == i]$star_rating)) > 4){
            p <- p +
              scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(2,4,6,8)]),
                                guide = guide_legend(reverse = TRUE))
          } else if(length(unique(temporary[age_of_exposure == i]$star_rating)) <= 4 & min(unique(as.numeric(temporary[age_of_exposure == i]$star_rating))) > 0){
            p <- p +
              scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(1,2,3,4,6,8)]),
                                guide = guide_legend(reverse = TRUE)) 
          }
          p <- p + 
            scale_y_discrete(expand = c(0,0)) + 
            theme_bw() +
            theme(axis.text.x = element_text(face='bold'), 
                  plot.title = element_text(face="bold"), 
                  panel.grid.major = element_blank(), 
                  axis.text.y=element_text(face = 'bold'), 
                  text = element_text(size = 5), 
                  legend.title = element_text(face = "bold"),
                  strip.text.x = element_blank(), 
                  strip.text.y = element_text(size = 5, face = "bold")) +
            labs(title = paste0("Associations between ", i, " violence\nexposure and various health outcomes"), 
                 y = "", x = "", fill = fill_label) +
            facet_grid(cols = vars(study_version), scales = "free_x", space = "free_x") 
        } else {
          p <- ggplot(temporary, aes(x = violence_labels, y = full_cause_names, fill = as.factor(get(filling_value)))) + 
            geom_tile()
          p <- p + geom_text(aes(label = paste0(mean_rr_label)), size=1.5, colour = "white")
          p <- p +
            scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(1,2,4,6,8)]),
                              guide = guide_legend(reverse = TRUE)) +
            scale_y_discrete(position = "right", expand = c(0,0)) + 
            labs(title = paste0("Associations between GBV or VAC exposure and various health outcomes"), 
                 y = "", x = "", fill = "Strength of Evidence") +
            facet_grid(cols = vars(age_of_exposure), scales = "free_x", space = "free_x") +
            theme_bw() +
            theme(axis.text.x = element_text(face='bold'), 
                  plot.title = element_text(face="bold"), 
                  panel.grid.major = element_blank(), 
                  axis.text.y=element_text(face = 'bold'), 
                  text = element_text(size = 5), 
                  legend.title = element_text(face = "bold"),
                  legend.position = "bottom",
                  strip.text.x = element_text(size = 5, face = "bold"),
                  strip.text.y = element_text(size = 6, face = "bold"), 
                  strip.placement = "outside", 
                  strip.background.x = element_rect(fill = "white", color = "white"))
          
        }

          p <- p + 
          geom_hline(yintercept = 0.5 + 0:35, colour = "white", linewidth = 1) +
          geom_vline(xintercept = 0.5 + 0:5, colour = "white", linewidth = 1) +
          scale_x_discrete(expand = c(0,0), position = "top")
        
        assign(paste0(versions, "_", i, "_heatmap_update_", filling_value), p)
      }
    }
    
     grDevices::cairo_pdf(paste0(viz_fp, versions, '_heatmap_new.pdf'), width = 4.72441, height = 3.5433, onefile = T)     
     for(filling_value in fill_options){
       for(i in unique(temporary$age_of_exposure)){
         print(get(paste0(versions, "_", i, "_heatmap_update_", filling_value)))
       }
    }
    dev.off()

    options(scipen=999)
    for(i in unique(temp_data$age_of_exposure)){
      plot_data <- copy(temp_data[age_of_exposure == i])
      plot_data[, sex := ifelse(is.na(sex) | sex == "", "Combined Male and Female", sex)]
      plot_data[, study_label := factor(study_label, levels = rev(sort(unique(plot_data$study_label))))]
      plot_data[, ymax_val := uniqueN(study_label), by = c("full_cause_names", "violence_labels")]

      plot_data[, scale_to_using := 1:.N, by = c("full_cause_names", "violence_labels")]
      plot_data[scale_to_using != 1, `:=` (lower_with_gamma = NA, upper_with_gamma = NA, lower_without_gamma = NA, upper_without_gamma = NA)]
      plot_data[, cause_labels := ifelse(full_cause_names == "Sexually transmitted infections excluding HIV", "STIs, excluding HIV", 
                                         ifelse(full_cause_names == "Maternal abortion and miscarriage", "Maternal abortion\nand miscarriage", full_cause_names))]
      plot_data[, number_of_studies := length(unique(study_label)), by = c("full_cause_names", "violence_labels")]
      plot_data[, number_of_obs := .N, by = c("full_cause_names", "violence_labels")]
      plot_data[, cause_labels := paste0(cause_labels, "\n(N = ", number_of_studies, ")")]
      plot_data[, different_columns := 1:.N, by = c("violence_labels")]
      plot_data[, different_columns := max(different_columns), by = c("full_cause_names", "violence_labels")]
      plot_data[, column_cutoff := ceiling(max(different_columns)/2), by = c("violence_labels")]
      plot_data[, different_columns := factor(ifelse(different_columns >= column_cutoff, "First column", "Second column"), levels = c("First column", "Second column"))]
      
      for(m in unique(plot_data$violence_labels)){
        saving_label <- str_to_lower(gsub("\n", " ", m))
        x_minimum <- min(c(exp(plot_data[violence_labels == m]$ln_rr-1.96*plot_data[violence_labels == m]$ln_rr_se), as.numeric(temporary[age_of_exposure == i & violence_labels == m]$lower_with_gamma))) - 0.001
        x_maximum <- max(c(exp(plot_data[violence_labels == m]$ln_rr+1.96*plot_data[violence_labels == m]$ln_rr_se), as.numeric(temporary[age_of_exposure == i & violence_labels == m]$upper_with_gamma))) + 1
        
        x_min <- min(c(exp(plot_data[violence_labels == m]$ln_rr), as.numeric(temporary[age_of_exposure == i & violence_labels == m]$lower_with_gamma))) - 0.001
        x_max <- max(c(exp(plot_data[violence_labels == m]$ln_rr), as.numeric(temporary[age_of_exposure == i & violence_labels == m]$upper_with_gamma))) + 1
        
        print(paste0(i, " ", m, " ", x_minimum, " ", x_maximum))
        temp <- ggplot(data=plot_data[violence_labels == m & different_columns == "First column"])+
          geom_vline(xintercept=1, linetype='dashed')
        if(age_of_interest == "adulthood"){
          temp <- temp +
            geom_pointrange(aes(x=exp(ln_rr), xmin=exp(ln_rr-1.96*ln_rr_se), xmax=exp(ln_rr+1.96*ln_rr_se), 
                                y=study_label,
                                shape=sex,
                                color=as.factor(is_outlier)),
                            position=position_dodge2(0.75), size = 0.2, linewidth = 0.2) +
            geom_rect(aes(xmin = as.numeric(lower_with_gamma), xmax = as.numeric(upper_with_gamma),
                          ymin=0, ymax=ymax_val+1),
                      fill='#6DBCC3', alpha = 0.3, inherit.aes = FALSE)+
            geom_rect(aes(xmin = as.numeric(lower_without_gamma), xmax = as.numeric(upper_without_gamma),
                          ymin=0, ymax=ymax_val+1),
                      fill='#6DBCC3', alpha = 0.5, inherit.aes = FALSE)+
            geom_vline(aes(xintercept=as.numeric(bprf)), linetype='solid', color="#FF0000")+
            geom_vline(aes(xintercept=as.numeric(mean_rr)), linetype='solid', color="#6DBCC3")+
            scale_color_manual(values=c('1'='red', '0'='black'))+
            scale_shape_manual(values=c('Male'=4, 'Female'=17, "Combined Male and Female" = 16))+
            labs(title= '', 
                 shape='Sex',
                 color = 'Is Outlier',
                 y='Study ID', 
                 x=''
            )+
            theme_bw()+
            theme(plot.title=element_text(face='bold'), 
                  strip.text = element_text(face = "bold"), 
                  strip.text.x = element_blank(),
                  strip.background = element_rect(fill = "white", color = "black", size = 1, linetype = "solid"), 
                  strip.text.y.right = element_text(angle = 0), 
                  text = element_text(size = 5),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
            guides(linetype = "none", shape = "none", color = "none")+
            scale_x_log10(limits = c(x_minimum, x_maximum), labels = scales::number_format(accuracy = 0.1)) + 
            coord_cartesian(xlim = c(x_min, x_max))  +
            facet_grid(rows = vars(cause_labels), scales = "free_y", space = "free_y")
        } else {
          temp <- temp +
            geom_pointrange(aes(x=exp(ln_rr), xmin=exp(ln_rr-1.96*ln_rr_se), xmax=exp(ln_rr+1.96*ln_rr_se), 
                                y=study_label,
                                shape=as.factor(is_outlier),
                                color=as.factor(is_outlier)),
                            position=position_dodge2(0.75), size = 0.15, linewidth = 0.15) +
            geom_rect(aes(xmin = as.numeric(lower_with_gamma), xmax = as.numeric(upper_with_gamma),
                          ymin=0, ymax=ymax_val+1),
                      fill='#6DBCC3', alpha = 0.3, inherit.aes = FALSE)+
            geom_rect(aes(xmin = as.numeric(lower_without_gamma), xmax = as.numeric(upper_without_gamma),
                          ymin=0, ymax=ymax_val+1),
                      fill='#6DBCC3', alpha = 0.5, inherit.aes = FALSE)+
            geom_vline(aes(xintercept=as.numeric(bprf)), linetype='solid', color="#FF0000")+
            geom_vline(aes(xintercept=as.numeric(mean_rr)), linetype='solid', color="#6DBCC3")+
            scale_color_manual(values=c('1'='red', '0'='black'))+
            scale_shape_manual(values=c('1'=4, '0'=16))+
            labs(title= '', 
                 shape='Is Outlier',
                 color = 'Is Outlier',
                 y='Study ID', 
                 x=''
            )+
            theme_bw()+
            theme(plot.title=element_text(face='bold'), 
                  strip.text = element_text(face = "bold"), 
                  strip.text.x = element_blank(),
                  strip.background = element_rect(fill = "white", color = "black", size = 1, linetype = "solid"), 
                  strip.text.y.right = element_text(angle = 0), 
                  text = element_text(size = 5),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
            guides(linetype = "none", shape = "none", color = "none")+
            scale_x_log10(limits = c(x_minimum, x_maximum), labels = scales::number_format(accuracy = 0.1)) + 
            coord_cartesian(xlim = c(0.09, 20)) +
            facet_grid(rows = vars(cause_labels), scales = "free_y", space = "free_y")
        }
     
        
        if(nrow(plot_data[violence_labels == m & different_columns == "Second column"]) > 0){
          temp2 <- ggplot(data=plot_data[violence_labels == m & different_columns == "Second column"])+
            geom_vline(linetype = "dashed", xintercept=1)
          if(age_of_interest == "adulthood"){
            temp2 <- temp2 +
            geom_pointrange(aes(x=exp(ln_rr), xmin=exp(ln_rr-1.96*ln_rr_se), xmax=exp(ln_rr+1.96*ln_rr_se), 
                                y=study_label,
                                shape=sex,
                                color=as.factor(is_outlier)),
                            position=position_dodge2(0.75), linewidth = 0.2, size = 0.2) +
            geom_rect(aes(xmin = as.numeric(lower_with_gamma), xmax = as.numeric(upper_with_gamma),
                          ymin=0, ymax=ymax_val+1, alpha = "95% UI\nwith gamma"),
                      fill='#6DBCC3')+
            geom_rect(aes(xmin = as.numeric(lower_without_gamma), xmax = as.numeric(upper_without_gamma),
                          ymin=0, ymax=ymax_val+1, alpha = "95% UI\nwithout gamma"),
                      fill='#6DBCC3')+
            geom_vline(aes(xintercept=as.numeric(bprf), linetype = "BPRF"), color="#FF0000")+
            geom_vline(aes(xintercept=as.numeric(mean_rr), linetype = "Estimated mean RR"), linetype='solid', color="#6DBCC3")+
            scale_color_manual(limits = c('1', '0'), values=c('1'='red', '0'='black'), labels = c('1' = "Trimmed", '0' = "Included"))+
            scale_linetype_manual(limits = c("Estimated mean RR", "BPRF", "Null RR"), values = c("BPRF" = "solid", "Estimated mean RR" = "solid", "Null RR" = "dashed")) +
            scale_shape_manual(limits = c('Male', 'Female', 'Combined Male and Female'), values=c('Male'=4, 'Female'=17, "Combined Male and Female" = 16), labels = c('Male only', 'Female only', 'Combined Male\nand Female'))+
            scale_alpha_manual(values = c("95% UI\nwithout gamma" = 0.5, "95% UI\nwith gamma" = 0.3)) +
            labs(title= '', 
                 alpha = "Model Uncertainty",
                 shape='Sex of Sample',
                 color = 'Effect Sizes',
                 linetype = "Model Results",
                 y='', 
                 x='')+
            theme_bw()+
            theme(plot.title=element_text(face='bold'), 
              strip.text = element_text(face = "bold"), 
              strip.text.x = element_blank(),
              strip.background = element_rect(fill = "white", color = "black", size = 1, linetype = "solid"), 
              strip.text.y.right = element_text(angle = 0), 
              # text = element_text(size = 15), 
              legend.title = element_text(size = 7, face = "bold"), 
              legend.text = element_text(size = 5),
              text = element_text(size = 5),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
            # scale_x_continuous(trans = "log10", limits = c(x_minimum, x_maximum)) + 
            scale_x_log10(limits = c(x_minimum, x_maximum), labels = scales::number_format(accuracy = 0.1)) +
            facet_grid(rows = vars(cause_labels), scales = "free_y", space = "free_y") +
              coord_cartesian(xlim = c(x_min, x_max)) +
            guides(shape = guide_legend(order = 2), color = guide_legend(order = 1),
                    alpha = guide_legend(override.aes = list(fill = c("#6DBCC3", "#6DBCC3")), order = 3), 
                   linetype = guide_legend(override.aes = list(color = c("#6DBCC3", "#FF0000", "#000000")), order = 4)
                   )
          } else {
            temp2 <- temp2 +
              geom_pointrange(aes(x=exp(ln_rr), xmin=exp(ln_rr-1.96*ln_rr_se), xmax=exp(ln_rr+1.96*ln_rr_se), 
                                  y=study_label,
                                  shape=as.factor(is_outlier),
                                  color=as.factor(is_outlier)),
                              position=position_dodge2(0.75), linewidth = 0.15, size = 0.15) +
              geom_rect(aes(xmin = as.numeric(lower_with_gamma), xmax = as.numeric(upper_with_gamma),
                            ymin=0, ymax=ymax_val+1, alpha = "95% UI\nwith gamma"),
                        fill='#6DBCC3')+
              geom_rect(aes(xmin = as.numeric(lower_without_gamma), xmax = as.numeric(upper_without_gamma),
                            ymin=0, ymax=ymax_val+1, alpha = "95% UI\nwithout gamma"),
                        fill='#6DBCC3')+
              geom_vline(aes(xintercept=as.numeric(bprf), linetype = "BPRF"), color="#FF0000")+
              geom_vline(aes(xintercept=as.numeric(mean_rr), linetype = "Estimated mean RR"), linetype='solid', color="#6DBCC3")+
              scale_color_manual(limits = c('1', '0'), values=c('1'='red', '0'='black'), labels = c('1' = "Trimmed", '0' = "Included"))+
              scale_linetype_manual(limits = c("Estimated mean RR", "BPRF", "Null RR"), values = c("BPRF" = "solid", "Estimated mean RR" = "solid", "Null RR" = "dashed")) +
              scale_shape_manual(limits = c('1', '0'), values=c('1'=4, '0'=16), labels = c('1' = "Trimmed", '0' = "Included"))+
              scale_alpha_manual(values = c("95% UI\nwithout gamma" = 0.5, "95% UI\nwith gamma" = 0.3)) +
              labs(title= '', 
                   alpha = "Model Uncertainty",
                   shape='Effect Sizes',
                   color = 'Effect Sizes',
                   linetype = "Model Results",
                   y='', 
                   x='')+
              theme_bw()+
              theme(plot.title=element_text(face='bold'), 
                    strip.text = element_text(face = "bold"), 
                    strip.text.x = element_blank(),
                    strip.background = element_rect(fill = "white", color = "black", size = 1, linetype = "solid"), 
                    strip.text.y.right = element_text(angle = 0), 
                    legend.title = element_text(size = 7, face = "bold"), 
                    legend.text = element_text(size = 5),
                    text = element_text(size = 5),
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
              scale_x_log10(limits = c(x_minimum, x_maximum), labels = scales::number_format(accuracy = 0.1)) +
              facet_grid(rows = vars(cause_labels), scales = "free_y", space = "free_y") +
              coord_cartesian(xlim = c(0.09, 20)) +
              guides(shape = guide_legend(order = 1), color = guide_legend(order = 1),
                     alpha = guide_legend(override.aes = list(fill = c("#6DBCC3", "#6DBCC3")), order = 2), 
                     linetype = guide_legend(override.aes = list(color = c("#6DBCC3", "#FF0000", "#000000")), order = 3)
              )
          }
          assign(paste0(versions, "_", saving_label, "_", i, "_plot2"), temp2)
          
        }
        
        
        assign(paste0(versions, "_", saving_label, "_", i, "_plot"), temp)

      }
    }
    
    
    pdf(paste0(viz_fp, versions, '_treeplot_new.pdf'), width = 7.08661, height = 4.72)
    
    for(i in unique(temp_data$age_of_exposure)){
      for(m in unique(temp_data[age_of_exposure == i]$violence_labels)){
        saving_label <- str_to_lower(gsub("\n", " ", m))
        
        title1 <- text_grob(paste0("Assocations between ", i, " ", saving_label, " and various health outcomes"), size = 7, face = "bold")
        bottom1 <- text_grob("Effect Size (95% UI)", size = 6)
        if(exists(paste0(versions, "_", saving_label, "_", i, "_plot2"))){
          legend <- get_legend(get(paste0(versions, "_", saving_label, "_", i, "_plot2")))
          temp <- get(paste0(versions, "_", saving_label, "_", i, "_plot2"))
          temp <- temp +
            theme(legend.position="none")
          layout_matrix <- rbind(c(1,1,1,1,1,1,2,2,2,2,2,3,3))
          grid.arrange(grobs = list(get(paste0(versions, "_", saving_label, "_", i, "_plot")), temp, legend), layout_matrix = layout_matrix, bottom = bottom1, 
                       top = title1, padding = unit(0.25, "line"))
        } else {
          layout_matrix <- rbind(c(1))
          grid.arrange(grobs = list(get(paste0(versions, "_", saving_label, "_", i, "_plot"))), layout_matrix = layout_matrix, bottom = bottom1, 
                       top = title1)
        }
        
      }
    }
    dev.off()
    
}

# Now lets get a heatmap comparing all of the different options
temporary <- data.table()
for(versions in c("main_models", "notrim", alts)){
  if(versions == "main_models"){
    main_models_results$analysis <- "Primary models"
    temporary <- rbindlist(list(main_models_results, temporary), fill = T)
  } else if(versions == "notrim") {
    notrim_results$analysis <- "notrim"
    temporary <- rbindlist(list(notrim_results, temporary), fill = T)
  } else {
    if(exists(paste0(versions, "_trim_results"))){
      temporary_trim <- get(paste0(versions, "_trim_results"))
      temporary_trim$analysis <- versions
      if(!exists(paste0(versions, "_notrim_results"))){
        temporary <- rbindlist(list(temporary_trim, temporary), fill = T)
      }
    }
    if(exists(paste0(versions, "_notrim_results"))){
      temporary_notrim <- get(paste0(versions, "_notrim_results"))
      temporary_notrim$analysis <- versions
      
      if(!exists(paste0(versions, "_trim_results"))){
        temporary <- rbindlist(list(temporary_notrim, temporary), fill = T)

      }
    }
    if(exists(paste0(versions, "_notrim_results")) & exists(paste0(versions, "_trim_results"))){
      temporary <- rbindlist(list(temporary_notrim, temporary_trim, temporary), fill = T)
    }
    
  }
}
  
temporary[, star_rating := as.numeric(star_rating)]
  temporary$stars <- sapply(temporary$star_rating, generate_stars)
  temporary[violence_type == "ipv", violence_labels := "Intimate partner\nsexual and/or\nphysical violence"]
  
  temporary[, violence_labels := factor(violence_labels, levels = c("Neglect", "Psychological\nabuse", "Physical\nabuse", "Sexual\nabuse", "Intimate partner\nsexual and/or\nphysical violence"))]
  temporary[, mean_rr_label := ifelse(star_rating == 0 | is.na(star_rating), "", mean_rr)]
  temporary[, full_cause_names := factor(full_cause_names, levels = cause_order)]
  temporary[, analysis_labels := ifelse(analysis == "perp_specific" & age_of_exposure == "childhood", "Only Family/\nHousehold\nPerpetrators", 
                                        ifelse(analysis == "perp_specific" & age_of_exposure == "adulthood", "Only Partner\nPerpetrators", 
                                               ifelse(analysis == "anyperp_specific", "Only Anyone/\nUnspecified\nPerpetrators", 
                                                      ifelse(analysis == "nopregnancyrecall", "No Pregnancy\nRecall", 
                                                             ifelse(analysis == "onlypregnancy", "Only Pregnancy\nRecall",
                                                             ifelse(analysis == "maleonly", "Only Male\nObservations", 
                                                                    ifelse(analysis == "femaleonly", "Only Female\nObservations", 
                                                                           ifelse(analysis == "noadjustment", "No SE\nAdjustment", 
                                                                                  ifelse(analysis == "outcome_sens", "Outcome-Specific\nAnalyses", 
                                                                                         ifelse(analysis == "notrim", "No\nTrimming", "Primary\nAnalysis"))))))))))]
  
  temporary[, count_by_outcome := 1:.N, by = c("full_cause_names", "violence_labels", "analysis_labels", "age_of_exposure")]
  
  extra_analyses <- temporary[count_by_outcome != 1]
  temporary <- temporary[count_by_outcome == 1]
  
  extra_analyses[, analysis_labels := "Extra"]
  temporary <- rbindlist(list(temporary, extra_analyses), use.names = T)
  
  temporary[, analysis_labels := factor(analysis_labels, levels = c("Primary\nAnalysis", "No SE\nAdjustment", "No\nTrimming", "Only Anyone/\nUnspecified\nPerpetrators", 
                                                                    "Only Family/\nHousehold\nPerpetrators", "Only Partner\nPerpetrators", "No Pregnancy\nRecall",
                                                                    "Only Pregnancy\nRecall", "Only Male\nObservations", "Only Female\nObservations", "Outcome-Specific\nAnalyses", "Extra"))]
  
  fill_options <- c("star_rating")
  for(filling_value in fill_options){
    temporary[, star_rating := ifelse(star_rating == 0, NA, star_rating)]
    for(types in unique(temporary$violence_labels)){
      for(i in unique(temporary[violence_labels == types]$age_of_exposure)){
        if(filling_value == "ros"){
          fill_label <- "Risk Outcome Score"
        } else if(filling_value == "star_rating"){
          fill_label <- "Star rating"
          temporary[, star_rating := ifelse(is.na(star_rating), 0, star_rating)]
        } else if(filling_value == "mean_rr"){
          fill_label <- "Mean Relative Risk"
        }
        p <- ggplot(temporary[age_of_exposure == i & violence_labels == types], aes(x = analysis_labels, y = full_cause_names, fill = as.factor(get(filling_value)))) + 
          geom_tile() 
        if(filling_value == "ros"){
          p <- p + geom_tile(data = temporary[age_of_exposure == i & is.na(star_rating)& violence_labels == types], aes(x = analysis_labels, y = full_cause_names), fill = "black") +
            geom_text(aes(label = paste0(mean_rr)), size=6, colour = "white", vjust = -0.75) + 
            geom_text(aes(label = paste0(stars)), size=10, colour = "white", vjust = 0.75)
        } else if(filling_value == "star_rating"){
          p <- p + geom_text(aes(label = paste0(mean_rr_label)), size=2.5, colour = "white")
        } else if(filling_value == "mean_rr"){
          p <- p + 
            geom_text(aes(label = paste0(mean_rr, " (", ui_with_gamma, ")")), size=6, colour = "white", vjust = -0.75) + 
            geom_text(aes(label = paste0(stars)), size = 10, colour = "white", vjust = 0.75)
        }
        p <- p + 
          geom_hline(yintercept = 0.5 + 0:35, colour = "white", linewidth = 1.5) +
          geom_vline(xintercept = 0.5 + 0:5, colour = "white", linewidth = 1.5) +
          scale_x_discrete(expand = c(0,0), position = "top") +
          scale_y_discrete(expand = c(0,0)) 
        if(length(unique(temporary[age_of_exposure == i & violence_labels == types]$star_rating)) <= 4 & min(unique(as.numeric(temporary[age_of_exposure == i & violence_labels == types]$star_rating))) == 0){
          p <- p +
            scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(1,2,4,6,8)]),
                              guide = guide_legend(reverse = TRUE))
        } else if(length(unique(temporary[age_of_exposure == i & violence_labels == types]$star_rating)) > 4){
          p <- p +
            scale_fill_manual(values = (RColorBrewer::brewer.pal(9,'Blues')[-c(2,4,6,8)]),
                              guide = guide_legend(reverse = TRUE))
        } else if(length(unique(temporary[age_of_exposure == i & violence_labels == types]$star_rating)) <= 4 & min(unique(as.numeric(temporary[age_of_exposure == i & violence_labels == types]$star_rating))) > 0){
          p <- p +
            scale_fill_manual(values = c((RColorBrewer::brewer.pal(9,'Blues')[-c(1,2,3,4,6,8)]), "#06224a"),
                              guide = guide_legend(reverse = TRUE))
        }
        p <- p +
          facet_grid(cols = vars(analysis_labels), scales = "free_x", space = "free_x") +
          theme_bw() +
          theme(plot.title = element_text(face="bold", 
                                          margin = margin(10, 0, 10, 0)), 
                panel.grid.major = element_blank(), 
                axis.text.y=element_text(face = 'bold'), 
                text = element_text(size = 5.5), 
                strip.text.x = element_blank(), 
                axis.text = element_text(size = 5)) +
          labs(title = paste0("Sensitivity analyses for the associations between ", i, " exposure to ", str_to_lower(types)," and various health outcomes"), 
               y = "", x = "", fill = fill_label)
        
        assign(paste0("combined_", i, "_heatmap_", types), p)
    }
    
    }
  }
  
  grDevices::cairo_pdf(paste0(viz_fp, 'heatmap_final_test.pdf'), width = 7.08661, height = 4.72, onefile = T)
  
  for(filling_value in fill_options){
    for(types in unique(temporary$violence_labels)){
      for(i in unique(temporary[violence_labels == types]$age_of_exposure)){
      print(get(paste0("combined_", i, "_heatmap_", types)))
    }
    }
  }
  dev.off()


