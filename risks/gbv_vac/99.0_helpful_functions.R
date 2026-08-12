## SCRIPT ##################################################################
outcome_heatmaps <- function(total_dataset, plotting_dataset = NA, total_only = F, perp_type, outcome_level, perp_levels, viol_type, fill_version = "binary"){
  if(total_only == T){
    plotting_data <- copy(total_dataset)
  } else {
    plotting_data <- copy(plotting_dataset)
  }
  
  if(outcome_level == "main"){
    if(total_only){
      setnames(plotting_data, c("perp_group", "violence_type"), c("summary", "perp_group"))
      title_version <- paste0("Number of studies by outcome and violence type")
      caption_version <- paste0("Main outcomes")
    } else {
      setnames(plotting_data, c(paste0("count_by_outcome_", perp_type), paste0("perp_group_", perp_type)), c("count_by_outcome", "perp_group"))
      title_version <- paste0("Number of studies with ", viol_type, " data by outcome and perpetrator type")
      caption_version <- paste0("Perpetrator groupings used: ", perp_type, " v. non-", perp_type, "\n Main outcomes")
    }

  } else {
    if(total_only){
      plotting_data[, `:=` (count_by_outcome = NULL, outcome = NULL)]
      setnames(plotting_data, c(paste0("count_by_level", outcome_level), paste0("outcome_grouping_level", outcome_level), "perp_group", "violence_type"), c("count_by_outcome", "outcome", "summary", "perp_group"))
      title_version <- paste0("Number of studies by outcome and violence type")
      caption_version <- paste0("Level ", outcome_level, " GBD outcomes and non-GBD outcomes")
      
    } else {
      plotting_data[, `:=` (count_by_outcome = NULL, outcome = NULL)]
      total_dataset[, `:=` (count_by_outcome = NULL, outcome = NULL)]
      setnames(plotting_data, c(paste0("count_by_level", outcome_level, "_", perp_type), paste0("perp_group_", perp_type), paste0("outcome_grouping_level", outcome_level)), c("count_by_outcome", "perp_group", "outcome"))
      setnames(total_dataset, c(paste0("count_by_level", outcome_level), paste0("outcome_grouping_level", outcome_level)), c("count_by_outcome", "outcome"))
      title_version <- paste0("Number of studies with ", viol_type, " data by outcome and perpetrator type")
      caption_version <- paste0("Perpetrator groupings used: ", perp_type, " v. non-", perp_type, "\n Level ", outcome_level, " GBD outcomes and non-GBD outcomes")
    }
    
  }

  if(total_only){
    xlabel <- "Violence type"
  } else {
    plotting_data <- rbindlist(list(plotting_data, total_dataset), fill = T)
    xlabel <- "Perpetrator group"
  } 

  if(fill_version != "binary"){
    p <- ggplot(plotting_data[!is.na(outcome)]) + 
      geom_tile(aes(x = perp_group, y = outcome, fill = count_by_outcome)) +
      geom_text(aes(x = perp_group, y = outcome, label = count_by_outcome), color = "white") + 
      theme_bw() + 
      scale_x_discrete(limits = perp_levels) +
      labs(title = title_version, caption = caption_version, x = xlabel, y = "Outcome", fill = "Unique studies") + 
      geom_vline(xintercept = 1.5, color = "white")
  } else if(fill_version == "binary"){
    plotting_data[!is.na(count_by_outcome) & count_by_outcome != "", outcome_limit := ifelse(count_by_outcome >= 3, "Yes", "No")]
    
    p <- ggplot(plotting_data[!is.na(outcome)]) + 
      geom_tile(aes(x = perp_group, y = outcome, fill = outcome_limit)) +
      geom_text(aes(x = perp_group, y = outcome, label = count_by_outcome), color = "white") + 
      theme_bw() + 
      scale_x_discrete(limits = perp_levels) +
      labs(title = title_version, caption = caption_version, x = xlabel, y = "Outcome", fill = "Enough studies to analyze") 
    if(total_only){
      p <- p + 
      geom_vline(xintercept = 1.5, color = "white") + 
        geom_vline(xintercept = 2.5, color = "white") + 
        geom_vline(xintercept = 3.5, color = "white")
    } else {
      p <- p + 
        geom_vline(xintercept = 1.5, color = "white")
    }
  }
  return(p)
}

removing_untestable_covariates <- function(data, verbose = T){
  #find covariates with minimum studies represented
  candidate_covs_pot <- names(data)[names(data) %like% 'cov_']
  
  candidate_covs_one <- c()
  candidate_covs <- c()
  for (c in candidate_covs_pot) {
    if (length(unique(data[get(c)==1]$Study_ID))>1 & length(unique(data[get(c)==0]$Study_ID))>1) {
      candidate_covs_one <- c(candidate_covs_one, c)
    }
  }
  
  covs_to_remove <- setdiff(candidate_covs_pot, candidate_covs_one)
  if(verbose){
    print(paste0("Potential_covariates: ", paste0(candidate_covs_pot, collapse = ", ", sep = "")))
    print(paste0("Covariates with minimum number of studies: ", paste0(candidate_covs_one, collapse = ", ", sep = "")))
  }
  data[, c(covs_to_remove):=NULL]
  
  #check for identical covariates
  if (length(candidate_covs_one)>0){
    sub2 <- copy(as.data.frame(data[, .SD, .SDcols=candidate_covs_one]))
    duplicated_columns <- duplicated(t(sub2))
    
    # Show the Names of the Duplicated Columns
    dup_covs <- colnames(sub2[duplicated_columns])
  } else {dup_covs <- c()}
  
  # Remove the Duplicated Columns
  if(length(dup_covs)>0){
    data[, c(dup_covs):=NULL]
    if(verbose){
      print(paste0("Covariates removed as duplicates: ", paste0(dup_covs, collapse = ", ", sep = "")))
    }
  }
  return(data)
}

