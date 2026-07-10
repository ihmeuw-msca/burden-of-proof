##' *************************************************************************************
##' Title: update_location.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Steps: 
##' 1. get location ids
##' 2. update for GBD 2023 location IDs (includes fix for UK/location_id=95)
##' 3. reassign location IDs that aren't in GBD 2023 
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)


# get location IDs --------------------------------------------------------

get_location_id <- function(input_df, gbd_rel_id){
  df <- copy(input_df)
  
  loc_df <- ihme::get_location_metadata(
    location_set_id = 35,
    release_id = gbd_rel_id
  )
  
  loc_df <- loc_df[,.(location_id, ihme_loc_id)]
  
  df$ihme_loc_id <- clean_ihme_loc_id(df$ihme_loc_id)
  
  i_vec <- which(is.na(df$location_id) & !(is.na(df$ihme_loc_id)))
  inverse_i_vec <- setdiff(seq_len(nrow(df)), i_vec)
  
  to_update_df <- df[i_vec, ]
  df <- df[inverse_i_vec, ]
  to_update_df$location_id <- NULL
  
  to_update_df <- merge.data.table(
    x = to_update_df,
    y = loc_df,
    by = "ihme_loc_id",
    all.x = TRUE
  )
  
  df <- rbindlist(
    list(df, to_update_df),
    use.names = TRUE, 
    fill = TRUE
  )
  
  if(any(is.na(df$location_id))){
    stop("Not all locations are defined. Please check and rerun.")
  }
  
  return(df)
}

clean_ihme_loc_id <- function(vec){
  vec <- str_remove_all(string = vec, pattern = "\\n")
  vec <- str_remove_all(string = vec, pattern = "\\r")
  vec <- str_remove_all(string = vec, pattern = "\\t")
  return(vec)
}

# update for GBD 2023 location IDs ----------------------------------------

update_gbd2023_location_ids <- function(input_df){
  df <- copy(input_df)
  
  df$orig_location_id <- df$location_id
  df<-df[orig_location_id == 95 & (smaller_site_unit == 0 | 
                                     grepl("eng", tolower(field_citation_value)) | 
                                     grepl("eng", tolower(site_memo))), location_id := 4749]
  df$location_id <- nch::convert_2021_locations_to_2023(df$location_id)
  df$location_name <- nch::name_for('location', id = df$location_id)
  
  return(df)
}

split_ethiopia_subnats <- function(input_df){
  df <- copy(input_df)
  
  i_vec <- which(df$location_id == 44858)
  inverse_i_vec <- setdiff(seq_len(nrow(df)), i_vec)
  
  eth_df <- df[i_vec, ]
  df <- df[inverse_i_vec, ]
  
  subnat_id_vec <- c(60908, 95069, 94364)
  for(i in subnat_id_vec){
    temp_df <- copy(eth_df)
    temp_df$location_id <- i
    temp_df$sample_size <- temp_df$sample_size / length(subnat_id_vec)
    temp_df$cases <- temp_df$cases / length(subnat_id_vec)
    df <- rbindlist(
      list(df, temp_df),
      use.names = TRUE, 
      fill = TRUE
    )
  }
  return(df)
}

# reassign location IDs that aren't in GBD 2023 ---------------------------

reassign_gbd_locations <- function(input_df, release_id){
  df <- copy(input_df)
  
  loc_df <- ihme::get_location_metadata(location_set_id = 35, release_id = release_id)
  
  i_vec <- which(!(df$location_id %in% loc_df$location_id))
  inverse_i_vec <- setdiff(seq_len(nrow(df)), i_vec)
  
  to_update_df <- df[i_vec, ]
  df <- df[inverse_i_vec, ]
  
  for(r in seq_len(nrow(to_update_df))){
    temp_row <- to_update_df[r, ]
    loc_index <- NA_integer_
    if(temp_row$location_name %in% loc_df$location_name){
      loc_index <- match(temp_row$location_name, loc_df$location_name)
    }else if(temp_row$ihme_loc_id %in% loc_df$ihme_loc_id){
      loc_index <- match(temp_row$ihme_loc_id, loc_df$ihme_loc_id)
    }else if(substr(temp_row$ihme_loc_id, 1, 3) %in% loc_df$ihme_loc_id){
      loc_index <- match(substr(temp_row$ihme_loc_id, 1, 3), loc_df$ihme_loc_id)
    }
    
    if(!(is.na(loc_index))){
      temp_row$location_id <- loc_df$location_id[loc_index]
      df <- rbindlist(
        list(df, temp_row),
        use.names = TRUE,
        fill = TRUE
      )
    }
  }
  return(df)
}
