
# define directory with final results -------------------------------------

parent_bop_dir <- 'PATH TO PARENT FOLDER'

bop_dirs <- list.dirs(
  path = parent_bop_dir,
  full.names = TRUE,
  recursive = FALSE
)

names(bop_dirs) <- sapply(bop_dirs, basename, simplify = TRUE)

# helper functions --------------------------------------------------------

get_rr_curve <- function(bop_dir) {
  reticulate::use_condaenv("PATH TO BOP ENVIRONMENT")
  reticulate::source_python(
    file = '~/repos/bop_pipeline_nch/tmrel_ros/get_signal.py'
  )
  adjust_data_for_tmrel(bop_dir) |>
    data.table::setDT()
}

shift_curve <- function(dat, reference_val) {
  ref_index <- dplyr::near(
    x = reference_val,
    y = dat$risk,
    tol = 0.1
  ) |> which() |> median()
  cols_to_shift <- setdiff(colnames(dat), 'risk')
  for(c in cols_to_shift) {
    dat[[c]] <- dat[[c]] - dat[[c]][ref_index]
  }
  return(dat)
}

get_tmrel <- function(dat) {
  dat |>
    dplyr::filter(mean == min(.data$mean)) |>
    dplyr::pull("risk") |>
    unique()
}

get_ros <- function(risk_score_bounds, dat) {
  integrate(
    f = approxfun(
      x = dat$risk,
      y = dat$lower_re
    ),
    lower = risk_score_bounds[1],
    upper = risk_score_bounds[2]
  ) |>
    purrr::chuck('value') / (risk_score_bounds[2] - risk_score_bounds[1])
}

get_star_rating <- function(bprf_val) {
  if(is.na(bprf_val)) return(NA_integer_)
  score_index_vec <- c(-Inf, 0, 0.1398, 0.4055, 0.6152, Inf)
  for(i in seq_len(length(score_index_vec))) {
    next_index <- i + 1
    if (bprf_val > score_index_vec[i] && bprf_val <= score_index_vec[next_index]) {
      return(i)
    }
  }
  return(NA_integer_)
}

mean_rr_for <- function(dat, risk_range, ui_type) {
  cols_to_check <- list(
    mean = 'mean',
    lower = paste0('lower_', ui_type),
    upper = paste0('upper_', ui_type)
  )
  if(any(is.na(risk_range))) {
    return(sapply(cols_to_check, \(c) { NA_real_ }, USE.NAMES = TRUE, simplify = FALSE))
  }
  risk_index_vals <- c(
    dplyr::near(x = risk_range[1], y = dat$risk, tol = 0.1) |> which() |> min(),
    dplyr::near(x = risk_range[2], y = dat$risk, tol = 0.1) |> which() |> max()
  )
  sapply(cols_to_check, \(c) {
    dat |>
      dplyr::slice(
        risk_index_vals[1]:risk_index_vals[2]
      ) |>
      purrr::chuck(c) |>
      mean()
  }, USE.NAMES = TRUE, simplify = FALSE)
}

# get new ros based on 121 g/L tmrel --------------------------------------

draw_sets <- c('fe', 're')
names(draw_sets) <- draw_sets

tmrel_vec <- c(NA_integer_, 121)
names(tmrel_vec) <- c('from_bop', 'shifted_tmrel')

anemia_pregnancy_thresholds <- list(
  mild = c(100, 110),
  moderate = c(70, 100), 
  severe = c(0, 70)
)

bop_list <- sapply(bop_dirs, \(i) {
  print(i)
  
  rr_curve <- get_rr_curve(bop_dir = i)
  x_bop <- bophf::bop(directory = i)
  
  sapply(tmrel_vec, \(TMREL) {
    
    TMREL <- if(is.na(TMREL)) {
      get_tmrel(dat = rr_curve)
    } else {
      TMREL
    }
    
    rr_curve <- if(is.na(TMREL)) {
      rr_curve
    } else {
      shift_curve(dat = rr_curve, reference_val = TMREL)
    }
    
    brpf_value <- if(is.na(TMREL)) {
      x_bop$overall$score
    } else {
      get_ros(
        risk_score_bounds = x_bop$overall$risk_score_bounds,
        dat = rr_curve
      )
    }
    
    star_rating <- if(is.na(TMREL)) {
      x_bop$overall$star_rating
    } else {
      get_star_rating(brpf_value)
    }
    
    min_risk_value <- min(rr_curve$risk)
    
    sapply(draw_sets, \(d) {
      
      min_to_tmrel_rr <- mean_rr_for(
        dat = rr_curve,
        risk_range = c(
          min_risk_value,
          TMREL
        ), 
        ui_type = d
      )
      
      non_anemic_range_rr <- mean_rr_for(
        dat = rr_curve,
        risk_range = c(
          110,
          if(TMREL > 110) TMREL else NA_real_
        ),
        ui_type = d
      )
      
      anemia_cutoff_rrs <- sapply(anemia_pregnancy_thresholds, \(a) {
        if(a[1] < min_risk_value) a[1] <- min_risk_value
        if(a[2] < min_risk_value) a[2] <- NA_real_
        mean_rr_for(
          dat = rr_curve,
          risk_range = c(a[1], a[2]),
          ui_type = d
        )
      }, USE.NAMES = TRUE, simplify = FALSE)
      
      list(
        root_dir = i,
        tmrel = TMREL,
        min_risk_value = min_risk_value,
        ros = brpf_value,
        star_rating = star_rating,
        draw_set = d,
        min_risk_to_tmrel_rr = min_to_tmrel_rr,
        non_anemic_range_rr = non_anemic_range_rr,
        anemia_cutoff_rrs = anemia_cutoff_rrs
      )
    }, USE.NAMES = TRUE, simplify = FALSE)
  }, USE.NAMES = TRUE, simplify = FALSE)
}, USE.NAMES = TRUE, simplify = FALSE)

# get rr comparisons ------------------------------------------------------

comparison_df <- lapply(names(bop_list), \(x) {
  print(x)
  category_vec <- c('min_risk_to_tmrel_rr', 'non_anemic_range_rr', 'anemia_cutoff_rrs')
  lapply(category_vec, \(c) {
    bop_run <- c('from_bop', 'shifted_tmrel')
    lapply(bop_run, \(b) {
      dat <- data.table::data.table(
        ro_pair = x,
        analysis_category = c,
        bop_run = b,
        tmrel = bop_list[[x]][[b]][['re']][['tmrel']],
        min_risk_value = bop_list[[x]][[b]][['re']][['min_risk_value']],
        ros = bop_list[[x]][[b]][['re']][['ros']],
        star_rating = bop_list[[x]][[b]][['re']][['star_rating']],
        anemia_severity = NA_character_
      )
      if(c == 'anemia_cutoff_rrs') {
        dat <- lapply(names(anemia_pregnancy_thresholds), \(a) {
          temp_dat <- data.table::copy(dat)
          temp_dat$anemia_severity <- a
          temp_dat$lower_outer_rr_val <- bop_list[[x]][[b]][['re']][[c]][[a]][['lower']]
          temp_dat$lower_inner_rr_val <- bop_list[[x]][[b]][['fe']][[c]][[a]][['lower']]
          temp_dat$mean_rr_val <- bop_list[[x]][[b]][['re']][[c]][[a]][['mean']]
          temp_dat$upper_inner_rr_val <- bop_list[[x]][[b]][['fe']][[c]][[a]][['upper']]
          temp_dat$upper_outer_rr_val <- bop_list[[x]][[b]][['re']][[c]][[a]][['upper']]
          return(temp_dat)
        }) |> data.table::rbindlist(use.names = TRUE, fill = TRUE)
      } else {
        dat$lower_outer_rr_val <- bop_list[[x]][[b]][['re']][[c]][['lower']]
        dat$lower_inner_rr_val <- bop_list[[x]][[b]][['fe']][[c]][['lower']]
        dat$mean_rr_val <- bop_list[[x]][[b]][['re']][[c]][['mean']]
        dat$upper_inner_rr_val <- bop_list[[x]][[b]][['fe']][[c]][['upper']]
        dat$upper_outer_rr_val <- bop_list[[x]][[b]][['re']][[c]][['upper']]
      }
      return(dat)
    }) |> data.table::rbindlist(use.names = TRUE, fill = TRUE)
  }) |> data.table::rbindlist(use.names = TRUE, fill = TRUE)
}) |> data.table::rbindlist(use.names = TRUE, fill = TRUE)

# write out data ----------------------------------------------------------

output_dir <- 'PARENT TO OUTPUT FOLDER'

bop_list_file_name <- file.path(
  output_dir,
  'tmrel_comparision.rds'
)

saveRDS(
  object = bop_list,
  file = bop_list_file_name
)

bop_comparison_file_name <- file.path(
  output_dir,
  'tmrel_comparision.csv'
)

data.table::fwrite(
  x = comparison_df,
  file = bop_comparison_file_name
)

