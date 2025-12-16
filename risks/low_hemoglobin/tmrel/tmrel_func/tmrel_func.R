inner_file_list <- list.files(
  path = 'PATH TO PARENT FOLDER',
  recursive = TRUE,
  pattern = 'inner_draws.csv',
  full.names = TRUE
)

inner_file_list <- grep(
  pattern = "00_Archive",
  x = inner_file_list,
  value = TRUE,
  invert = TRUE
) 

inner_file_list <- grep(
  pattern = 'results_updated',
  x = inner_file_list,
  value = TRUE
)

# helper functions --------------------------------------------------------

get_directory <- function(file_name) {
  file_vec <- stringr::str_split(
    string = file_name,
    pattern = '/',
    simplify = TRUE
  )
  file_vec <- file_vec[-length(file_vec)]
  return(paste(file_vec, collapse = '/'))
}

# find tmrel values for all bop results -----------------------------------

tmrel_vec <- c()
for(f in inner_file_list) {
  print(f)
  parent_dir <- get_directory(f)
  
  summary_yaml <- yaml::read_yaml(
    file = file.path(parent_dir, 'summary.yaml')
  )
  
  bop_tmrel <- summary_yaml$tmrel_hb_value
  lower_hb_bound <- 110
  upper_hb_bound <- 135
  
  if(bop_tmrel < lower_hb_bound || bop_tmrel > upper_hb_bound) {
  
    a <- data.table::fread(f) 
    
    # summarize draws
    a_summary <- nch::pivot_draws_longer(a) |>
      nch::summarize_draws()
    
    # find the lowest ln(RR) value
    tmrel_index <- which.min(a_summary$mean)
    
    # normalize summarized draws to TMREL
    a_summary$normalized_mean <- a_summary$mean - a_summary$mean[tmrel_index]
    a_summary$normalized_lower <- a_summary$upper - a_summary$upper[tmrel_index]
    a_summary$normalized_upper <- a_summary$lower - a_summary$lower[tmrel_index]
    
    # create the risk function
    slope_func <- smooth.spline(
      x = a_summary$risk,
      y = a_summary$normalized_mean
    )
    
    # plot summarized draws vs approximated function
    plot(a_summary$risk, a_summary$normalized_mean)
    plot(a_summary$risk, predict(slope_func, deriv = 0)$y)
    
    # take the first derivative and evaluate the risks of the RR function
    first_deriv <- predict(slope_func, deriv = 1)
    a_summary$slope <- first_deriv$y

    # Plot first derivative with vertical line at bop_tmrel
    plot(a_summary$risk, a_summary$slope, main = "First Derivative of RR Function", xlab = "Hb (g/L)", ylab = "Slope")
    abline(v = bop_tmrel, col = "red", lwd = 2, lty = 2)
    mtext(paste("TMREL =", round(bop_tmrel, 2)), side = 3, line = 0.5, cex = 0.8)
    mtext(f, side = 1, line = 4, cex = 0.7, col = "blue")

    # take the second derivative and evaluate the risks of the RR function
    second_deriv <- predict(slope_func, deriv = 2)
    a_summary$inflection <- second_deriv$y

    # Plot second derivative with vertical line at bop_tmrel
    plot(a_summary$risk, a_summary$inflection, main = "Second Derivative of RR Function", xlab = "Hb (g/L)", ylab = "Inflection")
    abline(v = bop_tmrel, col = "red", lwd = 2, lty = 2)
    mtext(paste("TMREL =", round(bop_tmrel, 2)), side = 3, line = 0.5, cex = 0.8)
    mtext(f, side = 1, line = 4, cex = 0.7, col = "blue")
    
    # isolate the range where we expect a plausible TMREL
    to_check_index_vec <- which(a_summary$risk >= lower_hb_bound & a_summary$risk <= upper_hb_bound)
    non_anemic_index <- min(to_check_index_vec) # get the starting row of this range
    
    # define the "flat slope" as the slope "1 SD away from the mean of the absolute values of all slopes in this range"  
    # anything below this slope will be considered flat
    flat_slope <- as.numeric(quantile(abs(a_summary$slope[to_check_index_vec]), pnorm(-2)))
    
    # get the value of the second derivative at the starting index.
    # assuming the RR function is becoming asymptotically closer to 0, we expect the second derivative from the function calculated 
    # above to flatten out by bouncing below and above 0 before finally being "0".
    # so we get the current value of the second derivative at the starting index, and once we see it cross 0, we assume it starts
    # to become flat
    root_second_derivative <- sign(a_summary$inflection[(non_anemic_index - 1)])
    
    # create flags for identifying the new TMREL
    found_inflection <- FALSE
    found_good_point <- FALSE
    
    # loop through all exposure values in the plausible TMREL range
    for(i in to_check_index_vec) {
      current_second_derivative <- sign(a_summary$inflection[i])
      if(!found_inflection && root_second_derivative != current_second_derivative) {
        found_inflection <- TRUE
      }
      current_slope <- abs(a_summary$slope[i])
      
      # if we have found where the inflection point has crossed 0 and have a slope that is 'flat',
      # break the for loop and define the current risk as the TMREL
      if(found_inflection && current_slope < flat_slope) {
        tmrel_vec <- append(tmrel_vec, a_summary$risk[i])
        found_good_point <- TRUE
        print(a_summary$risk[i])
        break
      }
    }
  } else {
    found_good_point <- TRUE
    tmrel_vec <- append(tmrel_vec, bop_tmrel)
  }
  
  if(!found_good_point) {
    print('no tmrel found')
  }
  
  cat('\n\n')
}

library(ggplot2)

ggplot(data = data.frame(x = a_summary$risk, y = a_summary$normalized_mean), aes(x = x, y = y)) +
  geom_line() +
  geom_vline(xintercept = 116.85, colour = 'red')
