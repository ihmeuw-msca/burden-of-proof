# Load the necessary libraries
library(ggplot2)
library(scales)
library(data.table)


# Load data
# MATERNAL 
### WITH GAMMA
a <- fread('PATH TO DRAWS DATASET')
summary(a$risk)


# NEONATAL
### WITH GAMMA
a <- fread('PATH TO DRAWS DATASET')
summary(a$risk)


# summarize draws
a_summary <- nch::pivot_draws_longer(a) |>
  nch::summarize_draws()

# find the lowest ln(RR) value
tmrel_index <- which.min(a_summary$mean)

# normalize summarized draws to TMREL
a_summary$normalized_mean <- a_summary$mean - a_summary$mean[tmrel_index]
a_summary$normalized_lower <- a_summary$upper - a_summary$upper[tmrel_index]
a_summary$normalized_upper <- a_summary$lower - a_summary$lower[tmrel_index]


# calculate mean_rr (linear space)
a_summary$normalized_mean_rr <- exp(a_summary$normalized_mean)
a_summary$normalized_lower_rr <- exp(a_summary$normalized_lower)
a_summary$normalized_upper_rr <- exp(a_summary$normalized_upper)

# subset
b <- a_summary[normalized_upper_rr>=0.99 & normalized_upper_rr<=1.01]

# plot normalized_lower_rr
ggplot(b, aes(x = risk, y = normalized_upper_rr)) +
  geom_point() +
  labs(
    x = "Midpoint of alternate hemoglobin range (g/L)",
    y = "Lower RR (with gamma)",
    title = "All-cause total neonatal mortality\nHb level vs. lower RR (with gamma), where RR= 1 +/- 0.01"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))


# Load necessary libraries
title_text <- "All-cause total neonatal mortality\nHb level vs. lower RR (with gamma), where RR= 1 +/- 0.01"
title_text <- "All-cause total maternal mortality\nHb level vs. lower RR (with gamma), where RR= 1 +/- 0.01"


# Load necessary libraries
library(ggplot2)
library(scales)

# Identify the risk values where normalized_upper_rr equals 1, 1.01, and 1.001
risk_rr_1 <- b$risk[which.min(abs(b$normalized_upper_rr - 1))]
risk_rr_1_01 <- b$risk[which.min(abs(b$normalized_upper_rr - 1.01))]
risk_rr_1_001 <- b$risk[which.min(abs(b$normalized_upper_rr - 1.001))]

# Create the plot with vertical lines and annotations
ggplot(b, aes(x = risk, y = normalized_upper_rr)) +
  geom_point() +
  labs(
    x = "Midpoint of alternate hemoglobin range (g/L)",
    y = "Lower RR (with gamma)",
    title = "All-cause total neonatal mortality\nHb level vs. lower RR (with gamma), where RR= 1 +/- 0.01"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 
# geom_vline(xintercept = risk_rr_1, linetype = "dashed", color = "red") +
# geom_vline(xintercept = risk_rr_1_01, linetype = "dashed", color = "blue") +
# geom_vline(xintercept = risk_rr_1_001, linetype = "dashed", color = "purple") +
# annotate("text", x = risk_rr_1, y = 1.008, label = paste("Lower RR = 1 at", round(risk_rr_1, 2)), vjust = 2, color = "red") +
# annotate("text", x = 130, y = Inf, label = paste("Lower RR = 1.01 at", round(risk_rr_1_01, 2)), vjust = 2, color = "blue") +
# annotate("text", x = risk_rr_1_001, y = Inf, label = paste("Lower RR = 1.001 at", round(risk_rr_1_001, 2)), vjust = 2, color = "purple")


# calculate the slopes of the lower UI between each risk value
a_summary$slope <- NA_real_
a_summary$slope <- sapply(1:10000, \(index) {
  prev_index <- index - 1
  
  if(index == 1) return(NA_real_)
  
  y2 <- a_summary$normalized_lower[index]
  y1 <- a_summary$normalized_lower[prev_index]
  
  x2 <- a_summary$risk[index]
  x1 <- a_summary$risk[prev_index]
  
  (y2 - y1) / (x2 - x1)
}, simplify = TRUE)

# caluclate the difference between the slopes
a_summary$slope_diff <- NA_real_
a_summary$slope_diff <- sapply(1:10000, \(index) {
  prev_index <- index - 1
  
  if(index == 1) return(NA_real_)
  
  y2 <- a_summary$slope[index]
  y1 <- a_summary$slope[prev_index]
  
  y2 - y1
}, simplify = TRUE)

# set threshold values to defined a "flat, non-changing line"
threshold_slope_diff_val <- 1 * 10 ^ (-8) # non-changing = where the slopes have stabalized 
threshold_val_slope_val <- 2 * 10 ^ (-4) # flat = where the slops isn't increasing nor decreasing
a_summary$flat_slope <- NA
a_summary$flat_slope <- sapply(1:10000, \(index) {
  if(!(is.na(a_summary$slope_diff[index])) && abs(a_summary$slope_diff[index]) < threshold_slope_diff_val
     && !(is.na(a_summary$slope[index])) && abs(a_summary$slope[index]) < threshold_val_slope_val) {
    return(TRUE)
  }
  FALSE
}, simplify = TRUE)

# get the tmrel based on the above criteria
tmrel <- min(
  a_summary |>
    dplyr::filter(flat_slope == TRUE) |>
    purrr::chuck('risk')
)


# plot results
library(ggplot2)

title_text <- 'First Derivative of RR in Neonatal Mortality'
title_text <- 'First Derivative of RR in Maternal Mortality'

ggplot(a_summary[!(is.na(slope))], aes(x = risk, y = slope)) +
  geom_point() +
  geom_vline(xintercept = tmrel, colour = 'red', linetype = 'dashed') +
  labs(
    title = title_text,
    subtitle = paste('TMREL =', round(tmrel, 2)),
    x = 'Hemoglobin (g/L)',
    y = 'First Derivative of ln(RR) Lower UI'
  ) +
  theme_minimal()

