# load in libraries -------------------------------------------------------

library(data.table)
library(ggplot2)
library(stringr)

# plotting code -----------------------------------------------------------

bop_forest_plot_w_mrbrt_publication <- function(input_df,
                                    main_title =
                                      paste("MR-BRT Analysis - Meta Analyzed", 
                                            "Relative Risk"),
                                    num_decimal_points = 3,
                                    plot_linear_space = TRUE,
                                    output_directory = NULL) {
  df <- data.table::copy(input_df)
  df <- df[order(-y_axis_title)]
  df$new_seq <- seq_len(nrow(df))
  plot_df <- data.table::data.table(
    seq = df$new_seq,
    y_axis_title = df$y_axis_title
  )
  study_rr <- new_seq <- seq <- pred_rr <- fe_se <- study_se <- re_se <- NULL
  x_axis_title <- "ln(RR)"
  label_vec <- plot_df$y_axis_title
  names(label_vec) <- as.factor(plot_df$seq)
  print(label_vec)
  rr_bound <- 0
  z_val <- qnorm(0.975)
  if (plot_linear_space) {
    plot_df$study_rr <- df$original_rr
    
    plot_df$study_lower_se <- exp(df$ln_rr - df$ln_rr_se * z_val)
    plot_df$study_upper_se <- exp(df$ln_rr + df$ln_rr_se * z_val)
    
    plot_df$pred_rr <- df$pred_rr
    
    plot_df$lower_se_fe <- exp(df$pred_ln_rr - df$pred1_ln_se * z_val)
    plot_df$upper_se_fe <- exp(df$pred_ln_rr + df$pred1_ln_se * z_val)
    
    plot_df$lower_se_re <- exp(df$pred_ln_rr - df$pred2_ln_se * z_val)
    plot_df$upper_se_re <- exp(df$pred_ln_rr + df$pred2_ln_se * z_val)
    
    x_axis_title <- "RR"
    rr_bound <- 1
  } else {
    plot_df$study_rr <- df$ln_rr
    
    plot_df$study_lower_se <- df$ln_rr - df$ln_rr_se * z_val
    plot_df$study_upper_se <- df$ln_rr + df$ln_rr_se * z_val
    
    plot_df$pred_rr <- df$pred_ln_rr
    
    plot_df$lower_se_fe <- df$pred_ln_rr - df$pred1_ln_se * z_val
    plot_df$upper_se_fe <- df$pred_ln_rr + df$pred1_ln_se * z_val
    
    plot_df$lower_se_re <- df$pred_ln_rr - df$pred2_ln_se * z_val
    plot_df$upper_se_re <- df$pred_ln_rr + df$pred2_ln_se * z_val
  }
  
  format_decimal <- function(x, k) format(round(x, k), trim=T, nsmall=k)
  
  max_x_val <- ceiling(max(plot_df$study_upper_se, na.rm = T))
  
  f <- ggplot2::ggplot(data = plot_df, ggplot2::aes(y = as.factor(seq))) +
    ggplot2::geom_point(
      ggplot2::aes(x = study_rr),
      shape = 19,
      size = 3,
      colour = "black"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        xmin = study_lower_se,
        xmax = study_upper_se
      ),
      colour = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = rr_bound,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = plot_df$pred_rr[1],
      linetype = "solid",
      colour = "black"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        xmin = lower_se_fe,
        xmax = upper_se_fe,
        y = as.numeric(plot_df$seq)
      ),
      fill = "darkslategray4",
      alpha = .3
    ) +
    ggplot2::geom_ribbon(
      data = plot_df,
      ggplot2::aes(
        xmin = lower_se_re,
        xmax = upper_se_re,
        y = as.numeric(plot_df$seq)
      ),
      fill = "darkslategray2",
      alpha = .3
    ) +
    ggplot2::scale_y_discrete(
      labels = label_vec
    ) +
    ggplot2::scale_x_continuous(
      breaks=seq(0, max_x_val, 1)
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = paste0(
        "Mean Relative Risk = ",
        round(plot_df$pred_rr[1], num_decimal_points),
        " (",
        round(plot_df$lower_se_re[1], num_decimal_points),
        " - ",
        round(plot_df$upper_se_re[1], num_decimal_points),
        ")"
      ),
      x = x_axis_title,
      y = "Study Name"
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )
  plot(f)
  if (!(is.null(output_directory))) {
    ggplot2::ggsave(
      filename = file.path(output_directory,"final_mrbrt_plot.png"),
      width = 10,
      height = 8,
      units = 'in',
      plot = f
    )
  }
}
