# Helper function to plot diagnostics before running BoP ------------------
plot_diagnostics_before_bop <- function(
  dat,
  outcome,
  analysis,
  trial,
  outputs_dir,
  outcome_severity = NULL,
  trimester = NULL,
  cov_con_socio = NULL,
  cov_rev = NULL,
  incl_undef = NULL
) {
  box::use(
    diagnostic_code / mrbrt_analysis[run_mrbrt],
    risk_hemoglobin / lib / create_mrbrt_plot_title[create_mrbrt_plot_title],
    risk_hemoglobin / lib / params[risk_dir],
  )

  cli::cli_progress_step(
    "Generating pre-BoP diagnostic plots...",
    msg_failed = "Error generating plots.",
    msg_done = "Pre-BoP data prep and plotting tasks complete."
  )

  outcome_title <- dplyr::case_when(
    grepl("lbw", outcome) ~ "low birth weight",
    outcome == "lga" ~ "large-for-gestational-age birth",
    outcome == "mat_hem" ~ "postpartum hemorrhage",
    outcome == "antenatal_hem" ~ "antepartum hemorrhage",
    outcome == "mat_mort" ~ "all-cause maternal mortality",
    outcome == "mat_sepsis" ~ "maternal sepsis",
    outcome == "neo_mort" ~ "all-cause neonatal mortality",
    outcome == "neo_sepsis" ~ "neonatal sepsis",
    outcome == "ppd" ~ "peripartum depression",
    analysis == "reverse_causal_agnostic_incl_severe_preecl_ecl" ~ "preeclampsia, severe preeclampsia, and eclampsia",
    analysis == "reverse_causal_agnostic_incl_severe_preecl" ~ "preeclampsia and severe precclampsia",
    outcome == "preecl" & !(analysis %in% c("reverse_causal_agnostic_incl_severe_preecl_ecl", "reverse_causal_agnostic_incl_severe_preecl")) ~ "preeclampsia",
    grepl("ptb", outcome) ~ "preterm birth",
    outcome == "sga" ~ "small-for-gestational-age birth",
    outcome == "stillbirth" ~ "stillbirth",
    outcome == "ghtn" ~ "gestational hypertension",
    outcome == "hdop" ~ "hypertensive disorders of pregnancy"
  )

  # MR-BRT forest plots by risk type ----
  if (any(dat$hb_risk_type %in% "low_hb")) {
    title <- create_mrbrt_plot_title(
      outcome_title = outcome_title,
      hb_risk_type = "low",
      cov_con_socio = cov_con_socio,
      cov_rev = cov_rev,
      outcome_severity = outcome_severity,
      incl_undef = incl_undef,
      trimester = trimester
    )
    title_w_break <- sub("(", "\n(", title, fixed = TRUE)

    meta_analysis_result <- run_mrbrt(input_df = dat[hb_risk_type == "low_hb"])
    mr_brt_df <- data.table::copy(meta_analysis_result$mr_brt_df)

    f <- bophf::bop_forest_plot_w_mrbrt(
      input_df = mr_brt_df,
      plot_linear_space = T,
      num_decimal_points = 2,
      main_title = title_w_break
    ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                     plot.subtitle = ggplot2::element_text(size = 10)
      )

    max_ticks <- 10
    mr_brt_df$seq <- factor(mr_brt_df$seq)
    num_levels <- nlevels(mr_brt_df$seq)
    if (num_levels > 0) {
      if (num_levels > max_ticks) {
        break_interval <- ceiling(num_levels / max_ticks)
        y_breaks <- levels(mr_brt_df$seq)[seq(1, num_levels, by = break_interval)]
      } else {
        y_breaks <- levels(mr_brt_df$seq)
      }
      f <- f + ggplot2::scale_y_discrete(breaks = y_breaks)
    }
  
    path <- file.path(outputs_dir, "mrbrt_low_hb.png")
    ggplot2::ggsave(
      filename = path,
      plot = f,
      width = 11,
      height = 8
    )
  }

  if (any(dat$hb_risk_type %in% "high_hb")) {
    title <- create_mrbrt_plot_title(
      outcome_title = outcome_title,
      hb_risk_type = "high",
      cov_con_socio = cov_con_socio,
      cov_rev = cov_rev,
      outcome_severity = outcome_severity,
      incl_undef = incl_undef,
      trimester = trimester
    )
    title_w_break <- sub("(", "\n(", title, fixed = TRUE)

    meta_analysis_result <- run_mrbrt(input_df = dat[hb_risk_type == "high_hb"])
    mr_brt_df <- data.table::copy(meta_analysis_result$mr_brt_df)

    f <- bophf::bop_forest_plot_w_mrbrt(
      input_df = mr_brt_df,
      plot_linear_space = T,
      num_decimal_points = 2,
      main_title = title_w_break
    ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                     plot.subtitle = ggplot2::element_text(size = 10)
      )

    max_ticks <- 10
    mr_brt_df$seq <- factor(mr_brt_df$seq)
    num_levels <- nlevels(mr_brt_df$seq)
    if (num_levels > 0) {
      if (num_levels > max_ticks) {
        break_interval <- ceiling(num_levels / max_ticks)
        y_breaks <- levels(mr_brt_df$seq)[seq(1, num_levels, by = break_interval)]
      } else {
        y_breaks <- levels(mr_brt_df$seq)
      }
      f <- f + ggplot2::scale_y_discrete(breaks = y_breaks)
    }
  
    path <- file.path(outputs_dir, "mrbrt_high_hb.png")
    ggplot2::ggsave(
      filename = path,
      plot = f,
      width = 11,
      height = 8
    )
  }

  # Plot RR vs. alt ----
  p <- bophf::plot_rr_vs_alt(df = dat)
  path <- file.path(outputs_dir, "rr_vs_alt.png")
  ggplot2::ggsave(
    filename = path,
    plot = p,
    width = 11,
    height = 8
  )

  # Plot frequency or domain spline knots depending on value of trial arg ----
  trials <- as.numeric(strsplit(as.character(trial), ",")[[1]])
  for (i in trials) {
    if (i %in% 1:4) {
      p <- bophf::plot_spline_knots(
        df = dat,
        knot_type = "frequency"
      )
      p <- plotly::ggplotly(p$knot_plot)
      html_path <- glue::glue("{outputs_dir}/trial_{i}_freq.html")
      htmlwidgets::saveWidget(p, html_path)
      unlink(
        stringr::str_replace(html_path, "\\.html", "_files"),
        recursive = TRUE
      )
    }
    if (i %in% 5:8) {
      # Calculate domain spline knots
      checkmate::assert_true(max(dat$ref_risk_upper) >= 110)
      knot_value_vec <- c(70, 100, 110)
      domain_spline_knots <- bophf::get_spline_knots(
        min_exposure = min(dat$alt_risk_lower),
        max_exposure = max(dat$ref_risk_upper),
        knot_value_vec = knot_value_vec
      )

      # If lowest internal knot isn't at least 0.03, then set to 0.031
      if (domain_spline_knots[2] <= 0.03) {
        domain_spline_knots[2] <- 0.031
      }

      p <- bophf::plot_spline_knots(
        df = dat,
        knot_type = "domain",
        knot_vec = domain_spline_knots
      )
      p <- plotly::ggplotly(p$knot_plot)
      html_path <- glue::glue("{outputs_dir}/trial_{i}_dom.html")
      htmlwidgets::saveWidget(p, html_path)
      unlink(
        stringr::str_replace(html_path, "\\.html", "_files"),
        recursive = TRUE
      )
    }
  }

  # Distribution of midpoints by risk type ----
  title <- paste("Distribution of Midpoints by Risk Type -", outcome_title)

  dat$alt_risk_midpoint <- (dat$alt_risk_lower + dat$alt_risk_upper) / 2
  dat$ref_risk_midpoint <- (dat$ref_risk_lower + dat$ref_risk_upper) / 2

  dat$hb_risk_type <- factor(dat$hb_risk_type, levels = c("low_hb", "high_hb"))

  df_long <- dat |>
    tidyr::pivot_longer(
      cols = c(alt_risk_midpoint, ref_risk_midpoint),
      names_to = "risk_type_midpoint",
      values_to = "midpoint"
    )

  d <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = midpoint, fill = risk_type_midpoint)
  ) +
    ggplot2::geom_histogram(
      position = "stack",
      bins = 20,
      alpha = 0.7
    ) +
    ggplot2::facet_wrap(~hb_risk_type) +
    ggplot2::labs(
      title = title,
      x = "Hb Midpoint",
      y = "Count",
      fill = "Risk Type"
    ) +
    ggplot2::theme_minimal()

  path <- file.path(outputs_dir, "hb_dist_by_risk_type.png")
  ggplot2::ggsave(
    filename = path,
    plot = d,
    width = 12,
    height = 7,
    bg = "white"
  )
}


# Run function ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 0) {
  sandbox_data_path <- args[1]
  outcome <- args[2]
  analysis <- args[3]
  trial <- args[4]
  outcome_severity <- args[5]
  trimester <- args[6]
  cov_con_socio <- args[7]
  cov_rev <- args[8]
  incl_undef <- args[9]
  outputs_dir <- args[10]

  outcome_severity <- ifelse(outcome_severity == 999, "", outcome_severity)
  trimester <- ifelse(trimester == 999, "", trimester)
  cov_con_socio <- ifelse(cov_con_socio == 999, "", cov_con_socio)
  cov_rev <- ifelse(cov_rev == 999, "", cov_rev)
  incl_undef <- ifelse(incl_undef == 999, "", incl_undef)

  dat <- data.table::fread(sandbox_data_path)
  plot_diagnostics_before_bop(
    dat = dat,
    outcome = outcome,
    analysis = analysis,
    trial = trial,
    outcome_severity = outcome_severity,
    trimester = trimester,
    cov_con_socio = cov_con_socio,
    cov_rev = cov_rev,
    incl_undef = incl_undef,
    outputs_dir = outputs_dir
  )
} else {
  message("No commandArgs() found; no plot_diagnostics_before_bop() call run")
}
