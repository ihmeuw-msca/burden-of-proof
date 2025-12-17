#' Use various BoP results files to create a summary PowerPoint
#'
#' @param map A map containing columns `outcome`, `analysis`, `trial_dir`, `include_mrbrt` (1 or 0), `inlier` (1.0 or 0.9), and `notes`
#' @param risk A string containing the risk name to be used in slide titles.

create_summary_slides <- function(map, risk) {
  # Add columns to map
  map <- map |>
    dplyr::mutate(
      outputs_dir = file.path(trial_dir),
      outcome_title = dplyr::case_when(
        grepl("lbw", outcome) ~ "Low Birth Weight",
        outcome == "lga" ~ "Large-for-Gestational-Age",
        category == "maternal" & grepl("hem", outcome) ~ "Postpartum Hemorrhage",
        outcome == "mat_mort" ~ "Maternal Mortality",
        outcome == "mat_sepsis" ~ "Maternal Sepsis",
        outcome == "neo_mort" ~ "Neonatal Mortality",
        outcome == "neo_sepsis" ~ "Neonatal Sepsis",
        outcome == "ppd" ~ "Peripartum Depression",
        outcome == "preecl" ~ "Preeclampsia",
        grepl("ptb", outcome) ~ "Preterm Birth",
        outcome == "sga" ~ "Small-for-Gestational-Age",
        outcome == "stillbirth" ~ "Stillbirth",
        outcome == "ghtn" ~ "Gestational Hypertension",
        outcome == "hdop" ~ "Hypertensive Disorders of Pregnancy",
        TRUE ~ outcome
      ),
      slide_title = dplyr::case_when(
        outcome %in% c("vlbw", "vptb") ~ paste("Very", outcome_title),
        outcome %in% c("elbw", "eptb") ~ paste("Extremely", outcome_title),
        outcome == "mptb" ~ paste("Moderate", outcome_title),
        outcome == "antenatal_hem" ~ "Antepartum Hemorrhage ",
        outcome == "postnatal_dep" ~ paste("Postnatal Depression", outcome_title),
        outcome == "prenatal_dep" ~ paste("Prenatal Depression", outcome_title),
        TRUE ~ outcome_title
      ),
      slide_title = paste0(
        slide_title,
        " | Analysis: ", analysis,
        " | Trial ", trial_no,
        " | Inlier = ", sprintf("%.1f", inlier),
        "\n Notes: ", notes
      ),
      image_path = NA_character_  # Initialize column
    )
  
  shown_mrbrt <- list()
  
  for (i in seq_len(nrow(map))) {
    map_row <- map[i, ]
    
    # Check if trial_dir exists and summary_updated.yaml exists
    has_valid_trial_dir <- isTRUE(!is.na(map_row$trial_dir) &&
                                    file.exists(file.path(map_row$trial_dir, "summary_updated.yaml")))
    
    out_analysis <- paste(map_row$outcome, map_row$analysis, sep = "_")
    output_path <- file.path(path.expand("~/test_images"), paste0("slide_", i, ".png"))
    
    # Helper: check MR-BRT image availability
    include_mrbrt_available <- isTRUE(
      map_row$include_mrbrt == 1 &&
        !(out_analysis %in% shown_mrbrt) &&
        file.exists(file.path(map_row$base, "mrbrt_low_hb.png"))
    )
    
    if (!has_valid_trial_dir) {
      if (include_mrbrt_available) {
        image_dir <- map_row$base
        
        mrbrt_low_img <- magick::image_read(file.path(image_dir, "mrbrt_low_hb.png"), density = 300)
        mrbrt_high_img <- tryCatch({
          magick::image_read(file.path(image_dir, "mrbrt_high_hb.png"), density = 300)
        }, error = function(e) NULL)
        
        low_grob <- grid::rasterGrob(mrbrt_low_img, interpolate = TRUE)
        mrbrt_layout <- if (!is.null(mrbrt_high_img)) {
          high_grob <- grid::rasterGrob(mrbrt_high_img, interpolate = TRUE)
          cowplot::plot_grid(NULL, low_grob, high_grob, NULL, NULL,
                             ncol = 5, nrow = 1,
                             rel_widths = c(0.2, 1, 1, 0.25, 0.05),
                             align = "hv")
        } else {
          cowplot::plot_grid(low_grob, ncol = 1, align = "hv")
        }
        
        save_and_crop_plot(mrbrt_layout, output_path, crop = TRUE, title = map_row$slide_title)
        shown_mrbrt <- c(shown_mrbrt, out_analysis)
        map$image_path[i] <- output_path
        next
      } else {
        # No trial_dir and no MR-BRT images; skip slide
        next
      }
    }
    
    # Normal processing if trial_dir valid
    p1 <- create_summary_table(map_row$outputs_dir) |> flextable::gen_grob(fit = "fixed")
    p2 <- create_forest_img(map_row$outputs_dir) |> grid::rasterGrob(width = grid::unit(1, "npc"), hjust = 0.48, vjust = -0.3)
    p3 <- create_tmrel_table(map_row$outputs_dir) |> flextable::gen_grob(fit = "fixed")
    
    if (include_mrbrt_available) {
      image_dir <- map_row$base
      
      mrbrt_low_img <- magick::image_read(file.path(image_dir, "mrbrt_low_hb.png"), density = 300)
      mrbrt_high_img <- tryCatch({
        magick::image_read(file.path(image_dir, "mrbrt_high_hb.png"), density = 300)
      }, error = function(e) NULL)
      
      low_grob <- grid::rasterGrob(mrbrt_low_img, interpolate = TRUE)
      mrbrt_layout <- if (!is.null(mrbrt_high_img)) {
        high_grob <- grid::rasterGrob(mrbrt_high_img, interpolate = TRUE)
        cowplot::plot_grid(NULL, low_grob, high_grob, NULL, NULL,
                           ncol = 5, nrow = 1,
                           rel_widths = c(0.2, 1, 1, 0.25, 0.05),
                           align = "hv")
      } else {
        cowplot::plot_grid(low_grob, ncol = 1, align = "hv")
      }
      
      p2 <- create_forest_img(map_row$outputs_dir) |> grid::rasterGrob(width = grid::unit(1, "npc"), hjust = 0.48, vjust = 0.3)
      grid_layout <- cowplot::plot_grid(NULL, p1, p2, NULL, p3, align = "hv", ncol = 5, nrow = 1,
                                        rel_widths = c(0.05, 1, 1.4, 0.05, 1), expand = FALSE)
      combined_plot <- cowplot::plot_grid(grid_layout, mrbrt_layout, ncol = 1, rel_heights = c(2, 1.7), align = "v")
      save_and_crop_plot(combined_plot, output_path, crop = TRUE, title = map_row$slide_title)
      shown_mrbrt <- c(shown_mrbrt, out_analysis)
    } else {
      grid_layout <- cowplot::plot_grid(NULL, p1, p2, NULL, p3, align = "hv", ncol = 5, nrow = 1,
                                        rel_widths = c(0.05, 1, 1.3, 0.05, 1), expand = FALSE)
      save_and_crop_plot(grid_layout, output_path, crop = TRUE, title = map_row$slide_title)
    }
    
    map$image_path[i] <- output_path
  }
  
  return(map)
}

# Helper to save/crop flextable as image
save_and_crop_plot <- function(plot, path, crop = TRUE, title) {
  # Save to disk
  cowplot::save_plot(filename = path, plot = plot, base_width = 10, base_height = 4)
  
  if (crop) {
    # Read and crop image using magick
    img <- magick::image_read(path)
    img_cropped <- magick::image_trim(img) # Automatically trims whitespace
    magick::image_write(img_cropped, path)
  }
}

create_summary_table <- function(trial_dir) {
  bop_summary_fp <- file.path(trial_dir, "summary_updated.yaml")
  
  # Generate table of ROS, RR, and star rating using values from summary YAML
  summary_values_list <- yaml::yaml.load_file(bop_summary_fp)
  
  # Helper to extract values or return NA
  safe_extract <- function(list, path) {
    purrr::pluck(list, !!!path, .default = NA)
  }
  
# Add line break before the parentheses
  insert_linebreak <- function(x) {
    if (is.null(x)) return(NA_character_)
    sub(" \\(", "\n(", x)
  }
  
  # Paste risk bounds as "min - max"
  format_bounds <- function(bounds) {
    if (is.null(bounds) || length(bounds) < 2) return(NA_character_)
    paste(round(as.numeric(bounds[1]), 2), "-", round(as.numeric(bounds[2]), 2))
  }
  
  summary_df <- tibble::tibble(
    ` ` = c(
      "ROS bound (g/L)",
      "Mean RR (95% CI), fe",
      "Mean RR (95% CI), re",
      "Star rating"
    ),
    `Overall` = c(
      format_bounds(safe_extract(summary_values_list, c("overall", "risk_score_bounds"))),
      insert_linebreak(safe_extract(summary_values_list, c("overall", "mean_rr_fe"))),
      insert_linebreak(safe_extract(summary_values_list, c("overall", "mean_rr_re"))),
      safe_extract(summary_values_list, c("overall", "star_rating"))
    ),
    `<TMREL` = c(
      format_bounds(safe_extract(summary_values_list, c("lower", "risk_score_bounds"))),
      insert_linebreak(safe_extract(summary_values_list, c("lower", "mean_rr_fe"))),
      insert_linebreak(safe_extract(summary_values_list, c("lower", "mean_rr_re"))),
      safe_extract(summary_values_list, c("lower", "star_rating"))
    ),
    `>TMREL` = c(
      format_bounds(safe_extract(summary_values_list, c("upper", "risk_score_bounds"))),
      insert_linebreak(safe_extract(summary_values_list, c("upper", "mean_rr_fe"))),
      insert_linebreak(safe_extract(summary_values_list, c("upper", "mean_rr_re"))),
      safe_extract(summary_values_list, c("upper", "star_rating"))
    )
  )
  
  summary_table <- summary_df |>
    flextable::flextable() |>
    flextable::height_all(height = .4) |>
    flextable::hrule(rule = "exact") |>
    flextable::width(width = c(0.75, 0.75, 0.75, 0.75)) |>
    flextable::align(align = "left", part = "body") |>
    flextable::align(align = "center", part = "header") |>
    flextable::fontsize(size = 7, part = "all") |>
    flextable::set_table_properties(layout = "fixed")

  summary_table
}

#' Helper function to prepare plot image
create_forest_img <- function(trial_dir) {
  bop_forest_fp <- list.files(
    path = trial_dir,
    pattern = ".*-.*\\.pdf$",
    full.names = TRUE
  )
  bop_forest_img <- magick::image_read_pdf(bop_forest_fp, density = 300)
  bop_forest_img <- magick::image_crop(bop_forest_img, geometry = "1630x750+0+0")
  
  bop_forest_img_fp <- sub(".pdf", ".png", bop_forest_fp)
  magick::image_write(bop_forest_img, bop_forest_img_fp)
  bop_forest_img <- png::readPNG(bop_forest_img_fp)

  bop_forest_img
}


#' Helper function to prepare TMREL and NID values values from summary YAML
create_tmrel_table <- function(trial_dir) {
  bop_summary_fp <- file.path(trial_dir, "summary_updated.yaml")
  summary_values_list <- yaml::yaml.load_file(bop_summary_fp)

  tmrel_and_nid_counts <- tibble::tibble(
    `TMREL & NID Obs counts` = c(
      paste("TMREL:", round(summary_values_list$tmrel)),
      paste("Total NID#:", summary_values_list$unique_nids_total),
      paste("NID# (low-Hb):", summary_values_list$unique_nids_lowhb),
      paste("NID# (high-Hb):", summary_values_list$unique_nids_highhb),
      paste("Total Obs#:", summary_values_list$obs_total),
      paste("Obs# (low-Hb):", summary_values_list$obs_lowhb),
      paste("Obs# (high-Hb):", summary_values_list$obs_highhb)
    )
  )

  tmrel_and_nid_table <- tmrel_and_nid_counts |>
    flextable::flextable(col_keys = names(tmrel_and_nid_counts)) |>
    flextable::height(height = .2) |>
    flextable::hrule(rule = "exact", part = "body") |>
    flextable::width(width = 1.2) |>
    flextable::align(align = "left", part = "body") |>
    flextable::align(align = "center", part = "header") |>
    flextable::fontsize(size = 7) |>
    flextable::set_table_properties(layout = "fixed")

  tmrel_and_nid_table
}
