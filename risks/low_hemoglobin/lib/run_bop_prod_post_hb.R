#' Helper function to run post-processing pipeline steps for bop-prod
run_bop_prod_post_hb <- function(
  dir,
  ro_pair,
  analysis,
  update_output
) {
  box::use(
    risk_hemoglobin / lib / update_bop_summary,
    risk_hemoglobin /
      tmrel /
      tmrel_func /
      run_find_tmrel_func_quantiles[run_find_tmrel_func_quantiles, find_tmrel],
  )

  # Update summary.yaml files in each of the BoP-dev trial folders
  if (update_output == TRUE) {
    split_method <- "simple"
    # Run with split at TMREL--if unsuccessful, tries again w/out `split_at` arg
    update_bop_summary$update_summary(
      parent_dir = dir,
      split_at = "tmrel",
      split_method = split_method
    )

    # Write additional files with risk split at 120 g/L
    update_bop_summary$update_summary(
      parent_dir = dir,
      split_at = 120,
      split_method = split_method
    )
  }

  # Calculate fTMREL and write risk curve plot to BoP outputs folder
  run_find_tmrel_func_quantiles(outputs_dir = dir,
                                analysis = analysis)

  # Run python script to generate 3-panel plot
  base_dir <- dirname(dir)
  shell_path <- "PATH TO SHELL SCRIPT"

  system(paste(shell_path, base_dir, dir, ro_pair))
}
