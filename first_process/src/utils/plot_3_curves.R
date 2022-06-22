#
# plot_3_curves.R
#
# Reed Sorensen
# June 2020
#
CURDIR <- "/ihme/homes/jiaweihe/msca/mrbrt/evidence_score_diet"
CODE_PATH <- paste0(CURDIR, "/parallel/")
source(paste0(CODE_PATH, "00_globals.R"))

library(dplyr)
library(mrbrt001, lib.loc = "/ihme/code/mscm/R/packages/")

pairs_without_selectedcovs <- readRDS(paste0(J_DIR, "/pairs_with_no_selectedcovs.RDS"))
output_dir <- paste0(RESULTS_DIR, VERSION_ID, "/04_monospline_pkl_files/")
ro_pairs <- gsub(".pkl", "", list.files(output_dir))
rds_dir <- paste0(RESULTS_DIR, VERSION_ID, "/04_monospline_models/")
old_dir <- "/home/j/temp/jiaweihe/red_meat_paper/predict_draws/"

pdf(paste0(RESULTS_DIR, VERSION_ID, "/04_monospline_pdf1/compare_curves.pdf"))
for (pair in ro_pairs) {
# for (pair in ro_pairs[ro_pairs %in% pairs_without_selectedcovs]) {
  dev <- FALSE
  if (dev) {
    pair <- "redmeat_diabetes"
  }
  cat(pair, "\n")
  
  # try({
    # Sys.sleep(2)
    x <- readRDS(paste0(rds_dir, pair, ".RDS"))
    mod2 <- py_load_object(paste0(output_dir, pair, ".pkl"))
    
    get_knots_ensemble <- function(model) {
      cov_name_tmp <- model$ensemble_cov_model_name
      tmp <- model$sub_models[[1]]
      tmp2 <- tmp$get_cov_model(name = cov_name_tmp)
      tmp2$spline_knots
    }
    
    draws_path <- paste0(
      paste0(RESULTS_DIR, VERSION_ID, "/04_monospline_pdf1/"),
      x$ro_pair, "_y_draws_fe.pkl"
    )

    draws_dat <- py_load_object(draws_path)
    draws_mean <- exp(apply(draws_dat, 1, mean))
    # draws_mean <- exp(apply(draws_dat, 1, function(x) quantile(x, 0.5)))
    pred_lo <- exp(apply(draws_dat, 1, function(x) quantile(x, 0.025)))
    pred_hi <- exp(apply(draws_dat, 1, function(x) quantile(x, 0.975)))
    
    df_pred2 <- x$df_pred2
    
    old_results <- readRDS(paste0(J_DIR, "old_results_fe.RDS"))
    names(old_results) <- sapply(old_results, function(x) x$ro_pair)
    tmp <- old_results[[x$ro_pair]]

    tmp_vec <- c(draws_mean, tmp$pred_fe)
    x_data <- seq(min(tmp$exp_tmp), max(tmp$exp_tmp), length.out = length(draws_mean))
    new_df <- data.frame(exposure=x_data, pred_lo=pred_lo, pred_hi=pred_hi)

    bias_covs <- x$selected_covs[x$selected_covs != "exposure_linear"]

    old_draws_path <- paste0(old_dir, x$ro_pair, "_y_draws_fe.pkl")
    old_draws_dat <- py_load_object(old_draws_path)
    old_draws_mean <- apply(old_draws_dat, 1, mean)
    old_pred_lo <- apply(old_draws_dat, 1, function(x) quantile(x, 0.025))
    old_pred_hi <- apply(old_draws_dat, 1, function(x) quantile(x, 0.975))
    old_x_data <- seq(min(tmp$exp_tmp), max(tmp$exp_tmp), length.out = length(old_draws_mean))
    old_df <- data.frame(exposure=old_x_data, pred_lo=old_pred_lo, pred_hi=old_pred_hi)
    
    if (length(bias_covs) == 0) {
      covlabel <- "[None]"
    }else {
      covlabel <- paste(bias_covs, collapse = ", ")
    }
    
    covlabel1 <- paste0("c(", paste(bias_covs, collapse = ", "), ")")
    covlabel2 <- paste0("[", paste(tmp$x_covs, collapse = ", "), "]")

    min_y <- min(min(pred_lo), min(tmp_vec), min(old_pred_lo))
    max_y <- max(max(pred_hi), max(tmp_vec), max(old_pred_hi))

    with(df_pred2, plot(
      x_data, draws_mean,
      lwd = 3, type = "l", col="blue",
      main = paste0(x$ro_pair),
      ylim = c(min_y, max_y)
    ))
    abline(h = 1, lwd = 2, lty = 2)
    mtext(text = paste0(paste0(covlabel1, "; ", covlabel2)), side = 3, line = 0.4)
    
    lines(tmp$exp_tmp, tmp$pred_fe, col = adjustcolor("red", 0.6), lwd = 2)

    # function for plotting uncertainty intervals
    add_ui <- function(dat, x_var, lo_var, hi_var, color = "darkblue", opacity = 0.1) {
      polygon(
        x = c(dat[, x_var], rev(dat[, x_var])),
        y = c(dat[, lo_var], rev(dat[, hi_var])),
        col = adjustcolor(col = color, alpha.f = opacity), 
        border = FALSE
      )
    }

    add_ui(new_df, 'exposure', 'pred_lo', 'pred_hi')

    add_ui(old_df, 'exposure', 'pred_lo', 'pred_hi', color="firebrick1")

    if (min(pred_hi) < 1){
        legend_pos <- c("bottomleft")
    }else if (max(pred_lo) > 1){
        legend_pos <- c("topleft")
    }
     
    legend(legend_pos, 
           legend = c("New model", "Old model, 20 iter."),
           lwd = 2, 
           col = c("blue", "red"),
           cex = 0.85
           )
  # })

}

dev.off()
