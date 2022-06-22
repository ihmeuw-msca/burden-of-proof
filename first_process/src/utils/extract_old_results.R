#
# extract_old_results.R
#
# Reed Sorensen
# June 2020
#


library(reticulate)
library(dplyr)
use_condaenv(condaenv="mr_brt_refactor_env", conda="/ihme/code/evidence_score/miniconda3/bin/conda", required = TRUE)

py_cmds <- c(
  "import sys",
  "import os",
  "import dill as pickle",
  "import argparse",
  "import numpy as np",
  "import pandas as pd",
  "sys.path.append(os.path.dirname('/home/j/temp/reed/prog/repos/mr_brt_ihme/refactor/'))",
  "from mrbrt.__init__ import MR_BRT, MR_BeRT",
  "from mrbrt.utils import ratioInit, sampleKnots"
)





for (cmd in py_cmds) py_run_string(cmd)

path1 <- "/home/j/temp/rmbarber/red_meat_paper/diet_model_pipeline_2020_01_23"
dirs1 <- list.dirs(path1)[-c(1:2)]

path2 <- "/home/j/temp/reed/jiawei/red_meat_paper/diet_model_pipeline_2020_06_29_100_iters"
dirs2 <- list.dirs(path2)[-c(1)]

path3 <- "/home/j/temp/reed/jiawei/red_meat_paper/diet_model_pipeline_2020_06_29_200_iters"
dirs3 <- list.dirs(path3)[-c(1)]


get_old_diet_results <- function(dir, verbose = TRUE) {
  
  dev <- FALSE
  if (dev) {
    dir <- dirs1[27]
  }
  if (verbose) cat(dir, "\n")
  
  try({
    path_stage1 <- paste0(dir, "/stage1.pkl")
    
    if (file.exists(path_stage1)) {
      py_run_string(paste0("with open('", path_stage1, "', 'rb') as fopen:  model1 = pickle.load(fopen)"))
      x_covs <- py$model1$ratio_x_covs
      z_covs <- py$model1$ratio_z_covs
    } else {
      x_covs <- z_covs <- NA
    }
    
    path_mod_mono <- paste0(dir, "/ratio_mod_mono.pkl")
    if (file.exists(path_mod_mono)) {
      py_run_string(paste0("with open('", path_mod_mono, "', 'rb') as fopen:  model2 = pickle.load(fopen)"))
      
      knots_tmp = py$model2$mr$spline_list[[1]]$knots
      k0 = knots_tmp[1]
      k1 = knots_tmp[length(knots_tmp)]
      # pred_x_cov_list, pred_z_cov_list, y_samples, y_samples_fe
      py_run_string(paste0("a, b, c, d = model2.mr_predict(domain = [", k0, ", ", k1, "])"))
      # pred = y_samples.mean(axis = 1)
      pred = apply(py$c, 1, mean)
      pred_fe = apply(py$d, 1, mean)
      # exp_tmp = pred_x_cov_list[0]['mat']
      exp_tmp = py$a[[1]][['mat']]
      df <- py$model2$df
    } else {
      knots_tmp <- pred <- exp_tmp <- df <- NA
    }
    
    ro_pair_tmp <- strsplit(dir, "\\/")[[1]]
    ro_pair <- ro_pair_tmp[length(ro_pair_tmp)]
    
    out <- list(
      ro_pair=ro_pair, x_covs=x_covs, z_covs=z_covs, 
      knots_tmp=knots_tmp, pred=pred, pred_fe=pred_fe, exp_tmp=exp_tmp, df=df
    )
    return(out)
    
  })
}

#####

old_results <- lapply(dirs1, get_old_diet_results)
old_results <- old_results[!sapply(old_results, function(x) class(x) == "try-error")]
names(old_results) <- sapply(old_results, function(x) x$ro_pair)
# saveRDS(old_results, "/home/j/temp/reed/misc/old_results_fe.RDS")

old_results <- readRDS("/home/j/temp/reed/misc/old_results_fe.RDS")
tmp1 <- lapply(old_results, function(x) {
  # x <- old_results[[1]] # dev
  list(
    n_rows = nrow(x$df),
    col_names = names(x$df),
    x_covs = x$x_covs
  )
})
names(tmp11) <- names(old_results)

new_dir <- "/ihme/scratch/users/rsoren/evidence_score_diet/v5_test/04_monospline_models/"
new_results <- lapply(list.files(new_dir, full.names = TRUE), readRDS)
tmp2 <- lapply(new_results, function(x) {
  # x <- new_results[[1]] # dev
  list(
    n_rows = nrow(x$df),
    col_names = names(x$df),
    x_covs = x$selected_covs[x$selected_covs != "exposure_linear"]
  )
})
names(tmp2) <- gsub(".RDS", "", list.files(new_dir))

tmp3 <- do.call("rbind", lapply(gsub(".RDS", "", list.files(new_dir)), function(x) {
  # x <- gsub(".RDS", "", list.files(new_dir))[4] # dev
  data.frame(
    pair = x,
    nrows_old = tmp1[[x]]$n_rows,
    nrows_new = tmp2[[x]]$n_rows,
    xcovs_old = paste(tmp1[[x]]$x_covs, collapse = ","),
    xcovs_new = paste(tmp2[[x]]$x_covs, collapse = ","),
    colnames_old = paste(tmp1[[x]]$col_names, collapse = ","),
    colnames_new = paste(tmp2[[x]]$col_names, collapse = ",")
  )
}))
write.csv(tmp3, "/home/j/temp/reed/misc/comparison_with_ryans_data_prep.csv")


# get pairs without selected covs in Ryan's model
tmp <- sapply(old_results, function(x) x$x_covs)
tmp2 <- names(tmp)[sapply(tmp, function(x) length(x) == 0)]
saveRDS(tmp2, "/home/j/temp/reed/misc/pairs_with_no_selectedcovs.RDS")

##### jiawei's results with 100 iterations
old_results2 <- lapply(dirs2, get_old_diet_results)
old_results2 <- old_results2[!sapply(old_results2, function(x) class(x) == "try-error")]
names(old_results2) <- sapply(old_results2, function(x) x$ro_pair)
saveRDS(old_results2, "/home/j/temp/reed/misc/old_results_fe2.RDS")

##### jiawei's results with 200 iterations
old_results3 <- lapply(dirs3, get_old_diet_results)
old_results3 <- old_results3[!sapply(old_results3, function(x) class(x) == "try-error")]
names(old_results3) <- sapply(old_results3, function(x) x$ro_pair)
saveRDS(old_results3, "/home/j/temp/reed/misc/old_results_fe3.RDS")

