# Drives
os <- Sys.info()["sysname"]

j <- if (os == "Linux") "/home/j/" else if (os == "Windows") "J:/"
h <- if (os == "Linux") paste0("/homes/", user, "/") else if (os == "Windows") "H:/"

out_dir <- "FILEPATH"
WORK_DIR <- "FILEPATH"
central_model_output_folder <- "FILEPATH"

library(dplyr)
library(ggplot2)
library(data.table)


df <- fread("FILEPATH")

# FILE COMPARISON GROUPS IF LOW INTAKE GROUP WAS CODED AS ALTERNATIVE #
df[a_0>b_0, b_0 := a_0]
df[a_1>b_1, b_1 := a_1]

df[, index := 1:.N, by = "nid"]

alt_exp_study <- unique(df[a_0!=b_0, ][,.(nid, b_1,b_0)])
alt_exp_study[, b_midpoint := b_0 + (b_1-b_0)/2]

lower <- as.numeric(quantile(alt_exp_study$b_0, 0.85))
upper <- as.numeric(quantile(alt_exp_study$b_midpoint, 0.85))


tmrel_draws <- data.table("tmrel" = runif(1000, min = lower, max = upper))
tmrel_draws[, `:=`(ro = p, draw = paste0("draw_",0:999), lower = lower, upper = upper)]
write.csv(tmrel_draws, "FILEPATH", row.names = F)
