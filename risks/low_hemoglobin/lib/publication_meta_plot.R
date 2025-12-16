
# load in libraries -------------------------------------------------------

library(data.table)
library(ggplot2)
library(stringr)
source("publication_code/publication_mr_brt.R")
source("publication_code/publication_plot_code.R")

# load in data sets -------------------------------------------------------

df <- fread("PATH TO DATASET")

# step 1 - only keep hb values and OR and RR values for testing
#df <- df[biomarker_type_1=="hemoglobin" & (`Parameter Measure`=="Odds ratio (OR)" | `Parameter Measure`=="Relative risk (RR)") & Adjusted==1 & analysis_v2==1]
df <- df[analysis_exclcrx == 1]
author_df <- unique(df[,.(nid, Reference, Publication_Year)])
author_df$y_axis_title <- unlist(lapply(seq_len(nrow(author_df)), function(x){
  author_df$Reference[x] <- str_remove_all(author_df$Reference[x], ",")
  lead_author <- unlist(str_split(author_df$Reference[x], " "))[1]
  return(paste(lead_author, author_df$Publication_Year[x], sep = ", "))
}))

author_df$y_axis_title[author_df$nid==565199] <- "La Verde, 2014"

fwrite(author_df, 'PATH TO AUTHOR DATASET')

bop_df <- fread("PATH TO BOP DATASET")

# run mr brt meta analysis ------------------------------------------------

mr_brt_meta_analysis <- run_mrbrt_publication(input_df = bop_df)
mrbrt_df <- copy(mr_brt_meta_analysis$mr_brt_df)
mrbrt_df <- merge.data.table(
  x = mrbrt_df,
  y = author_df[,.(nid, y_axis_title)],
  by.x = "study_id",
  by.y = "nid",
  all.x = TRUE
)

# plot results ------------------------------------------------------------

bop_forest_plot_w_mrbrt_publication(
  input_df = mrbrt_df,
  main_title = "MR-BRT Analysis - Meta Analyzed Relative Risk - Peripartum Depression",
  num_decimal_points = 2,
  plot_linear_space = T,
  output_directory = 'PATH TO OUTPUT'
)
