
# define tmrel bop dir ----------------------------------------------------

parent_bop_dir <- 'PATH TO PARENT FOLDER'

# load in data ------------------------------------------------------------

bop_comparison_df <- read.csv(
  file = file.path(parent_bop_dir, 'tmrel_comparision.csv')
)

bop_list <- readRDS(
  file = file.path(parent_bop_dir, 'tmrel_comparision.rds')
)
