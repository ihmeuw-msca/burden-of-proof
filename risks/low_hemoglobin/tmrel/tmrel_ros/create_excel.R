library(xlsx)

# define tmrel bop dir ----------------------------------------------------

parent_bop_dir <- 'PATH TO PARENT FOLDER'

# load in data ------------------------------------------------------------

bop_comparison_df <- read.csv(
  file = file.path(parent_bop_dir, 'tmrel_comparision.csv')
)

# get analysis tabs -------------------------------------------------------

bop_comparison_df$sheet_name <- ifelse(
  bop_comparison_df$anemia_severity == '',
  bop_comparison_df$analysis_category,
  paste(
    bop_comparison_df$analysis_category,
    bop_comparison_df$anemia_severity,
    sep = '_'
  )
)
sheet_names <- unique(bop_comparison_df$sheet_name)

# create excel ------------------------------------------------------------

wb <- createWorkbook()

for(s in sheet_names) {
  dat <- bop_comparison_df |>
    dplyr::filter(sheet_name == s)
  dat <- merge.data.frame(
    x = dat |> dplyr::filter(bop_run == 'from_bop'),
    y = dat |> dplyr::filter(bop_run == 'shifted_tmrel'),
    by = c('ro_pair', 'analysis_category', 'anemia_severity', 'sheet_name'),
    suffixes = c('.from_bop', '.shifted_tmrel')
  ) |>
    dplyr::mutate(
      tmrel_diff = tmrel.shifted_tmrel - tmrel.from_bop,
      rr_diff = mean_rr_val.shifted_tmrel - mean_rr_val.from_bop,
      ros_diff = ros.shifted_tmrel - ros.from_bop,
      star_diff = star_rating.shifted_tmrel - star_rating.from_bop,
      sheet_name = NULL,
      bop_run.from_bop = NULL,
      bop_run.shifted_tmrel = NULL
    )
  
  excel_sheet <- createSheet(
    wb = wb,
    sheetName = s
  )
  
  addDataFrame(
    x = dat,
    sheet = excel_sheet,
    row.names = FALSE
  )
}

saveWorkbook(
  wb = wb,
  file = file.path(
    parent_bop_dir,
    'tmrel_comparison.xlsx'
  )
)
