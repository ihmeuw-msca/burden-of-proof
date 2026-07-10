# Import updated raw dataset
dat_new <- data.table::fread(
  'IMPORT PATH TO DATASET'
)

dat_new_exclcrx <- dat_new |> # used for new exclcrx analysis
  dplyr::filter(analysis_exclcrx==1)

dat_new_exclcrx_postnataldep <- dat_new_exclcrx |> # used for new exclcrx & postnataldep analysis
  dplyr::filter(postnataldep==1)

data.table::fwrite(dat_new_exclcrx, "OUTPUT PATH TO MAIN DATASET")
data.table::fwrite(dat_new_exclcrx_postnataldep, "OUTPUT PATH TO POSTATAL DEPRESSION DATASET")
