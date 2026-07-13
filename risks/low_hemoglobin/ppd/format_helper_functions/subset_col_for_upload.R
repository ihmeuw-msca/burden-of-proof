##' *************************************************************************************
##' Title: subset_col_for_upload.R
##' Author: Corey Teply (modified by Taylor Noyes)
##' Purpose: Child function for the pre-bundle upload checks
##' Last updated: 2/14/24
##' Notes: 
##' *************************************************************************************

# source libraries --------------------------------------------------------

library(data.table)
library(stringr)

# subset columns in bundle df ---------------------------------------------

subset_columns <- function(input_df) {
  df <- copy(input_df)
  
  cv_columns <- grep("^cv_|^cov_", names(df), value = TRUE)
  
  bun_cols <- c('nid',
                # 'survey_name',
                'ihme_loc_id',
                'year_start',
                'year_end',
                # 'survey_module',
                'file_path',
                'sex_id',
                'age_start',
                'age_end',
                # 'cv_pregnant',
                'sample_size',
                # 'nclust',
                # 'nstrata',
                # 'var',
                'standard_error',
                # 'design_effect',
                'underlying_nid',
                'field_citation_value',
                'source_type',
                'location_name',
                'location_id',
                'smaller_site_unit',
                'site_memo',
                'sex',
                'sex_id',
                'sex_issue',
                'year_issue',
                'age_issue',
                'age_demographer',
                'measure',
                'mean',
                'lower',
                'upper',
                'cases',
                'unit_type',
                'unit_value_as_published',
                'measure_issue',
                'measure_adjustment',
                'uncertainty_type_value',
                'representative_name',
                'recall_type',
                'case_name',
                'case_definition',
                'case_diagnostics',
                'group',
                'specificity',
                'group_review',
                'note_modeler',
                'note_sr',
                'extractor',
                'seq',
                # 'origin_seq',
                # 'origin_id',
                'is_outlier',
                # 'sample_size_description',
                # 'measure_method',
                # 'lower_cutoff',
                # 'upper_cutoff',
                # 'me_type',
                'underlying_field_citation_value',
                'page_num',
                'table_num',
                'urbanicity_type',
                'recall_type_value',
                'sampling_type',
                'response_rate',
                'year_id',
                'age_group_id',
                # 'standard_deviation',
                # 'cv_biased',
                # 'unit_type_value',
                'input_type',
                'effective_sample_size',
                # 'meanstd1',
                # 'standard_deviation_se',
                'uncertainty_type',
                'val',
                # 'mean_outlier',
                'orig_location_id',
                # 'from_2021',
                cv_columns)
  
  df <- df[, bun_cols, with = FALSE]
  
  return(df)
}