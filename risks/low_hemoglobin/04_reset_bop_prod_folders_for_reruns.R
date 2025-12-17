# Define the root path containing the parent folders
root_path <- "PARENT FOLDER PATH"

# List all immediate subdirectories (i.e., parent folders)
parent_folders <- list.dirs(root_path, recursive = FALSE, full.names = TRUE)

# Loop through each parent folder
for (folder in parent_folders) {
  
  # 1. Delete the "results_updated" folder if it exists
  results_updated_path <- file.path(folder, "results_updated")
  if (dir.exists(results_updated_path)) {
    unlink(results_updated_path, recursive = TRUE, force = TRUE)
    message(paste("Deleted folder:", results_updated_path))
  }
  
  # 2. Delete all contents within the "results" folder (but not the folder itself)
  results_path <- file.path(folder, "results")
  if (dir.exists(results_path)) {
    results_files <- list.files(results_path, full.names = TRUE)
    unlink(results_files, recursive = TRUE, force = TRUE)
    message(paste("Cleared contents of:", results_path))
  }
  
  # 3. Delete .csv files starting with "low" in the "data" folder
  data_path <- file.path(folder, "data")
  if (dir.exists(data_path)) {
    low_csvs <- list.files(data_path, pattern = "^low.*\\.csv$", full.names = TRUE)
    unlink(low_csvs, force = TRUE)
    if (length(low_csvs) > 0) {
      message(paste("Deleted .csv files in", data_path, ":", paste(basename(low_csvs), collapse = ", ")))
    }
  }
}
