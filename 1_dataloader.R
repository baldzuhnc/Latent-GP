if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(tidyverse)
p_load(vroom)

#' Process Panel Data CSV to RDS
#' Reads a fixed CSV path, renames unique_id to pid, and saves to country subfolder.
process_data <- function(file_path, country_code) {
  if (!file.exists(file_path)) {
    message(sprintf("File not found: %s", file_path))
    return(NULL)
  }

  message(sprintf("Processing %s data...", country_code))
  df <- vroom(file_path, guess_max = Inf)

  if ("unique_id" %in% names(df)) {
    df <- df %>% rename(pid = unique_id)
  }

  # Detect wave and ensure output directory
  wv <- max(df$wave, na.rm = TRUE)
  out_dir <- file.path("data", country_code)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_file <- file.path(out_dir, sprintf("W%s_%s.rds", wv, country_code))

  if (!file.exists(out_file)) {
    message(sprintf("Saving standardized RDS: %s", out_file))
    saveRDS(df, out_file)
  } else {
    message(sprintf("Metadata: RDS for wave %s already exists.", wv))
  }
}

# Process the standard files
process_data("data/overall_extended_germany_data.csv", "de")
process_data("data/overall_extended_us_data.csv", "us")
