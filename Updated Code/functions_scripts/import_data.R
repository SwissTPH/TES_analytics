# This file defines a function to import and clean the genotyping data.
# It takes a file path as input and returns a list with two data frames:
# 'late_failures' and 'additional'.

import_and_clean_data <- function(filepath) {
  
  message("INFO: Importing and cleaning data from: ", filepath)
  
  # Read Late Treatment Failures Data
  late_failures_df <- as.data.frame(read_excel(filepath, sheet = "Late Treatment Failures", skip = 3))
  missing_markers <- c(0, "0", "N/A", "-", "NA")
  late_failures_df[late_failures_df %in% missing_markers] <- NA
  late_failures_df$Sample.ID <- sub("D0$", " Day 0", late_failures_df$Sample.ID)
  late_failures_df$Sample.ID <- sub("D[0-9]+$", " Day Failure", late_failures_df$Sample.ID)
  day0_ids <- unique(unlist(strsplit(late_failures_df$Sample.ID[grepl("Day 0", late_failures_df$Sample.ID)], " Day 0")))
  dayF_ids <- unique(unlist(strsplit(late_failures_df$Sample.ID[grepl("Day Failure", late_failures_df$Sample.ID)], " Day Failure")))
  
  if (sum(!paste(day0_ids, "Day Failure") %in% late_failures_df$Sample.ID) > 0) {
    stop("FATAL ERROR: Some 'Day 0' samples are missing their 'Day Failure' pair. Please check your data.")
  }
  if (sum(!paste(dayF_ids, "Day 0") %in% late_failures_df$Sample.ID) > 0) {
    stop("FATAL ERROR: Some 'Day Failure' samples are missing their 'Day 0' pair. Please check your data.")
  }

  # Read Additional Data
  additional_df <- as.data.frame(read_excel(filepath, sheet = "Additional", skip = 3))
  
  if (nrow(additional_df) > 0) {
    additional_df[additional_df %in% missing_markers] <- NA
    additional_df$Sample.ID <- sub("_D0", " Day 0", additional_df$Sample.ID)
    additional_df$Sample.ID <- sub("_D[0-9]*", " Day Failure", additional_df$Sample.ID)
  }

  cols_to_convert_late <- colnames(late_failures_df)[-c(1, 2)]
  late_failures_df[, cols_to_convert_late] <- lapply(late_failures_df[, cols_to_convert_late], as.numeric)
  
  if (nrow(additional_df) > 0) {
    cols_to_convert_add <- colnames(additional_df)[-c(1, 2)]
    additional_df[, cols_to_convert_add] <- lapply(additional_df[, cols_to_convert_add], as.numeric)
  }
  return(
    list(
      late_failures = late_failures_df,
      additional = additional_df
    )
  )
}