##########################
# Process datasets and group alleles based on bin sizes then return each 
# dataset with the newly defined allele lengths
#
# 10.07.2024
# monica.golumbeanu@unibas.ch
##########################
library(stringr)
library(readxl)
library(tidyr)
library(dplyr)

# Load the functions that combine the alleles
source("~/GitRepos/STPHrepos/TES_analytics/define_alleles_markers/combine_alleles.R")

# TO BE ADAPTED BY THE USER:
# Folder with the datasets. Each dataset is named MSP_marker.xlsx
data_folder = "~/GitRepos/STPHrepos/TES_analytics/define_alleles_markers/"
# Folder to save the results
folder_output = "~/genotyping/tests_alleles/"

# Initialization
final_table = NULL
list_data = list.files(data_folder, full.names = TRUE, pattern = ".xlsx")

# Build the long data table with all the alleles
for (f in list_data) {
  # Extracts the marker name based on the data file name, assuming 
  # the naming convention is MSP_marker
  marker_name = sub("MSP_(.*)\\.xlsx", "\\1", basename(f))
  data_alleles = read_excel(f, sheet = "RecurrentInfections")
  
  data_alleles_long = data_alleles %>% pivot_longer(-c(PatientID, Day, Site), 
                                                    names_to = "allele_name", 
                                                    values_to = "allele_length")
  data_alleles_long$allele_id = substr(data_alleles_long$allele_name, 1, regexpr("_", data_alleles_long$allele_name)-1)
  data_alleles_long$Third_marker = marker_name
  glurp_idx = which(data_alleles_long$allele_id == "glurp")
  data_alleles_long[glurp_idx, "allele_id"] = marker_name
  data_alleles_long[glurp_idx, "allele_name"] = gsub("glurp", marker_name, data_alleles_long$allele_name[glurp_idx])
  
  final_table = rbind.data.frame(final_table, data_alleles_long)
}

# Merge alleles for each marker based on corresponding bin sizes
list_markers = unique(final_table$allele_id)
merged_allele_table = NULL
for (marker_name in unique(list_markers)) {
  print(marker_name)
  
  marker_drug_tab = final_table %>% dplyr::filter(allele_id == marker_name)
  allele_vector = marker_drug_tab %>% select(allele_length)
  bin_size = marker_bins %>% 
    dplyr::filter(marker_id == marker_name) %>% 
    select(bin_size)
  
  merged_alleles = merge_alleles(allele_vector = allele_vector$allele_length, 
                                 bin_size = bin_size$bin_size)
  
  marker_drug_tab$final_allele_length = sapply(marker_drug_tab$allele_length, 
                                               function(x) {
                                                 if (is.na(x)) {
                                                   return(NA)
                                                 } else {
                                                   return(merged_alleles[which.min(abs(merged_alleles - x))])
                                                 }
                                               })
  
  merged_allele_table = rbind.data.frame(merged_allele_table, marker_drug_tab)
  
}

# Generate the dataset files
for (third_marker in unique(merged_allele_table$Third_marker)) {
  # Select dataset by marker 
  output_file = paste0(folder_output, "combined_MSP_", third_marker, ".csv")
  merged_allele_table_dataset = merged_allele_table %>% dplyr::filter(Third_marker == third_marker)
  
  # Transform back from long format to wide format
  merged_allele_table_dataset$allele_length = merged_allele_table_dataset$allele_id = NULL
  output_table = merged_allele_table_dataset %>% 
    pivot_wider(names_from = allele_name, 
                values_from = final_allele_length)
  
  write.csv2(output_table, output_file, row.names = FALSE, quote = FALSE)
}




