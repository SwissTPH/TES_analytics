########################
# Merge length-polymorphic alleles
# monica.golumbeanu@unibas.ch
#
# 21.06.2024
########################

# Define the bin sizes
source("~/GitRepos/STPHrepos/TES_analytics/define_alleles_markers/define_marker_bin_sizes.R")

merge_alleles = function(allele_vector, bin_size) {
  
  # Remove any NAs from the allele vector
  allele_vector = allele_vector[which(!is.na(allele_vector))]
  
  # Build the underlying histogram
  breaks = seq(from = floor(min(allele_vector)) - 0.5, to = (max(allele_vector) + 1), by = 1)
  allele_values = round((breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2)
  hist_alleles = hist(allele_vector, breaks = breaks, plot = FALSE)
  
  counts_by_offset = sapply(1:min(length(hist_alleles$counts), bin_size), function (x) sum(hist_alleles$counts[seq(from = x, to = length(hist_alleles$counts), by = bin_size)]))
  possible_alleles = allele_values[seq(from = which.max(counts_by_offset), to = length(allele_values), by = bin_size)]
  
  if (min(allele_vector) <= (min(possible_alleles) - bin_size/2)) {
    possible_alleles = c(min(possible_alleles - bin_size), possible_alleles)
  }
  if (max(allele_vector) > (max(possible_alleles) + bin_size/2)) {
    possible_alleles = c(possible_alleles,max(possible_alleles + bin_size))
  }
  
  # assign clusters
  clusters = sapply(allele_vector, function (x) which.min(abs(possible_alleles - x)))
  k = length(unique(clusters))

  return(sort(possible_alleles[unique(clusters)]))
}

process_table_merge = function(input_allele_table, marker_name) {
  
  data_alleles_long = data_alleles %>% pivot_longer(-c(PatientID, Day, Site), 
                                                    names_to = "allele_name", 
                                                    values_to = "allele_length")
  data_alleles_long$allele_id = substr(data_alleles_long$allele_name, 1, regexpr("_", data_alleles_long$allele_name)-1)
  data_alleles_long$Drug = str_extract(data_alleles_long$PatientID, "(?<=_)[A-Z]+(?=_)")

  glurp_idx = which(data_alleles_long$allele_id == "glurp")
  data_alleles_long[glurp_idx, "allele_id"] = marker_name
  data_alleles_long[glurp_idx, "allele_name"] = gsub("glurp", marker_name, data_alleles_long$allele_name[glurp_idx])
  
  # Merge alleles based on bin sizes
  list_markers = unique(data_alleles_long$allele_id)
  merged_allele_table = NULL
  for (marker_name in unique(list_markers)) {
    marker_drug_tab = data_alleles_long %>% filter(allele_id == marker_name)
    allele_vector = marker_drug_tab %>% select(allele_length)
    bin_size = marker_bins %>% 
                filter(marker_id == marker_name) %>% 
                select(bin_size)
    
    print(paste("for", marker_name, bin_size$bin_size))
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
  
  # Transform back from long format to wide format
  merged_allele_table$allele_length = merged_allele_table$allele_id = NULL
  output_table = merged_allele_table %>% pivot_wider(names_from = allele_name, values_from = final_allele_length)
  
  return(output_table)
}

