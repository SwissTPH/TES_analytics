###########################
# Script for classifying recrudescence and new infection with the WHO 3/3 algorithm
# using msp1, msp2 and glurp
#
# According to the algorithm, first a classification is done at the allelic famility level.
# There are three gene markers: msp1, msp2 and glurp/microsatellite.
# For each marker, a recrudescence is called if we observe a recrudescence for at least one allelic family.
# The sample is classified as a recrudescence if all three markers are classified in part as recrudescence.
#
# allele naming: markerName_allelicFamily_alleleNumber
# sample naming: sampleName dayNumber
# Sample name should have two rows and these will correspond to day 0 and day X
# No check is done on the day numbers
#
# monica.golumbeanu@swisstph.ch
###########################
library(readxl)
library(dplyr)
library(numbers)


################## USER INPUT #############################
# Specify the path to your input and output data files
input_data_file = "~/genotyping/checking_samples/Input_CKAF156A2202_msp1_msp2_glurp_results.xlsx"
output_data_file = "~/genotyping/checking_samples/Output_CKAF156A2202_msp1_msp2_glurp_results.xlsx"

# Specify bin sizes for the markers
# Allele lengths within bin sizes are considered the same
bin_sizes = c(3, 3, 3)

###########################################################

# Load the data
input_data = as.data.frame(read_xlsx(input_data_file))

# Sample id name
s_id_name = colnames(input_data)[1]

# Data parsing and checking
cols = colnames(input_data)
l_cols = length(cols)
cols = cols[2:l_cols]
input_data[, 2:l_cols] = as.data.frame(lapply(input_data[, 2:l_cols], function(x) as.double(x)))
# Create a list with all the markers and their allelic families
dict = NULL
dict$marker = sapply(strsplit(cols, "_"), `[`, 1)
dict$af = sapply(strsplit(cols, "_"), `[`, 2)
dict = as.data.frame(dict)
dict$bin_size = 0

marker_list = unique(dict$marker)

#check if bin sizes are specified for all markers
if (length(bin_sizes) < length(marker_list)) {
  stop(pste0("Not all markers have bin sizes specified! Bins: ", length(bin_sizes), ", markers: ", length(markers)))
} else {
  # Assign bin sizes
  for (i in 1:length(marker_list)) {
    m = marker_list[i]
    dict[which(dict$marker == m), "bin_size"] = bin_sizes[i]
  }
}

# Print correspondence between markers, allelic families and bin sizes
print("Identified correspondences:")
print(dict)

# Final results table
results_tab = data.frame(matrix(ncol = length(marker_list) + 1, nrow = nrow(input_data)/2, 
                                dimnames = list(NULL, c(s_id_name, marker_list))))

# Marker classification
for (i in seq(1, nrow(input_data), 2)) {
  # Check if the sample names are the same
  s1 = strsplit(input_data[i, s_id_name], " ")[[1]][1]
  s2 = strsplit(input_data[i+1, s_id_name], " ")[[1]][1]
  if (s1 != s2) {
    stop(paste0("Error processing the samples. Day 0 and day x entries not found for sample ", input_data[i, s_id_name]))
  }
  
  # Initialization
  results_tab[div(i, 2) + 1, s_id_name] = s1 
  
  print(s1)
  
  # Loop through markers and classify them
  for (m in unique(dict$marker)) {
    # Initialization
    results_tab[div(i, 2) + 1, m] = "NI"
    # Extract the allelic families for the selected markers
    afs_m = unlist(unique(dict %>% filter(marker == m) %>% select(af)))
    # Check all the pairwise distances between alleles
    for (j in 1:length(afs_m)) {
      af_m = afs_m[j]
      print(af_m)
      v1 = input_data[i, which(grepl(paste0(m, "_", af_m), colnames(input_data)))]
      v2 = input_data[i+1, which(grepl(paste0(m, "_", af_m), colnames(input_data)))]
      all_combinations = as.data.frame(expand.grid(c(v1), c(v2)))
      all_combinations$diff = abs(unlist(all_combinations$Var1) - unlist(all_combinations$Var2))
      if (any(!is.na(all_combinations$diff) & 
              all_combinations$diff <= unique(dict[which(dict$marker == m & dict$af == af_m), "bin_size"]))) {
        results_tab[div(i, 2) + 1, m] = "R"
      }
      
      # For debugging
      if(input_data[i, s_id_name] == "1021036 D0") {
        print("a")
      }
    }
  }
}
  
  # Sample classification
  results_tab$Recrudescence_WHO = (rowSums(results_tab[,2:4]=="R")==3)
  results_tab[which(results_tab$Recrudescence_WHO), "Recrudescence_WHO"] = "R"
  results_tab[which(results_tab$Recrudescence_WHO == FALSE), "Recrudescence_WHO"] = "NI"
  
  results_tab$`Recrudescence_2/3` = (rowSums(results_tab[,2:4]=="R")>=2)
  results_tab[which(results_tab$`Recrudescence_2/3`), "Recrudescence_2/3"] = "R"
  results_tab[which(results_tab$`Recrudescence_2/3` == FALSE), "Recrudescence_2/3"] = "NI"
  
  # Saving final results to file
  write.csv(results_tab, output_data_file)
  
  
  
  