###########################
# Example of merging alleles according to their bin size and plotting diversity
#
# 26.06.2024
# monica.golumbeanu@unibas.ch
###########################

library(readxl)
library(dplyr)
library(tidyr)
library(plotly)
library(ggplot2)
library(tidyverse)
library(ggpubr)

source("~/GitRepos/bayesian_analysis_Uganda/scripts/merge_alleles.R")
source("~/GitRepos/bayesian_analysis_Uganda/scripts/define_marker_bin_sizes.R")

# Function for processing the data for plotting
process_plot_data = function(plot_data_df) {
  # Make sure that the colnames are the ones expected later for plotting
  colnames(plot_data_df) = c("Haplotype", "Amount", "Frequency")
  
  # Remove rows with NA size
  plot_data_df = plot_data_df[!is.na(plot_data_df$Haplotype), c("Haplotype", "Amount", "Frequency")]
  
  # Process the data frame, define labels
  plot_data_df$Haplotype = factor(plot_data_df$Haplotype, levels = sort(plot_data_df$Haplotype))
  plot_data_df$Label = paste0("n = ", plot_data_df$Amount, "\n(", round(plot_data_df$Amount/sum(plot_data_df$Amount), digits = 2)*100,"%)")
  
  plot_data_df = plot_data_df[order(plot_data_df$Amount, decreasing = TRUE), ]
  
  # Create sector labels
  plot_data_df$pct = as.double(round(plot_data_df$Amount/sum(plot_data_df$Amount), 2) * 100)
  plot_data_df[4:nrow(plot_data_df), "Label"] = ""
  # plot_data_df[4:nrow(plot_data_df), "pct"] = ""
  # 
  # plot_data_df$Label[which(plot_data_df$pct <= 4)] = 0 
  # plot_data_df$pct[which(plot_data_df$pct <= 4)] = 0  # Anything less than 4% should be blank
  # plot_data_df$pct = paste0(plot_data_df$pct, "%")
  # plot_data_df$pct[which(plot_data_df$pct == "0%")] = ""
  # plot_data_df$Label[which(plot_data_df$pct == "")] = ""
  
  return(plot_data_df)
}

plot_pie_ggplot = function(data_df, color_marker, marker_name) {
  plot_data_df = process_plot_data(data_df)
  title_marker = paste0(marker_name, "\n(", nrow(plot_data_df), " genotypes)")
  
  df2 <- plot_data_df %>% 
    mutate(csum = rev(cumsum(rev(Amount))), 
           pos = Amount/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Amount/2, pos))
  
  # plot_data_df <- plot_data_df[order(plot_data_df$Amount), ]
  
  pie_chart = ggplot(plot_data_df, aes(x = "", y = Amount, fill = fct_inorder(as.character(Amount)))) +
    geom_col(width = 1, color = "black") +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "Group")) +
    scale_fill_manual(values = rep(color_marker, 100)) +
    scale_y_continuous(breaks = df2$pos, labels = plot_data_df$Label) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 17),  
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white")) +
    ggtitle(title_marker) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  return(pie_chart)
}

# Read the data
data_alleles = read_excel(f, sheet = "RecurrentInfections")

# Transform into long format
data_alleles_long = data_alleles %>% pivot_longer(-c(PatientID, Day, Site), 
                                                  names_to = "allele_name", 
                                                  values_to = "allele_length")
data_alleles_long$allele_id = substr(data_alleles_long$allele_name, 1, regexpr("_", data_alleles_long$allele_name)-1)
data_alleles_long = data_alleles_long %>% filter(!is.na(allele_length))
data_alleles_long$Day = NULL
data_alleles_long = distinct(data_alleles_long)

# Final table with unmerged alleles
final_table = distinct(data_alleles_long)
list_markers = unique(final_table$allele_id)

# Merge alleles based on bin sizes
merged_allele_table = NULL
for (marker_name in unique(list_markers)) {
  marker_drug_tab = final_table %>% 
    filter(allele_id == marker_name)
  allele_vector = marker_drug_tab %>%
    select(allele_length)
  bin_size = marker_bins %>% 
    filter(marker_id == marker_name) %>% 
    select(bin_size)
  merged_alleles = merge_alleles(allele_vector = allele_vector$allele_length, 
                                 bin_size = bin_size$bin_size)
  
  marker_drug_tab$true_alleles = sapply(marker_drug_tab$allele_length, 
                                        function(x) merged_alleles[which.min(abs(merged_alleles - x))])
  marker_drug_tab$allele_length = marker_drug_tab$allele_name = NULL
  merged_allele_table = rbind.data.frame(merged_allele_table, marker_drug_tab)
}

merged_allele_table = unique(merged_allele_table)
n_patients = length(unique(merged_allele_table$PatientID))

final_table_plot = merged_allele_table %>% 
  group_by(allele_id, true_alleles) %>% 
  summarize(Amount = n())

final_table_plot$Frequency = final_table_plot$Amount/n_patients

colors_markers = c("#fc8d59", "#fc8d59", "#fc8d59", 
                   "#F692EF","#F692EF", "#74a9cf", 
                   "grey", "#ef3b2c", "#807dba",
                   "#f768a1", "#fd8d3c", "#41b6c4", "#CD5C5C")
list_markers = unique(final_table_plot$allele_id)
list_markers = c("K1", "MAD20", "R033","3D7", "FC27", "GLURP", "NULL",
                 "313", "383", "POLY a", "TA1", "TA109", "PFPK2", "2490")
p_array = vector('list', length(list_markers))

for (i in 1:length(list_markers)) {
  marker_name = list_markers[i]
  print(marker_name)
  final_table_plot1 = final_table_plot %>% filter(allele_id == marker_name) 
  final_table_plot1 = final_table_plot1 %>% ungroup() %>% dplyr::select(true_alleles, Amount, Frequency)
  final_table_plot1 = distinct(final_table_plot1)
  
  if (marker_name != "NULL") {
    p_array[[i]] = plot_pie_ggplot(final_table_plot1, colors_markers[i], marker_name)
  } else {
    p_array[[i]] = NULL
  }
  
}

figure = ggarrange(plotlist = p_array, nrow = 2, ncol = 7)

ggsave("~/pies_diversity.pdf", plot = figure,  width = 20, height = 12)
