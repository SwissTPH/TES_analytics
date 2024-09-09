########################
# Merge length-polymorphic alleles
# monica.golumbeanu@unibas.ch
# 21.06.2024
########################

merge_alleles = function(allele_vector, bin_size) {
  
  # Build the underlying histogram
  breaks = seq(from = floor(min(allele_vector)) - 0.5, to = (max(allele_vector) + 1), by = 1)
  allele_values = round((breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2)
  hist_alleles = hist(allele_vector, breaks = breaks, plot = FALSE)
  
  counts_by_offset = sapply(1:bin_size, function (x) sum(hist_alleles$counts[seq(from = x, to = length(hist_alleles$counts), by = bin_size)]))
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