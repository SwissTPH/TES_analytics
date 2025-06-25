
define_alleles = function(genotypedata, locirepeats, maxk) {
  
  ids = genotypedata$Sample.ID
  
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  
  nids = length(ids)
  nloci = length(locinames)
  
  alleles = list()
  observed_data = list()
  for (j in 1:nloci) {
    locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
    raw_alleles = c(as.matrix(genotypedata[,locicolumns]))
    raw_alleles = raw_alleles[!is.na(raw_alleles)]
    
    if (diff(range(raw_alleles)) < locirepeats[j]) {
      alleles[[j]] = matrix(c(min(raw_alleles)-locirepeats[j]/2,max(raw_alleles)+locirepeats[j]/2,length(raw_alleles)),1,3)
    } else {
      
      
      breaks = seq(from = floor(min(raw_alleles))-0.5, to = (max(raw_alleles)+1), by = 1)
      allele_values = round((breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)]) / 2)
      hist_alleles = hist(raw_alleles, breaks = breaks, plot = FALSE)
      
      
      counts_by_offset = sapply(1:locirepeats[j], function (x) sum(hist_alleles$counts[seq(from = x, to = length(hist_alleles$counts), by = locirepeats[j])]))
      possible_alleles = allele_values[seq(from = which.max(counts_by_offset), to = length(allele_values), by = locirepeats[j])]
      
      if (min(raw_alleles) <= (min(possible_alleles)-locirepeats[j]/2)) {
        possible_alleles = c(min(possible_alleles-locirepeats[j]),possible_alleles)
      }
      if (max(raw_alleles) > (max(possible_alleles)+locirepeats[j]/2)) {
        possible_alleles = c(possible_alleles,max(possible_alleles+locirepeats[j]))
      }
      
      # assign clusters
      clusters = sapply(raw_alleles, function (x) which.min(abs(possible_alleles - x)))
      k = length(unique(clusters))
      
      colv = rep("white",length(possible_alleles))
      colv[1:length(possible_alleles) %in% unique(clusters)] = rainbow(k)
      
      lower_break_value = sort(possible_alleles[unique(clusters)] - locirepeats[j]/2)
      upper_break_value = sort(possible_alleles[unique(clusters)] + locirepeats[j]/2)
      counts = sapply(1:length(lower_break_value), function (x) sum(raw_alleles > lower_break_value[x] & raw_alleles <= upper_break_value[x]))
      alleles[[j]] = cbind(lower_break_value, upper_break_value, counts)
    }
  }
  
  #### compress
  # take maxk most frequent alleles
  
  alleles2 = list()
  for (j in 1:nloci) {
    sortedindex = sort.int(alleles[[j]][,3],decreasing = TRUE,index.return = TRUE)$ix[1:maxk[j]]
    if (length(alleles[[j]][,3]) <= maxk[j]) {
      sortedindex = sort.int(alleles[[j]][,3],decreasing = TRUE,index.return = TRUE)$ix
    }
    print(sum(alleles[[j]][sortedindex,3])/sum(alleles[[j]][,3]))
    alleles2[[j]] = cbind(alleles[[j]][sortedindex,1],alleles[[j]][sortedindex,2])
  }
  alleles2
}