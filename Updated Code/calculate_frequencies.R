# Call the recode_allele function
source("recode_alleles.R")

calculate_frequencies3 = function(genotypedata, alleles_definitions) {
  
  ids = genotypedata$Sample.ID
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  nids = length(ids)
  nloci = length(locinames)
  
  frequencies = list()
  
  variability = c()
  
  for (j in 1:nloci) {
    locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
    raw_alleles = c(as.matrix(genotypedata[,locicolumns]))
    raw_alleles = raw_alleles[!is.na(raw_alleles)]
    low = alleles_definitions[[j]][,1]
    high = alleles_definitions[[j]][,2]
    frequencies[[j]] = sapply(1:dim(alleles_definitions[[j]])[1],function (x) sum(raw_alleles > low[x] & raw_alleles <= high[x]))
    meanSD = mean(sapply(1:dim(alleles_definitions[[j]])[1],function (x) sd(raw_alleles[raw_alleles > low[x] & raw_alleles <= high[x]])),na.rm=TRUE)
    if(is.na(meanSD)) {meanSD = 0}
    variability[j] = meanSD
    frequencies[[j]] = frequencies[[j]] / length(raw_alleles)
  }
  freqmatrix = matrix(0,nloci,max(unlist(lapply(frequencies,length))))
  
  for (j in 1:nloci) {
    freqmatrix[j,1:length(frequencies[[j]])] = frequencies[[j]]
  }
  
  ret = list()
  ret[[1]] = unlist(lapply(frequencies,length))
  ret[[2]] = freqmatrix
  ret[[3]] = variability
  ret
}