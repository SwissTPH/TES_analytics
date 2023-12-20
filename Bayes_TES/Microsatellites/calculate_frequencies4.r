# calculate frequencies of alleles
# input:
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# alleles_definitions is list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound
# output:
# list of length number of loci
# each entry contains a vector with frequencies of each allele (might not sum to 1 if allele definitions do not cover all observed fragment lengths)
# output[[3]] is mean SD of within allele length
# output[[4]] is mode allele

############################ calculate mode of each allele NEW 9/1/2020
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#mode_allele_lengths = sapply(1:nloci, function (j) sapply(1:dim(alleles_definitions_RR[[j]])[1], function (y) getmode(c(alleles0[,((j-1)*maxMOI+1):(j*maxMOI)][recoded0[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
#																						 allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recodedf[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
#																						 allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recoded_additional_neutral[,((j-1)*maxMOI+1):(j*maxMOI)]==y]))))
mode_allele_lengths = sapply(1:nloci, function (j) sapply(1:dim(alleles_definitions_RR[[j]])[1], function (y) getmode(c(alleles0[,((j-1)*maxMOI+1):(j*maxMOI)][recoded0[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
                                                                                                                        allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recodedf[,((j-1)*maxMOI+1):(j*maxMOI)]==y]))))


calculate_frequencies4 = function(genotypedata, alleles_definitions) {
	
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