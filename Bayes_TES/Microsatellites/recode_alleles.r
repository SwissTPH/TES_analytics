# recode alleles
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# the first column has to have a unique ID followed by either " Day 0" or " Day Failure"
# alleles_definitions is list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound

# output
# [[1]] is list of length number of loci
# each entry in the length is a vector of length number of ids
# each entry in the vector is a string with the following format:
# A-B-C/D-E
# the letters represent the Ath, Bth, Cth etc allele (as defined in alleles_definitions)
# the letters before the "/" represent day 0 alleles, and after day of failure alleles
# Example: 5-4/2 represents an individual that had alleles 5 and 4 on day 0, and then allele 2 on day of failure
# when nothing appears, either no alleles were detected, or alleles fell outside ranges in alleles_definitions 
# row names are the ids
# [[2]] is a number of ids X 2 matrix of multiplicity of infection, where first column is day 0 MOI and second column is day of failure MOI

recodeallele = function(alleles_definitions_subset,proposed) {
	
	ret = which(proposed > alleles_definitions_subset[,1] & proposed <= alleles_definitions_subset[,2])
	if (length(ret) == 0) {
		ret = NA
	}
	ret
}

recode_alleles = function(genotypedata, alleles_definitions) {

########### generate MOI for each sample

ids = unique(unlist(strsplit(genotypedata$Sample.ID[grepl("Day 0",genotypedata$Sample.ID)]," Day 0")))
locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
nids = length(ids)
nloci = length(locinames)


MOI0 = rep(0,nids)
MOIf = rep(0,nids)

# for each individual,
# cycle through each locus and count number of separate alleles

for (i in 1:nids) {
	for (j in 1:nloci) {
		locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
		nalleles0 = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day 0"),genotypedata$Sample.ID),locicolumns]))
		nallelesf = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day Failure"),genotypedata$Sample.ID),locicolumns]))

		MOI0[i] = max(MOI0[i],nalleles0)
		MOIf[i] = max(MOIf[i],nallelesf)
	}
}



observeddatamatrix = list()
for (j in 1:nloci) {
	locus = locinames[j]
	locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata))
	oldalleles = as.vector(genotypedata[,locicolumns])
	if (length(dim(oldalleles)[2]) == 0) {
		oldalleles = matrix(oldalleles,length(oldalleles),1)
	}
	newalleles = oldalleles
	ncolumns = dim(oldalleles)[2]
	for (i in 1:ncolumns) {
		newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions[[j]],oldalleles[x,i])))
	}
	newalleles[is.na(newalleles)] = ""
	
	tempobservedata = c()
	for (i in 1:nids) {
		locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
		day0alleles = newalleles[grepl(paste(ids[i],"Day 0"),genotypedata$Sample.ID),]
		day0alleles = day0alleles[day0alleles != ""]
		dayfalleles = newalleles[grepl(paste(ids[i],"Day Failure"),genotypedata$Sample.ID),]
		dayfalleles = dayfalleles[dayfalleles != ""]
		tempobservedata[i] = paste(paste(sort(unique(as.numeric(day0alleles))),collapse="-"),paste(sort(unique(as.numeric(dayfalleles))),collapse="-"),sep="/")
	}
	observeddatamatrix[[j]] = tempobservedata
}
MOItemp = cbind(MOI0,MOIf)
rownames(MOItemp) = ids
list(observeddatamatrix = observeddatamatrix, MOI = MOItemp)
}
