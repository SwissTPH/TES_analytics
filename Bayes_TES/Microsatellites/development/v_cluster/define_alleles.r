# generate definitions of alleles (ie binning)
# input:
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# locirepeats is vector of length number of loci with type of locus (dinucleotide, trinucleotide etc repeats)
# maxk is a vector of length length number of loci with the maximum number of alleles for each locus 

# output
# list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound


define_alleles = function(genotypedata, locirepeats, maxk) {
	
ids = genotypedata$Sample.ID

locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))

nids = length(ids)
nloci = length(locinames)

#windows(16,8)
#par(mfrow = c(2,ceiling(nloci/2)))
alleles = list()
observed_data = list()
for (j in 1:nloci) {
	locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
	raw_alleles = c(as.matrix(genotypedata[,locicolumns]))
	raw_alleles = raw_alleles[!is.na(raw_alleles)]
	
	if (diff(range(raw_alleles)) < locirepeats[j]) {
		alleles[[j]] = matrix(c(min(raw_alleles)-locirepeats[j]/2,max(raw_alleles)+locirepeats[j]/2,length(raw_alleles)),1,3)
	} else {
	# remove outliers
	#raw_alleles_range = mean(raw_alleles,na.rm=TRUE) + c(-sd(raw_alleles,na.rm=TRUE)*3,sd(raw_alleles,na.rm=TRUE)*3)
	#raw_alleles = raw_alleles[raw_alleles < raw_alleles_range[2] & raw_alleles > raw_alleles_range[1]]

	breaks = seq(from = floor(min(raw_alleles))-0.5, to = (max(raw_alleles)+1), by = 1)
	allele_values = round((breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)]) / 2)
	hist_alleles = hist(raw_alleles, breaks = breaks, plot = FALSE)

	#k = 8 # number of clusters
	#fit = pam(raw_alleles,k)

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
	
	#alleles[[j]] = unique(clusters)

	#hist_by_cluster = sapply(1:length(possible_alleles), function (x) hist(raw_alleles[clusters == x], breaks = breaks,plot=FALSE)$counts)
	colv = rep("white",length(possible_alleles))
	colv[1:length(possible_alleles) %in% unique(clusters)] = rainbow(k)
	#bar = barplot(t(hist_by_cluster), names = floor(allele_values),cex.names =0.5,col = colv, main = paste(locinames[j],"; Number of alleles: ", k,sep=""),ylim = c(0,max(hist_by_cluster)+2))
	#labels = rep("",length(bar))
	#labels[allele_values %in% possible_alleles[unique(clusters)]] = "*"
	#text(bar, colSums(t(hist_by_cluster))+2,labels)

	# find break values (lower and upper)
	#lower_break_value = breaks[which(allele_values %in% possible_alleles[unique(clusters)])]
	#upper_break_value = breaks[which(allele_values %in% possible_alleles[unique(clusters)])+1]
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