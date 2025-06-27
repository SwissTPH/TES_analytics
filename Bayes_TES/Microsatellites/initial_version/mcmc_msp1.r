

ids = unique(unlist(strsplit(genotypedata_RR$Sample.ID[grepl( "Day 0",genotypedata_RR$Sample.ID)]," Day 0")))

arbitrarydistance = 1000
thresholddistance = 500

## recode msp1 and msp2 alleles

returnnonempty = function(x) {
	x[!is.na(x)]
}

msp1_mad20 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("MAD20",colnames(genotypedata_RR))]))
msp1_K1 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("K1",colnames(genotypedata_RR))])+arbitrarydistance)
#msp1_RO33 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("R033",colnames(genotypedata_RR))])+arbitrarydistance * 2+rnorm(1,0,sd=4))
msp1_RO33 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("R033",colnames(genotypedata_RR))])+arbitrarydistance * 2)
msp1_all = sapply(1:(dim(genotypedata_RR)[1]), function (x) c(msp1_mad20[[x]],msp1_K1[[x]],msp1_RO33[[x]]))

msp1_MOI = unlist(lapply(msp1_all,length))
temp = matrix(NA,(dim(genotypedata_RR)[1]),max(msp1_MOI))
sapply(which(msp1_MOI!=0), function (x) temp[x,1:(msp1_MOI[x])] <<-  msp1_all[[x]][1:(msp1_MOI[x])])
msp1=temp
colnames(msp1) = paste("MSP1_",1:max(msp1_MOI),sep="")

msp2_3D7 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("3D7",colnames(genotypedata_RR))]))
msp2_FC27 = sapply(1:(dim(genotypedata_RR)[1]), function (x) returnnonempty(genotypedata_RR[x,grep("FC27",colnames(genotypedata_RR))])+arbitrarydistance)
msp2_all = sapply(1:(dim(genotypedata_RR)[1]), function (x) c(msp2_3D7[[x]],msp2_FC27[[x]]))

msp2_MOI = unlist(lapply(msp2_all,length))
temp = matrix(NA,(dim(genotypedata_RR)[1]),max(msp2_MOI))
sapply(which(msp2_MOI!=0), function (x) temp[x,1:(msp2_MOI[x])] <<-  msp2_all[[x]][1:(msp2_MOI[x])])
msp2=temp
colnames(msp2) = paste("MSP2_",1:max(msp2_MOI),sep="")

colnames(genotypedata_RR) = gsub("glurp.","glurp_",colnames(genotypedata_RR))

temp = cbind(Sample.ID=genotypedata_RR$Sample.ID,msp1,msp2,genotypedata_RR[,grep("glurp_",colnames(genotypedata_RR))])

##### NORMALIZE to make sure that all 3 markers have similar mean length (so error rate applies equally to all three)
msp1_mean = mean(do.call(c,temp[,grepl("MSP1",colnames(temp))])%% 1000,na.rm=TRUE)
msp2_mean = mean(do.call(c,temp[,grepl("MSP2",colnames(temp))])%% 1000,na.rm=TRUE)
glurp_mean = mean(do.call(c,temp[,grepl("glurp",colnames(temp))])%% 1000,na.rm=TRUE)
temp[,grepl("MSP2",colnames(temp))] = round((temp[,grepl("MSP2",colnames(temp))]%% 1000)*msp1_mean/msp2_mean)+1000*(temp[,grepl("MSP2",colnames(temp))]%/% 1000)
temp[,grepl("glurp",colnames(temp))] = round((temp[,grepl("glurp",colnames(temp))])*msp1_mean/glurp_mean)

genotypedata_RR=temp


msp1_mad20 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("MAD20",colnames(additional_neutral))]))
msp1_K1 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("K1",colnames(additional_neutral))])+arbitrarydistance)
#msp1_RO33 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("R033",colnames(additional_neutral))])+arbitrarydistance * 2+rnorm(1,0,sd=4))
msp1_RO33 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("R033",colnames(additional_neutral))])+arbitrarydistance * 2)
msp1_all = sapply(1:(dim(additional_neutral)[1]), function (x) c(msp1_mad20[[x]],msp1_K1[[x]],msp1_RO33[[x]]))

msp1_MOI = unlist(lapply(msp1_all,length))
temp = matrix(NA,(dim(additional_neutral)[1]),max(msp1_MOI))
sapply(which(msp1_MOI!=0), function (x) temp[x,1:(msp1_MOI[x])] <<-  msp1_all[[x]][1:(msp1_MOI[x])])
msp1=temp
colnames(msp1) = paste("MSP1_",1:max(msp1_MOI),sep="")

msp2_3D7 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("3D7",colnames(additional_neutral))]))
msp2_FC27 = sapply(1:(dim(additional_neutral)[1]), function (x) returnnonempty(additional_neutral[x,grep("FC27",colnames(additional_neutral))])+arbitrarydistance)
msp2_all = sapply(1:(dim(additional_neutral)[1]), function (x) c(msp2_3D7[[x]],msp2_FC27[[x]]))

msp2_MOI = unlist(lapply(msp2_all,length))
temp = matrix(NA,(dim(additional_neutral)[1]),max(msp2_MOI))
sapply(which(msp2_MOI!=0), function (x) temp[x,1:(msp2_MOI[x])] <<-  msp2_all[[x]][1:(msp2_MOI[x])])
msp2=temp
colnames(msp2) = paste("MSP2_",1:max(msp2_MOI),sep="")

colnames(additional_neutral) = gsub("glurp.","glurp_",colnames(additional_neutral))

temp = cbind(Sample.ID=additional_neutral$Sample.ID,msp1,msp2,additional_neutral[,grep("glurp_",colnames(additional_neutral))])

##### NORMALIZE to make sure that all 3 markers have similar mean length (so error rate applies equally to all three)
msp1_mean = mean(do.call(c,temp[,grepl("MSP1",colnames(temp))])%% 1000,na.rm=TRUE)
msp2_mean = mean(do.call(c,temp[,grepl("MSP2",colnames(temp))])%% 1000,na.rm=TRUE)
glurp_mean = mean(do.call(c,temp[,grepl("glurp",colnames(temp))])%% 1000,na.rm=TRUE)
temp[,grepl("MSP2",colnames(temp))] = round((temp[,grepl("MSP2",colnames(temp))]%% 1000)*msp1_mean/msp2_mean)+1000*(temp[,grepl("MSP2",colnames(temp))]%/% 1000)
temp[,grepl("glurp",colnames(temp))] = round((temp[,grepl("glurp",colnames(temp))])*msp1_mean/glurp_mean)

additional_neutral=temp


locinames = unique(sapply(colnames(genotypedata_RR)[-1],function(x) strsplit(x,"_")[[1]][1]))
nloci = length(locinames)
nids = length(ids)
maxalleles=100
k = rep(maxalleles, nloci)
alleles_definitions_RR  = define_alleles(smartbind(genotypedata_RR,additional_neutral),rep(locirepeats[1],nloci),k) # make sure to use same locirepeats if normalizing



maxMOI = max(as.numeric(sapply(1:length(colnames(genotypedata_RR)), function (x) strsplit(colnames(genotypedata_RR)[x],"_")[[1]][2])),na.rm=TRUE)


##### calculate MOI
MOI0 = rep(0,nids)
MOIf = rep(0,nids)
for (i in 1:nids) {
	for (j in 1:nloci) {
		locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata_RR))
		nalleles0 = sum(!is.na(genotypedata_RR[grepl(paste("^",ids[i]," Day 0",sep=""),genotypedata_RR$Sample.ID),locicolumns]))
		nallelesf = sum(!is.na(genotypedata_RR[grepl(paste("^",ids[i]," Day Failure",sep=""),genotypedata_RR$Sample.ID),locicolumns]))

		MOI0[i] = max(MOI0[i],nalleles0)
		MOIf[i] = max(MOIf[i],nallelesf)
	}
}
#maxMOI = max(MOI0,MOIf)


##### define statevector

alleles0 = matrix(0,nids,maxMOI*nloci)
recoded0 = matrix(0,nids,maxMOI*nloci)
hidden0 = matrix(NA,nids,maxMOI*nloci)
hidden_crossfamily0 = matrix(NA,nids,maxMOI*nloci)
recr0 = matrix(NA,nids,nloci)
recr_repeats0 = matrix(NA,nids,nloci) # number of times recrudescing allele is repeated on day 0
recr_repeatsf = matrix(NA,nids,nloci) # number of times recrudescing allele is repeated on day 0
allelesf = matrix(0,nids,maxMOI*nloci)
recodedf = matrix(0,nids,maxMOI*nloci)
hiddenf = matrix(NA,nids,maxMOI*nloci)
hidden_crossfamilyf = matrix(NA,nids,maxMOI*nloci)
recrf = matrix(NA,nids,nloci)
if (length(additional_neutral) > 0) { if (dim(additional_neutral)[1] > 0) {
	recoded_additional_neutral = matrix(0,dim(additional_neutral)[1],maxMOI*nloci)
}}
mindistance = matrix(0,nids,nloci)
alldistance = array(NA,c(nids,nloci,maxMOI*maxMOI))
allrecrf = array(NA,c(nids,nloci,maxMOI*maxMOI))
classification = rep(0,nids)
##### create state 0

for (j in 1:nloci) {
	locus = locinames[j]
	locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata_RR))
	oldalleles = as.vector(genotypedata_RR[,locicolumns])
	if (length(dim(oldalleles)[2]) == 0) {
		oldalleles = matrix(oldalleles,length(oldalleles),1)
	}
	newalleles = oldalleles
	ncolumns = dim(oldalleles)[2]
	for (i in 1:ncolumns) {
		newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions_RR[[j]],oldalleles[x,i])))
	}
	newalleles = matrix(as.numeric(unlist(c(newalleles))),dim(newalleles)[1],dim(newalleles)[2])
	newalleles[is.na(newalleles)] = 0
	oldalleles = matrix(as.numeric(unlist(c(oldalleles))),dim(oldalleles)[1],dim(oldalleles)[2])
	oldalleles[is.na(oldalleles)] = 0

	oldalleles[newalleles == 0] = 0
	alleles0[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day 0",genotypedata_RR$Sample.ID),]
	allelesf[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day Failure",genotypedata_RR$Sample.ID),]
	recoded0[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(newalleles)[2])] = newalleles[grepl("Day 0",genotypedata_RR$Sample.ID),]
	recodedf[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(newalleles)[2])] = newalleles[grepl("Day Failure",genotypedata_RR$Sample.ID),]

}


# 
# if (length(additional_neutral) > 0) { if (dim(additional_neutral)[1] > 0) {
# recoded_additional_neutral = matrix(0,dim(additional_neutral)[1],maxMOI*nloci)
# ##### recode additional_neutral
# for (j in 1:nloci) {
# 	locus = locinames[j]
# 	locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata_RR))
# 	oldalleles = as.vector(additional_neutral[,locicolumns])
# 	if (length(dim(oldalleles)[2]) == 0) {
# 		oldalleles = matrix(oldalleles,length(oldalleles),1)
# 	}
# 	newalleles = oldalleles
# 	ncolumns = dim(oldalleles)[2]
# 	for (i in 1:ncolumns) {
# 		newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions_RR[[j]],oldalleles[x,i])))
# 	}
# 	newalleles = matrix(as.numeric(unlist(c(newalleles))),dim(newalleles)[1],dim(newalleles)[2])
# 	newalleles[is.na(newalleles)] = 0
# 	oldalleles = matrix(as.numeric(unlist(c(oldalleles))),dim(oldalleles)[1],dim(oldalleles)[2])
# 	oldalleles[is.na(oldalleles)] = 0
# 
# 	oldalleles[newalleles == 0] = 0
# 	recoded_additional_neutral[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = newalleles
# }
# } else {
# 	recoded_additional_neutral = c()
# }}

## estimate frequencies

frequencies_RR = calculate_frequencies3(smartbind(genotypedata_RR,additional_neutral),alleles_definitions_RR)


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#mode_allele_lengths = sapply(1:nloci, function (j) sapply(1:dim(alleles_definitions_RR[[j]])[1], function (y) getmode(c(alleles0[,((j-1)*maxMOI+1):(j*maxMOI)][recoded0[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
#																						 allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recodedf[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
#																						 allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recoded_additional_neutral[,((j-1)*maxMOI+1):(j*maxMOI)]==y]))))
mode_allele_lengths = sapply(1:nloci, function (j) sapply(1:dim(alleles_definitions_RR[[j]])[1], function (y) getmode(c(alleles0[,((j-1)*maxMOI+1):(j*maxMOI)][recoded0[,((j-1)*maxMOI+1):(j*maxMOI)]==y],
                                                                                                                        allelesf[,((j-1)*maxMOI+1):(j*maxMOI)][recodedf[,((j-1)*maxMOI+1):(j*maxMOI)]==y]))))

mode_allele_lengths_new = mode_allele_lengths
### for empty ones use median 
for (i in 1:nloci) {
  
  whichempty = which(is.na(mode_allele_lengths[[i]]))
  if (length(whichempty)>0) {
    for (j in 1:length(whichempty)){
      mode_allele_lengths_new[[i]][whichempty[j]]= mean(alleles_definitions_RR[[i]][whichempty[j],])
    }
  }
}
mode_allele_lengths=mode_allele_lengths_new



## assign random hidden alleles
for (i in 1:nids) {
	for (j in 1:nloci) {
		nalleles0 = sum(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0)
		nmissing0 = MOI0[i] - nalleles0
		whichnotmissing0 = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] != 0)]
		whichmissing0 = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] == 0)]

		if (nalleles0 > 0) {
			hidden0[i,whichnotmissing0] = 0
		}

		if (nmissing0 == MOI0[i]) { 
			possible_hidden_index0 = 1:(frequencies_RR[[1]][j])
		} else {
			currentalleles = alleles0[i,whichnotmissing0]
			pairwisecomp = expand.grid(currentalleles,mode_allele_lengths[[j]])
			possible_hidden_index0 = which(mode_allele_lengths[[j]] %in% pairwisecomp[which(abs(pairwisecomp[,1]-pairwisecomp[,2])<500),2])
		}
		if (nmissing0 > 0) {
			if (length(possible_hidden_index0)>1) {
			newhiddenalleles0 = sample(possible_hidden_index0,nmissing0,replace=TRUE,frequencies_RR[[2]][j,possible_hidden_index0])
			} else {
			  print(as.factor(possible_hidden_index0))
			newhiddenalleles0 = as.numeric(as.character(sample(as.factor(possible_hidden_index0),nmissing0,replace=TRUE,frequencies_RR[[2]][j,possible_hidden_index0])))
			}
			jitter0 = rnorm(nmissing0,mean=0,sd=frequencies_RR[[3]][j])
			jitter0[jitter0 > locirepeats[j]/2] = locirepeats[j]/2
			jitter0[jitter0 < -locirepeats[j]/2] = -locirepeats[j]/2
			recoded0[i,whichmissing0] = newhiddenalleles0
			alleles0[i,whichmissing0] = mode_allele_lengths[[j]][newhiddenalleles0] + jitter0# hidden alleles get mean allele length + jitter
			hidden0[i,whichmissing0] = 1
			hidden_crossfamily0[i,whichmissing0] = 0
			#print(alleles0[i,whichnotmissing0])
			#print(alleles0[i,whichmissing0])
		}
		nallelesf = sum(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0)
		nmissingf = MOIf[i] - nallelesf
		whichnotmissingf = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] != 0)]
		whichmissingf = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] == 0)]

		if (nallelesf > 0) {
			hiddenf[i,whichnotmissingf] = 0
		}
		if (nmissingf == MOIf[i]) { 
			possible_hidden_indexf = 1:(frequencies_RR[[1]][j])
		} else {
			currentalleles = allelesf[i,whichnotmissingf]
			pairwisecomp = expand.grid(currentalleles,mode_allele_lengths[[j]])
			possible_hidden_indexf = which(mode_allele_lengths[[j]] %in% pairwisecomp[which(abs(pairwisecomp[,1]-pairwisecomp[,2])<500),2])
		}

		if (nmissingf > 0) {
			if (length(possible_hidden_indexf)>1) {
			newhiddenallelesf = sample(possible_hidden_indexf,nmissingf,replace=TRUE,frequencies_RR[[2]][j,possible_hidden_indexf])
			} else {
			newhiddenallelesf = as.numeric(as.character(sample(as.factor(possible_hidden_indexf),nmissingf,replace=TRUE,frequencies_RR[[2]][j,possible_hidden_indexf])))
			}
			jitterf = rnorm(nmissingf,mean=0,sd=frequencies_RR[[3]][j])
			jitterf[jitterf > locirepeats[j]/2] = locirepeats[j]/2
			jitterf[jitterf < -locirepeats[j]/2] = -locirepeats[j]/2

			recodedf[i,whichmissingf] = newhiddenallelesf
			allelesf[i,whichmissingf] = mode_allele_lengths[[j]][newhiddenallelesf] + jitterf# hidden alleles get mean allele length + jitter
			hiddenf[i,whichmissingf] = 1
			hidden_crossfamilyf[i,whichmissingf] = 0

		}
	}
}

## initial estimate of q (probability of an allele being missed)
qq = mean(c(hidden0,hiddenf),na.rm=TRUE)
qq_crossfamily = 10^-3

## initial estimate of dvect (likelihood of error in analysis)
dvect = rep(0,1+round(max(sapply(1:nloci,function (x) diff(range(c(alleles_definitions_RR[[x]])))))))
dvect[1] = 0.75
dvect[2] = 0.2
dvect[3] = 0.05
dvect[-c(1,2,3)] = 10^-5
## randomly assign recrudescences/reinfections
for (i in 1:nids) {
	z = runif(1)
	if (z < 0.5) {
		classification[i] = 1
	}
	for (j in 1:nloci) { # determine which alleles are recrudescing (for beginning, choose closest pair)
		allpossiblerecrud = expand.grid(1:MOI0[i],1:MOIf[i])
		closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (x) abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]])))
		mindistance[i,j] = abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]])
		alldistance[i,j,1:dim(allpossiblerecrud)[1]] = sapply(1:dim(allpossiblerecrud)[1], function (x) abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]]))
		allrecrf[i,j,1:dim(allpossiblerecrud)[1]] = recodedf[i,maxMOI*(j-1)+allpossiblerecrud[,2]]
		recr0[i,j] = maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]
		recrf[i,j] = maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]
		recr_repeats0[i,j] = sum(recoded0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recoded0[i,recr0[i,j]])
		recr_repeatsf[i,j] = sum(recodedf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recodedf[i,recrf[i,j]])
	}
}





#### correction factor (reinfection)
correction_distance_matrix = list() # for each locus, matrix of distances between each allele
correction_distance_matrix_constant = list() # for each locus, matrix of distances between each allele

for (i in 1:nloci) {
	correction_distance_matrix_constant[[i]] = as.matrix(dist((mode_allele_lengths[[i]])))
	### correct for variation in allele lengths within the same bin
	## correction_distance_matrix[[i]] = correction_distance_matrix_constant[[i]] + frequencies_RR[[3]][i]*sqrt(2/pi)
	correction_distance_matrix[[i]] = correction_distance_matrix_constant[[i]] + abs(rnorm(dim(correction_distance_matrix_constant[[i]])[1],mean=0,frequencies_RR[[3]][i]))
}


state_classification = matrix(NA,nids,(nruns-burnin)/record_interval)
state_alleles0 = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
state_allelesf = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
state_parameters = matrix(NA,3+2*nloci,(nruns-burnin)/record_interval)

count = 1
dposterior = 0.75
runmcmc = function() {
#for (mis in 1:10000) {
	# propose new classification
	# rellikelihood_reinfection = sapply(1:nids, function (x) (sum(log(frequencies_RR[[2]][cbind(1:nloci,recoded0[x,recrf[x,]])]))))
	#rellikelihood_recr = sapply(1:nids, function (x) (sum(log(dvect[round(mindistance[x,]+1)]))))
	# likelihoodratio = exp(rellikelihood_recr - rellikelihood_reinfection)
	# adjust for multiple corrections (ratio of multinomial coefficients)
	#likelihoodratio = sapply(1:nids, function (x) likelihoodratio[x]/exp(nloci*log(MOI0[x])+nloci*log(MOIf[x])-sum(log(recr_repeats0[x,]))-sum(log(recr_repeatsf[x,]))))
	#likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/frequencies_RR[[2]][y,allrecrf[x,y,]],na.rm=TRUE))))))
	#likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/colSums(frequencies_RR[[2]][y,1:frequencies_RR[[1]][y]]*matrix(dvect[correction_distance_matrix[[y]][,allrecrf[x,y,]]+1],frequencies_RR[[1]][y],frequencies_RR[[1]][y])),na.rm=TRUE))))))
	likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/
																		sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][y,1:frequencies_RR[[1]][y]]*dvect[round(correction_distance_matrix[[y]][,allrecrf[x,y,z]]+1)])),na.rm=TRUE))))))

	z = runif(nids)
	newclassification = classification
	newclassification[classification == 0 & z < likelihoodratio] = 1
	newclassification[classification == 1 & z < 1/likelihoodratio] = 0
	classification <<- newclassification
	
	# propose new hidden states
	sapply(1:nids, function (x) switch_hidden(x))
	
	# propose q (beta distribution is conjugate distribution for binomial process)
	q_prior_alpha = 0;
	q_prior_beta = 0;
	q_posterior_alpha = q_prior_alpha + sum(c(hidden0,hiddenf) == 1,na.rm=TRUE)
	q_posterior_beta = q_prior_beta + sum(c(hidden0,hiddenf)==0,na.rm=TRUE)
	if (q_posterior_alpha == 0) {
		q_posterior_alpha =1
	}
	qq <<- rbeta(1, q_posterior_alpha , q_posterior_beta)

	# propose q (beta distribution is conjugate distribution for binomial process)
	q_crossfamily_prior_alpha = 1;
	q_crossfamily_prior_beta = 1000;
	q_crossfamily_posterior_alpha = q_crossfamily_prior_alpha + sum(c(hidden_crossfamily0,hidden_crossfamilyf) == 1,na.rm=TRUE)
	q_crossfamily_posterior_beta = q_crossfamily_prior_beta + sum(c(hidden_crossfamily0,hidden_crossfamilyf)==0,na.rm=TRUE)
	if (q_crossfamily_posterior_alpha == 0) {
		q_crossfamily_posterior_alpha = q_crossfamily_posterior_beta * 10^-3
	}
	qq_crossfamily <<- rbeta(1, q_crossfamily_posterior_alpha , q_crossfamily_posterior_beta)

	
	#  update dvect (approximate using geometric distribution)
	# only if there is at least 1 recrudescent infection
	if (sum(classification==1) >= 1) {
	d_prior_alpha = 0;
	d_prior_beta = 0;
	mindistance_without_crossfamily = mindistance
	mindistance_without_crossfamily[mindistance_without_crossfamily > thresholddistance] = NA

	d_posterior_alpha = d_prior_alpha + sum(!is.na(c(mindistance_without_crossfamily[classification==1,])))
	d_posterior_beta = d_prior_beta + sum(c(round(mindistance_without_crossfamily[classification==1,])),na.rm=TRUE)
	if (d_posterior_beta == 0) {
		d_posterior_beta = sum(c((mindistance_without_crossfamily[classification==1,])),na.rm=TRUE)
	}
	if (d_posterior_beta == 0) { ## algorithm will get stuck if dposterior is allowed to go to 1
		d_posterior_beta = 1
	}	


	dposterior <<- rbeta(1, d_posterior_alpha , d_posterior_beta)

	### account for case that dposterior == 1
	if (dposterior > 1-1e-6) {
		dposterior = 1-1e-6
	}

	dvect = (1-dposterior) ^ (1:length(dvect)-1) * dposterior
	dvect <<- dvect / (sum(dvect))
	}


	# update frequencies
	# remove recrudescing alleles from calculations
	tempdata = recoded0
	sapply(which(classification == 1), function (x) tempdata[x,recr0[x,]] <<- 0)
	tempdata = rbind(tempdata, recodedf)
	sapply(1:nloci, function (x) findposteriorfrequencies(x,rbind(tempdata,recoded_additional_neutral)))


	### update correction distance matrix
	for (i in 1:nloci) {
		correction_distance_matrix[[i]] <<- correction_distance_matrix_constant[[i]] + abs(rnorm(dim(correction_distance_matrix_constant[[i]])[1],mean=0,frequencies_RR[[3]][i]))
	}

	# record state
	if (count > burnin & count %% record_interval == 0) {
		print(count)
		state_classification[,(count-burnin)/record_interval] <<- classification
		state_alleles0[,,(count-burnin)/record_interval] <<- alleles0
		state_allelesf[,,(count-burnin)/record_interval] <<- allelesf
		state_parameters[1,(count-burnin)/record_interval] <<- qq
		state_parameters[2,(count-burnin)/record_interval] <<- qq_crossfamily
		state_parameters[3,(count-burnin)/record_interval] <<- dposterior
		state_parameters[4:(4+nloci-1),(count-burnin)/record_interval] <<- apply(frequencies_RR[[2]],1,max)
		state_parameters[(4+nloci):(4+2*nloci-1),(count-burnin)/record_interval] <<- sapply(1:nloci,function (x) sum(frequencies_RR[[2]][x,]^2))
	}
	count <<- count + 1
}

replicate(nruns,runmcmc())

## make sure no NAs in result matrices
state_parameters = state_parameters[,!is.na(colSums(state_parameters))]
state_classification = state_classification[,!is.na(colSums(state_classification))]

## find mode of hidden alleles
modealleles = matrix("",2*nids,maxMOI*nloci)
for (i in 1:nids) {
	for (j in 1:nloci) {
		modealleles[2*(i-1)+1,((j-1)*maxMOI+1):(j*maxMOI)] = sapply(1:maxMOI, function (x) names(table(state_alleles0[i,(j-1)*maxMOI+x,]))[table(state_alleles0[i,(j-1)*maxMOI+x,])== max(table(state_alleles0[i,(j-1)*maxMOI+x,]))][1])
		modealleles[2*(i-1)+2,((j-1)*maxMOI+1):(j*maxMOI)] = sapply(1:maxMOI, function (x) names(table(state_allelesf[i,(j-1)*maxMOI+x,]))[table(state_allelesf[i,(j-1)*maxMOI+x,])== max(table(state_allelesf[i,(j-1)*maxMOI+x,]))][1])
	}
}

rowMeans2 = function(x){
	if (length(dim(x)) == 0) {
		ret = mean(x)
	} else {
		ret = rowMeans(x)
	}
	ret
}

temp_combined = c(sapply(1:length(ids), function (x) rep(rowMeans2(state_classification)[x],2)))
outputmatrix = cbind(temp_combined,modealleles)
colnames(outputmatrix) = c("Prob Rec",c(sapply(1:nloci, function (x) paste(locinames[x],"_",1:maxMOI,sep=""))))
hist(as.numeric(temp_combined))
write.csv(outputmatrix, paste(output_folder, "/", jobname,"_posterior",".csv",sep=""))


# summary statistics of parameters
write.csv(state_parameters, paste(output_folder, "/", jobname,"_state_parameters",".csv",sep=""))

summary_statisticsmatrix = cbind(format(rowMeans(state_parameters),digits=2),
					   apply(format(t(sapply(1:dim(state_parameters)[1], function (x) quantile(state_parameters[x,],c(0.25,0.75)))),digits=2),1, function (x) paste(x,collapse="?")))
summary_statisticsmatrix = rbind(summary_statisticsmatrix, c(format(mean(state_parameters[(3+nloci):(3+2*nloci-1),]),digits = 2),paste(format(quantile(state_parameters[(3+nloci):(3+2*nloci-1),],c(0.25,0.75)),digits=2),collapse="?")))
summary_statisticsmatrix = as.matrix(sapply(1:dim(summary_statisticsmatrix)[1], function (x) paste(summary_statisticsmatrix[x,1], " (",summary_statisticsmatrix[x,2],")",sep="")))
rownames(summary_statisticsmatrix) = c("q","q_crossfamily","d",locinames,locinames,"Mean diversity")
write.csv(summary_statisticsmatrix, paste(output_folder, "/", jobname,"_summarystatistics",".csv",sep=""))
