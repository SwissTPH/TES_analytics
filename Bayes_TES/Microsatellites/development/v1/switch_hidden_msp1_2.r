switch_hidden = function(x) {
	z = runif(1)
	if (sum(hidden0[x,], hiddenf[x,],na.rm=TRUE) > 0) { # if hidden alleles exist
		if (length(which(c(hidden0[x,], hiddenf[x,])==1))>1) {
			chosen = sample(which(c(hidden0[x,], hiddenf[x,])==1),1)
		} else {
			chosen = which(c(hidden0[x,], hiddenf[x,])==1)
		}
		if (classification[x] == 0) { # reinfection
			if (chosen <= nloci*maxMOI) { # day 0 hidden allele
				chosenlocus = ceiling(chosen/maxMOI)
				old = recoded0[x,chosen]
				new = sample(1:frequencies_RR[[1]][chosenlocus],1)


				oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
				repeatedold = qq
				repeatednew = qq
				if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
					repeatedold = 1;
				}
				if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
					repeatednew = 1;
				}
				old_alleles_lengths = alleles0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]

				crossfamily_old = hidden_crossfamily0[x,chosen]
				crossfamily_new = 0
				if (sum(!is.na(old_alleles_lengths))>0){
				if (min(abs((mode_allele_lengths[[chosenlocus]][new]) - old_alleles_lengths))>thresholddistance) {
					crossfamily_new = 1
				}}
				
				#alpha = (frequencies_RR[[2]][chosenlocus,new] * repeatednew) / (frequencies_RR[[2]][chosenlocus,old] * repeatedold)
				alpha = (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1]) * repeatednew)*qq_crossfamily^crossfamily_new / 
					   (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1]) * repeatedold)*qq_crossfamily^crossfamily_old
				if (z < alpha) { # switch made
					recoded0[x,chosen] <<- new
					####### new allele should have some variability (ie what's being see in real data)
					newallele_length = (mode_allele_lengths[[chosenlocus]][new]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
					#possible_allele_lengths = c(alleles0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recoded0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hidden0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0] ,
					#					   allelesf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recodedf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hiddenf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0])
					#if (length(possible_allele_lengths) == 1) {
					#	newallele_length = possible_allele_lengths 
					#} else {
					#	newallele_length = sample(possible_allele_lengths,1)
					#}
					alleles0[x,chosen] <<- newallele_length
					hidden_crossfamily0[x,chosen] <<- crossfamily_new
					# if (recr0[x,chosenlocus] == chosen) { # have chosen an allele that is marked as possibly recrudescing, so need to update mindistance
											    # actually need to do this for all substitutions
						allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
						closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
						mindistance[x,chosenlocus] <<- abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
						alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
						allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
						recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
						recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
						recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
						recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
					# }
				}
			} else { # day f hidden allele
				chosen = chosen - nloci*maxMOI
				chosenlocus = ceiling(chosen/maxMOI)
				old = recodedf[x,chosen]
				new = sample(1:frequencies_RR[[1]][chosenlocus],1)
				oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]
				repeatedold = qq
				repeatednew = qq
				if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
					repeatedold = 1;
				}
				if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
					repeatednew = 1;
				}

				old_alleles_lengths = allelesf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]

				crossfamily_old = hidden_crossfamilyf[x,chosen]
				crossfamily_new = 0
				if (sum(!is.na(old_alleles_lengths))>0){
				if (min(abs((mode_allele_lengths[[chosenlocus]][new]) - old_alleles_lengths))>thresholddistance) {
					crossfamily_new = 1
				}}

				#alpha = (frequencies_RR[[2]][chosenlocus,new] * repeatednew) / (frequencies_RR[[2]][chosenlocus,old] * repeatedold)
				alpha = (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1]) * repeatednew *qq_crossfamily^crossfamily_new) / 
					   (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1]) * repeatedold *qq_crossfamily^crossfamily_old)
				if (z < alpha) { # switch made
					recodedf[x,chosen] <<- new
					#newallele_length = (mode_allele_lengths[[chosenlocus]][new])
					newallele_length = (mode_allele_lengths[[chosenlocus]][new]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
					hidden_crossfamilyf[x,chosen] <<- crossfamily_new

					####### new allele should have some variability (ie what's being see in real data)
					
					#possible_allele_lengths = c(alleles0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recoded0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hidden0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0] ,
					#					   allelesf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recodedf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hiddenf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0])
					#if (length(possible_allele_lengths) == 1) {
					#	newallele_length = possible_allele_lengths 
					#} else {
					#	newallele_length = sample(possible_allele_lengths,1)
					#}

					allelesf[x,chosen] <<- newallele_length
					# if (recrf[x,chosenlocus] == chosen) { # update mindistance
						allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
						closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
						mindistance[x,chosenlocus] <<- abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
						alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
						allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
						recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
						recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
						recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
						recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
					# }
				}
			}
		} else { # recrudescence
			if (chosen <= nloci*maxMOI) { # day 0 hidden allele
				chosenlocus = ceiling(chosen/maxMOI)
				old = recoded0[x,chosen]
				new = sample(1:frequencies_RR[[1]][chosenlocus],1)
#				newallele_length = (mode_allele_lengths[[chosenlocus]][new])
				newallele_length = (mode_allele_lengths[[chosenlocus]][new]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])

				#possible_allele_lengths = c(alleles0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recoded0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hidden0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0] ,
				#						   allelesf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recodedf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hiddenf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0])
				#if (length(possible_allele_lengths) == 1) {
				#	newallele_length = possible_allele_lengths 
				#} else {
				#	newallele_length = sample(possible_allele_lengths,1)
				#}
				oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
				repeatedold = qq
				repeatednew = qq
				if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
					repeatedold = 1;
				}
				if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
					repeatednew = 1;
				}

				old_alleles_lengths = alleles0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]

				crossfamily_old = hidden_crossfamily0[x,chosen]
				crossfamily_new = 0
				if (sum(!is.na(old_alleles_lengths))>0){
				if (min(abs(newallele_length - old_alleles_lengths))>thresholddistance) {
					crossfamily_new = 1
				}}


				allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
				tempalleles = alleles0[x,maxMOI*(chosenlocus-1)+1:maxMOI]
				tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length 
				temprecoded = recoded0[x,maxMOI*(chosenlocus-1)+1:maxMOI]
				temprecoded[chosen-(chosenlocus-1)*maxMOI] = new

				newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
				newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]])
				newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
				newallrecrf = recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]

				# calculate new multiple-comparisons coefficient
				newrecr0 = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
				newrecrf = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
				newrecr_repeats0 = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud,1]],na.rm=TRUE)
				newrecr_repeatsf = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,newrecrf])

				#likelihoodnew = mean(dvect[round(newalldistance)+1]/frequencies_RR[[2]][chosenlocus,newallrecrf],na.rm=TRUE) * repeatednew
				likelihoodnew = mean(dvect[round(newalldistance)+1]/sapply(1:length(newallrecrf), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,newallrecrf[z]]+1])),na.rm=TRUE) * repeatednew *qq_crossfamily^crossfamily_new
				#likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/frequencies_RR[[2]][chosenlocus,allrecrf[x,chosenlocus,]],na.rm=TRUE) * repeatedold 
				likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,allrecrf[x,chosenlocus,z]]+1])),na.rm=TRUE) * repeatedold  *qq_crossfamily^crossfamily_old

				if (is.na(likelihoodnew) | is.na(likelihoodold)) { # debug
					write.csv(alleles0,"alleles0.csv")
					write.csv(allelesf,"allelesf.csv")
					write.csv(hidden0,"hidden0.csv")
					write.csv(hiddenf,"hiddenf.csv")
					write.csv(dvect,"dvect.csv")
					print(chosen)
					print(chosenlocus)
					print(chosenlocus)
					print(old)
					print(new)
					print(newallele_length)
					print(tempalleles)
					write.csv(newalldistance,"newalldistance.csv")
					write.csv(alldistance,"alldistance.csv")
				}

				if (likelihoodnew  == likelihoodold) {
					# if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
					alpha = 1							 
				} else {
					alpha = likelihoodnew / likelihoodold 						 
				}

				if (z < alpha) { # switch made
					recoded0[x,chosen] <<- new
					alleles0[x,chosen] <<- newallele_length
					hidden_crossfamily0[x,chosen] <<- crossfamily_new
					mindistance[x,chosenlocus] <<- newmindistance
					alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newalldistance
					allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newallrecrf
					recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
					recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
					recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
					recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
				}
			} else { # day f hidden allele
				chosen = chosen - nloci*maxMOI
				chosenlocus = ceiling(chosen/maxMOI)
				old = recodedf[x,chosen]
				new = sample(1:frequencies_RR[[1]][chosenlocus],1)
				#newallele_length = (mode_allele_lengths[[chosenlocus]][new])
				newallele_length = (mode_allele_lengths[[chosenlocus]][new]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
				#possible_allele_lengths = c(alleles0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recoded0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hidden0[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0] ,
				#						   allelesf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] [recodedf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == new & hiddenf[,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == 0])
				#if (length(possible_allele_lengths) == 1) {
				#	newallele_length = possible_allele_lengths 
				#} else {
				#	newallele_length = sample(possible_allele_lengths,1)
				#}
				oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]
				repeatedold = qq
				repeatednew = qq
				if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
					repeatedold = 1;
				}
				if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
					repeatednew = 1;
				}

				old_alleles_lengths = allelesf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]

				crossfamily_old = hidden_crossfamilyf[x,chosen]
				crossfamily_new = 0
				if (sum(!is.na(old_alleles_lengths))>0){
				if (min(abs(newallele_length - old_alleles_lengths))>thresholddistance) {
					crossfamily_new = 1
				}}

				allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
				tempalleles = allelesf[x,maxMOI*(chosenlocus-1)+1:maxMOI]
				tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length 
				temprecoded = recodedf[x,maxMOI*(chosenlocus-1)+1:maxMOI]
				temprecoded[chosen-(chosenlocus-1)*maxMOI] = new
				newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]])))
				newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]])
				newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]]))
				newallrecrf = temprecoded[allpossiblerecrud[,2]]

				# calculate new multiple-comparisons coefficient
				newrecr0 = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
				newrecrf = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
				newrecr_repeats0 = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,newrecr0])
				newrecr_repeatsf = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud,2]],na.rm=TRUE)

				#likelihoodnew = mean(dvect[round(newalldistance)+1]/frequencies_RR[[2]][chosenlocus,newallrecrf],na.rm=TRUE) * repeatednew
				likelihoodnew = mean(dvect[round(newalldistance)+1]/sapply(1:length(newallrecrf), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,newallrecrf[z]]+1])),na.rm=TRUE) * repeatednew *qq_crossfamily^crossfamily_new
				#likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/frequencies_RR[[2]][chosenlocus,allrecrf[x,chosenlocus,]],na.rm=TRUE) * repeatedold 
				likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,allrecrf[x,chosenlocus,z]]+1])),na.rm=TRUE) * repeatedold *qq_crossfamily^crossfamily_old

				if (is.na(likelihoodnew) | is.na(likelihoodold)) { # debug
					write.csv(alleles0,"alleles0.csv")
					write.csv(allelesf,"allelesf.csv")
					write.csv(hidden0,"hidden0.csv")
					write.csv(hiddenf,"hiddenf.csv")
					write.csv(dvect,"dvect.csv")
					print(chosen)
					print(chosenlocus)
					print(chosenlocus)
					print(old)
					print(new)
					print(newallele_length)
					print(tempalleles)
					write.csv(newalldistance,"newalldistance.csv")
					write.csv(alldistance,"alldistance.csv")
				}


				if (likelihoodnew  == likelihoodold) {
					# if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
					alpha = 1							 
				} else {
					alpha = likelihoodnew / likelihoodold 						 
				}
				if (z < alpha) { # switch made
					recodedf[x,chosen] <<- new
					allelesf[x,chosen] <<- newallele_length
					hidden_crossfamilyf[x,chosen] <<- crossfamily_new
					mindistance[x,chosenlocus] <<- newmindistance
					alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newalldistance
					allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newallrecrf
					recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
					recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
					recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
					recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
				}
			}
		}
	}
}
