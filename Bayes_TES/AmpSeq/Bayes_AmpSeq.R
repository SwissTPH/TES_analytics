###########################
# Implementation of the Bayesian algorithm for AmpSeq data
#
# 28.02.2025
###########################

library(gtools)
library(readxl)

# input file
setwd("~/GitRepos/STPHrepos/TES_analytics/Bayes_TES/AmpSeq/")
inputfile = "Haplotype_data.xlsx"

# bin sizes (one for each marker)
# this is the tolerance for each marker, used during match counting
bin_sizes = c(0,0,0)

# this is the size of the repeat region for each microsatellite (ie. for dimers use 2, for trimers use 3)
locirepeats = c(3,3,3)

# number of runs for MCMC algorithm
nruns = 100#5000

##### read in data


hap_data = read_xlsx(inputfile)

vec_hap_original = unique(unlist(c(hap_data[,3:ncol(hap_data)])))
vec_hap_original = vec_hap_original[grepl("-",vec_hap_original)]
vec_hap_numeric = seq(1:length(vec_hap_original))*3+100

hap_data_numeric = apply(hap_data[, 3:ncol(hap_data)],2, function (x) vec_hap_numeric[match(x,vec_hap_original)])

genotypedata_latefailures = cbind(hap_data[,c(1,2)],hap_data_numeric)

genotypedata_latefailures[genotypedata_latefailures == 0] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "0"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "N/A"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "-"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "NA"] = NA # missing data has to be coded as NA

### recode sample names so that each pair has a " Day 0" and a " Day Failure"
genotypedata_latefailures$Sample.ID = sub("D0$"," Day 0",genotypedata_latefailures$Sample.ID)
genotypedata_latefailures$Sample.ID = sub("D[0-9]+$"," Day Failure",genotypedata_latefailures$Sample.ID)

genotypedata_latefailures = genotypedata_latefailures[,grepl("Sample", colnames(genotypedata_latefailures)) | grepl("_", colnames(genotypedata_latefailures))]

# each sample in genotypedata_latefailures has to have day 0 and day of Failure
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day 0",genotypedata_latefailures$Sample.ID)]," Day 0")))
if (sum(!paste(ids, "Day Failure") %in% genotypedata_latefailures$Sample.ID) > 0) {
    print("Error - each sample must have day 0 and day of failure data")
}
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day Failure",genotypedata_latefailures$Sample.ID)]," Day Failure")))
if (sum(!paste(ids, "Day 0") %in% genotypedata_latefailures$Sample.ID) > 0) {
    print("Error - each sample must have day 0 and day of failure data")
}

ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day 0",genotypedata_latefailures$Sample.ID)]," Day 0")))
locinames = unique(sapply(colnames(genotypedata_latefailures)[-1],function(x) strsplit(x,"_")[[1]][1]))
nloci = length(locinames)
nids = length(ids)


##### calculate MOI for each sample
MOI0 = rep(0,nids)
MOIf = rep(0,nids)
for (i in 1:nids) {
    for (j in 1:nloci) {
        locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata_latefailures))
        nalleles0 = sum(!is.na(genotypedata_latefailures[grepl(paste(ids[i],"Day 0"),genotypedata_latefailures$Sample.ID),locicolumns]))
        nallelesf = sum(!is.na(genotypedata_latefailures[grepl(paste(ids[i],"Day Failure"),genotypedata_latefailures$Sample.ID),locicolumns]))

        MOI0[i] = max(MOI0[i],nalleles0)
        MOIf[i] = max(MOIf[i],nallelesf)
    }
}
maxMOI = max(c(MOI0, MOIf),na.rm=TRUE)


alleles0 = matrix(NA,nids,maxMOI*nloci)
allelesf = matrix(NA,nids,maxMOI*nloci)
mindistance = matrix(NA,nids,nloci)
alldistance = array(NA,c(nids,nloci,maxMOI*maxMOI))

## read in allele data into usable R objects (arrays)
for (j in 1:nloci) {
    locus = locinames[j]
    locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata_latefailures))
    oldalleles = as.vector(genotypedata_latefailures[,locicolumns])
    if (length(dim(oldalleles)[2]) == 0) {
        oldalleles = matrix(oldalleles,length(oldalleles),1)
    }
    oldalleles = matrix(as.numeric(unlist(c(oldalleles))),dim(oldalleles)[1],dim(oldalleles)[2])
    oldalleles[is.na(oldalleles)] = NA

    alleles0[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day 0",genotypedata_latefailures$Sample.ID),]
    allelesf[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day Failure",genotypedata_latefailures$Sample.ID),]

}

number_matches = rep(NA,nids)
match_output = matrix("",nids*2,nloci)
colnames(match_output) = locinames
number_loci = rep(NA, nids)
## count matches
for (i in 1:nids) {
    nmatches_temp = 0
    nloci_temp = 0
    for (j in 1:nloci) { # determine which alleles are recrudescing (for beginning, choose closest pair)
        allpossiblerecrud = expand.grid(1:MOI0[i],1:MOIf[i])
        if (sum(!is.na(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + MOI0[i])])) > 0 & sum(!is.na(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + MOIf[i])])) > 0){
            nloci_temp = nloci_temp+1

            closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (x) abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]])))
            mindistance[i,j] = abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]])
            if (mindistance[i,j] <= bin_sizes[j])
            {
                nmatches_temp=nmatches_temp+1
                match_output[2*(i-1)+1,j] = "R"
            } else {
                match_output[2*(i-1)+1,j] = "NI"
            }
        } else {
            match_output[2*(i-1)+1,j] = "IND"
        }
    }
    number_matches[i]=nmatches_temp
    number_loci[i] = nloci_temp
}
