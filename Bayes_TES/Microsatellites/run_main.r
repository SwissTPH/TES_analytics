# Main function for running the Bayesian algorithm for length-polymorphic markers

# Necessary libraries
library(gtools)
library(readxl)
library(cluster)

# Additional functions
source("define_alleles.r")
source("calculate_frequencies3.r")
source("recode_alleles.r")
source("switch_hidden_msp1_2.r")
source("findposteriorfrequencies.r")

# Options
# memory.limit(size=50000)
options(java.parameters = "-Xmx4096m")

#########################
# USER INPUT
# To be modified according to user's file system
setwd("/scicore/home/pothin/golmon00/GitRepos/STPHrepos/TES_analytics/Bayes_TES/Microsatellites/")
set.seed(1)
# Name of input file and output folder
# If using a different file than the example, to be modified according to user's file system
inputdata = "~/genotyping/analysis_Uganda/datasets/MSP GLURP_AL.xlsx" #"PPQ_63Days_simulated.xlsx"
# To be modified according to user's file system
output_folder = "~/genotyping/analysis_Uganda/results/"
# Precision of fragment length measurement for each marker
# (2 or 3 for capillary electrophoresis for microsatelites/msp, 10-25 for gels)
locirepeats = c(3, 3, 3)

###########################

##### Read in the data
genotypedata_latefailures = read_excel(inputdata, sheet=1)

genotypedata_latefailures[genotypedata_latefailures == 0] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "0"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "N/A"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "-"] = NA # missing data has to be coded as NA
genotypedata_latefailures[genotypedata_latefailures == "NA"] = NA # missing data has to be coded as NA

genotypedata_latefailures$Day[is.na(genotypedata_latefailures$Day)] = 0

sampleid = paste(genotypedata_latefailures$PatientID,"D",genotypedata_latefailures$Day,sep="")
genotypedata_latefailures = cbind(Sample.ID = sampleid, genotypedata_latefailures[,-c(1,2)])
### recode sample names so that each pair has a " Day 0" and a " Day Failure"
genotypedata_latefailures$Sample.ID = sub("D0$"," Day 0", genotypedata_latefailures$Sample.ID)
genotypedata_latefailures$Sample.ID = sub("D[0-9]+$"," Day Failure", genotypedata_latefailures$Sample.ID)

# each sample in genotypedata_RR has to have day 0 and day of Failure
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day 0",genotypedata_latefailures$Sample.ID)]," Day 0")))
if (sum(!paste(ids, "Day Failure") %in% genotypedata_latefailures$Sample.ID) > 0) {
  print("Error - each sample must have day 0 and day of failure data")
}
ids = unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day Failure",genotypedata_latefailures$Sample.ID)]," Day Failure")))
if (sum(!paste(ids, "Day 0") %in% genotypedata_latefailures$Sample.ID) > 0) {
  print("Error - each sample must have day 0 and day of failure data")
}

additional_genotypedata = read_excel(inputdata,sheet=2)
sampleid = paste(additional_genotypedata$PatientID,"D",additional_genotypedata$Day,sep="")
additional_genotypedata = cbind(Sample.ID = sampleid,additional_genotypedata[,-c(1,2)])


nruns = 100000
record_interval = ceiling(nruns/1000);
burnin = ceiling(nruns * 0.25);
source("run.r")
