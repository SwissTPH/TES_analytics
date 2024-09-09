

# necessary libraries
library(gtools)
library(readxl)
library(cluster)

#########################
# USER INPUT
set.seed(1)
setwd("/scicore/home/pothin/golmon00/GitRepos/STPHrepos/TES_analytics/Bayes_TES/Microsatellites/v1/")

# Name of input file
inputdata = "~/genotyping/tests_Bayesian/Test2.xlsx"#"PPQ_63Days_simulated.xlsx"

# inputdata1 = "~/GitRepos/bayesian_analysis_Uganda/datasets/MSP GLURP_DP.xlsx"
# inputdata2 = "~/GitRepos/bayesian_analysis_Uganda/datasets/MSP 313_DP.xlsx"
# genotypedata_latefailures1 = read_excel(inputdata1, sheet=1)
# genotypedata_latefailures2 = read_excel(inputdata2, sheet=1)

# Number of hours to run, 1h minimum recommended
hours = 0.05
# Precision of fragment length measurement (2 or 3 for capillary electrophoresis for microsatelites/msp, 10-25 for gels)
locirepeats = c(10,10,2)

###########################


source("Import_msp1.r")


start_time = Sys.time()
nruns_sample = 100;
nruns = nruns_sample
record_interval = ceiling(nruns_sample / 200);
burnin = nruns_sample / 10;
source("run.r")
end_time = Sys.time()
totaltime = end_time - start_time

# nruns = hours * 60 * 60* nruns_sample / as.numeric(totaltime, units = "secs")
nruns=100
record_interval = ceiling(nruns / 1000);
burnin = ceiling(nruns * 0.25);
source("run.r")

