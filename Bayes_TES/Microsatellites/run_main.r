
# necessary libraries
library(gtools)
library(readxl)
library(cluster)

# input variables

inputdata = "~/genotyping/microsatellite_analysis/BayesianAlgorithm/Microsatellites_patient_samples_results_PfPk2_msp1_msp2_bayesian.xlsx" # name of input file
hours = 1 # number of hours to run
locirepeats = c(3,3,3) # precision of fragment length measurement (2 or 3 for capillary electrophoresis for microsatelites/msp, 10-25 for gels)

setwd("/scicore/home/pothin/golmon00/GitRepos/bayesian_algorithm/")

source("Import_msp1.r")


start_time = Sys.time()
nruns_sample = 1000;
nruns = nruns_sample
record_interval = ceiling(nruns_sample / 200);
burnin = nruns_sample / 10;
source("run.r")
end_time = Sys.time()
totaltime = end_time - start_time

nruns = hours * 60 * 60* nruns_sample / as.numeric(totaltime, units = "secs")
record_interval = ceiling(nruns / 1000);
burnin = ceiling(nruns * 0.25);
source("run.r")

