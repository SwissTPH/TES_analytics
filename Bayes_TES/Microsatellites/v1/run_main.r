

# necessary libraries
library(gtools)
library(readxl)
library(cluster)

#########################
# USER INPUT
set.seed(1)
setwd("/scicore/home/pothin/golmon00/GitRepos/STPHrepos/TES_analytics/Bayes_TES/Microsatellites/v1/")

# Name of input file
inputdata = "PPQ_63Days_simulated.xlsx"
# Number of hours to run, 1h minimum recommended
hours = 1
# Precision of fragment length measurement (2 or 3 for capillary electrophoresis for microsatelites/msp, 10-25 for gels)
locirepeats = c(3,3,3)

###########################


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

