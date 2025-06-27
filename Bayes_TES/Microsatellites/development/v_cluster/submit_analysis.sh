#!/bin/bash
#SBATCH --job-name=bayes_tes
#SBATCH --account=pothin
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/scicore/home/pothin/golmon00/commodities_forecasts/JOB_OUT/bayesian%a_%A.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/commodities_forecasts/JOB_OUT/bayesian%a_%A.o
#SBATCH --qos=6hours

########################
# Submit Bayesian algorithm
# 24.02.2024
# monica.golumbeanu@unibas.ch
#
# sbatch --array=1-16 submit_analysis.sh
#######################

module purge
module load R/4.1.2-foss-2018b-Python-3.6.6

ID=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
cd /scicore/home/pothin/golmon00/GitRepos/STPHrepos/TES_analytics/Bayes_TES/Microsatellites/v_cluster/

Rscript run_main.r $ID
