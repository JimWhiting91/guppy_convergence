#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=4
#PBS -A Research_Project-T110748
#PBS -N Coverage_ratio_est_run
#PBS -e logs/Coverage_ratio_est.err.txt
#PBS -o logs/Coverage_ratio_est.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

module load R/3.5.1-foss-2018b

Rscript /gpfs/ts0/home/jw962/HP_LP_trials/deeptools/R/01_compare_HP_LP_coverage.R
