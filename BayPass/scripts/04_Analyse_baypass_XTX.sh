#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=00:20:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N baypass_analyse
#PBS -e logs/baypass_analyse.err.txt
#PBS -o logs/baypass_analyse.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

module load ifort/2017.4.196-GCC-6.4.0-2.28 R/3.5.1-foss-2018b

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass

# Eventually this will be an array I assume
POP=$1

# Run R script to assess significance and output outlier SNPs
Rscript $MASTER/R/get_xtx_significance.R $POP
