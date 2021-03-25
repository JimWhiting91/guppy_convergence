#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -A Research_Project-T110748
#PBS -N baypass_run
#PBS -e logs/baypass_run.err.txt
#PBS -o logs/baypass_run.out.txt
#PBS -V
#PBS -t 1-16
#PBS -m e -M j.whiting2@exeter.ac.uk

module load ifort/2017.4.196-GCC-6.4.0-2.28

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass

# This script runs BayPass as an array job on the HPC, running each of the 16 subsets on their own CPU
# Assumes that the inputs have already been subsetted
POP=$1

# Run XtX
baypass -npop 2  -seed ${MOAB_JOBARRAYINDEX} -gfile $MASTER/data/subsets/five_aside_STAR_sub${MOAB_JOBARRAYINDEX}_${POP}.geno -omegafile $MASTER/data/five_aside_STAR_${POP}_CovMatrix.txt -outprefix $MASTER/outputs/${POP}_five_aside_STAR_sub${MOAB_JOBARRAYINDEX} -nthreads 8

# Run with Env covariate e.g. HP-LP BUT WITH AUX model that can account for spatial SNPs
#baypass -npop 10 -seed ${MOAB_JOBARRAYINDEX} -auxmodel -isingbeta 1.0 -gfile $MASTER/data/subsets/${POP}_sub${MOAB_JOBARRAYINDEX}.geno -omegafile $MASTER/data/${POP}_CovMatrix_Final -efile $MASTER/data/${POP}.env -outprefix $MASTER/outputs/${POP}_auxmodel_sub${MOAB_JOBARRAYINDEX} -nthreads 8
