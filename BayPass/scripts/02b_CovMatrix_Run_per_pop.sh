#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N baypass_core_run
#PBS -e logs/baypass_core_run.err.txt
#PBS -o logs/baypass_core_run.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk
#PBS -t 0-4

######################################
# This script filters the whole genome snp set for linkage and then runs the core model over this filtered dataset to output a covariance matrix
######################################

# We need these...
module load ifort/2017.4.196-GCC-6.4.0-2.28

# Environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass
POP_N=2
DATASET=five_aside_STAR

# Pops
pops=(G O AP MAD TAC)

POP=${pops[$MOAB_JOBARRAYINDEX]}

# Dry run of core model
baypass -npop $POP_N -gfile $MASTER/data/${DATASET}_LD_pruned_${POP}.geno -outprefix $MASTER/outputs/${DATASET}_${POP}_core -nthreads 16
