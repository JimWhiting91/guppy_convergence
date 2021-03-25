#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=4
#PBS -A Research_Project-T110748
#PBS -N chromopainter_test
#PBS -e logs/chromopainter_run.err.txt
#PBS -o logs/chromopainter_run.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk
#PBS -t 1-33

# Set up environment
MASTER=~/HP_LP_trials/fineSTRUCTURE/
DATASET="OH_to_paint_APHP_OL"

CHR=$(awk '$2 > 1000000' ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | sed "${MOAB_JOBARRAYINDEX}q;d")

# I have added this so that the libgsl.so.0 symlink is in the library path. Links to ~./linuxbrew/libs/libgsl.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/software/fs_4.0.1

# Load modules
module load GLib/2.53.5-GCCcore-6.4.0 GSL/2.4-GCCcore-6.4.0 Perl/5.28.0-GCCcore-7.3.0

# Make output DIR
mkdir $MASTER/outputs/chromopainter/${DATASET}

# This is a test of chromosome painting the sex chromosome - The -b flag gives us the assignment prob per locus that we want...
chromopainter -g $MASTER/data/chromopainter/${DATASET}/${DATASET}_${CHR}_chromoIN.phase \
-r $MASTER/data/chromopainter/${DATASET}/${DATASET}_${CHR}_chromoIN.rec \
-f $MASTER/data/chromopainter/${DATASET}/${DATASET}.donors \
-i 10 -in -im -ip -b \
-o $MASTER/outputs/chromopainter/${DATASET}/${DATASET}_${CHR}_RESULTS
