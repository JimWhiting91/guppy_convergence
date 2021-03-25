#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N fineSTRUC
#PBS -e logs/chromopainter.err.txt
#PBS -o logs/chromopainter.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

# -----------------------
# NOTE - THIS SCRIPT IS RUN FROM THE COMMANDLINE, THE QSUBBING OCCURS THROUGH AN EDITED VERSION OF FINESTRUCTURE'S qsub_run.sh script
# -----------------------

# Set up environment
MASTER=~/HP_LP_trials/fineSTRUCTURE/
FINE_S=/gpfs/ts0/home/jw962/software/fs_4.0.1/
DATASET="five_aside_STAR"

# I have added this so that the libgsl.so.0 symlink is in the library path. Links to ~./linuxbrew/libs/libgsl.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/software/fs_4.0.1

# Load modules
module load GLib/2.53.5-GCCcore-6.4.0
module load GSL/2.4-GCCcore-6.4.0
module load Perl/5.28.0-GCCcore-7.3.0
module load R/3.5.1-foss-2018b

echo '**** START PAINT/fineSTRUCTURE ****'

# remove duplicates from ids file
#uniq $MASTER/outputs/${DATASET}_chromoIN.ids > $MASTER/outputs/${DATASET}_chromoIN.uniq.ids

# Change
cd $MASTER/outputs

# -------------------------------------
# STAGES 1-2 TAKE APPROX 12 HOURS ON PQ
# -------------------------------------

echo '**** Write Command Files 1 ****'

fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp \
-hpc 1 \
-idfile ${DATASET}_chromoIN.ids \
-phasefiles ${DATASET}_chromoIN.phase \
-recombfiles ${DATASET}_chromoIN.rec \
-s1emits 10 -s1minsnps 10000 -go

# We need to edit the commandfiles to replace the fs function with where we store the new fs command
sed -i 's/fs /fineSTRUCTURE /g' $MASTER/outputs/${DATASET}_fineSTRUCTURE/commandfiles/commandfile1.txt

# We also have to edit the path to all of the inputs because the commandfiles assume that they are in the wrong place...

echo '**** Run Command Files 1 on HPC ****'
# NOTES -
# I have edited the qsub script a fair bit.
# The script must be run from the commandfiles directory, this is so the qsub_run script includes all the correct paths for outputs
# However the qsub_run script now includes a cd ../../ which brings the actual analysis back to the outputs directory so that all of the inputs are available (these do not have full paths!)
cd ${DATASET}_fineSTRUCTURE/commandfiles/

######################
# RUN STAGE 1 ON COMMAND LINE - Completes within XX hour walltime on pq
######################
$FINE_S/qsub_run.sh -f commandfile1.txt -n 16 -m 16 -w 3 -P parallel # 16 ppn per node, 16 jobs per node, 3 hours walltime
######################
echo '**** Finished Command Files 1 on HPC ****'

cd $MASTER/outputs

# Stage 1 successful, reset to stage 2
# Need to run a number of processes that doesn't overwhelm the RAM limitations. Inds per process of 20 yields 6 jobs which keeps RAM in check
fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -indsperproc 20 -go

# Reset if needs be
#fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -reset 2 -go

# Edit the command files again
sed -i 's/fs /fineSTRUCTURE /g'  $MASTER/outputs/${DATASET}_fineSTRUCTURE/commandfiles/commandfile2.txt

cd ${DATASET}_fineSTRUCTURE/commandfiles/
######################
# RUN ON STAGE 2 COMMAND LINE - TOTAL OF 9 JOBS TO RUN (INDS PER PROC = 20) Ran in ~72 hours
######################
${MASTER}/scripts/qsub_run_highmem.sh -f commandfile2.txt -n 16 -m 9 -w 168 -P parallel # 16 ppn per node, 9 jobs per node, 24 hours walltime
######################

cd $MASTER/outputs

#Stage 2 successful, reset to stage 3 (runs very quickly)
fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -reset 3 -go
fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -s3iters 1000000 -allowdep 1 -hpc 1 -go

# Edit and move
sed -i 's/fs/fineSTRUCTURE/g'  $MASTER/outputs/${DATASET}_fineSTRUCTURE/commandfiles/commandfile3.txt
cd ${DATASET}_fineSTRUCTURE/commandfiles/

######################
# RUN ON STAGE 3 COMMAND LINE - TOTAL OF 9 JOBS TO RUN (INDS PER PROC = 20) runs fast >6 hrs
######################
${MASTER}/scripts/qsub_run_highmem.sh -f commandfile3.txt -n 16 -w 96
######################

# five_aside = Conversion looks good
# holi_13 = Conversion looks good
# five_aside_v2 = Conversion looks good

# Stage 3 successful, reset to stage 4 (runs very quickly)
fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -reset 4 -go

#Need to edit commandfile4 so that 'fs fs' becomes 'fineSTRUCTURE fs'
sed -i 's/fs/fineSTRUCTURE/' $MASTER/outputs/${DATASET}_fineSTRUCTURE/commandfiles/commandfile4.txt

cd ${DATASET}_fineSTRUCTURE/commandfiles/
######################
# RUN ON STAGE 4 COMMAND LINE - I sometimes just run this in an interactive job, it is quick. NOTE - run from outputs/ obviously...
######################
${MASTER}/scripts/qsub_run.sh -f commandfile4.txt -n 16 -m 1 -w 6 -p
######################
echo '**** END PAINT/fineSTRUCTURE ****'

######################
# PLOT OUTPUT USING R
######################
Rscript ${MASTER}/R/Plot_fineSTRUCTURE_tree.R ${DATASET}
