#!/bin/bash
#PBS -q pq
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N run_smoove
#PBS -e logs/run_smoove.err.txt
#PBS -o logs/run_smoove.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk
#PBS -l excludenodes=comp164
##PBS -t 3-4

####################################################################################################
# SCRIPT TO RUN SMOOVE (ILLUMINA) FOR SV DETECTION FROM NATURAL DATA
####################################################################################################

# We need these...
module purge
module load Anaconda3/5.2.0 BCFtools/1.9-foss-2018b GCCcore/9.3.0

# Load the env
source activate smoove-env

# Set paths and general inputs
BAMS=/gpfs/ts0/home/jw962/GUPPYCon/PCR_free_bams/
REF=/gpfs/ts0/home/jw962/guppy_research/STAR/STAR.chromosomes.release.fasta
BAD=/gpfs/ts0/home/jw962/guppy_research/STAR/STAR.chromosomes.release.repeats.bed

# Set up working environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/SV/
cd $MASTER

# Choose pop from array
pops=(Tacarigua Guanapo Aripo Oropouche Madamas)
POP=${pops[$MOAB_JOBARRAYINDEX]}
POPMAP=/gpfs/ts0/home/jw962/HP_LP_trials/popmaps/five_aside_STAR/${POP}.popmap

# For small cohort, we run smoove over single command
rm -f outputs/${POP}_bams.txt
for ind in $(cat $POPMAP)
do
# ls $BAMS | grep "${ind}_" | grep -v ".bam." >> outputs/${POP}_bams.txt
ls $BAMS | grep "${ind}." | grep -v ".bai" >> outputs/${POP}_bams.txt
done

# Set up the command flag
genotype_flag=""
for bam in $(cat outputs/${POP}_bams.txt)
do
  genotype_flag+="--genotype ${BAMS}/${bam} "
done

# Run smoove
smoove call -x --name ${POP}_SV --fasta $REF --exclude $BAD --processes 16 --outdir $MASTER/outputs/$POP $genotype_flag

source deactivate smoove-env
