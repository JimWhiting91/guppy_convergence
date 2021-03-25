#!/bin/bash
#PBS -q sq
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1
#PBS -A Research_Project-T110748
#PBS -N coverage_to_bed
#PBS -e logs/coverage_to_bed.err.txt
#PBS -o logs/coverage_to_bed.out.txt
#PBS -V
##PBS -m e -M j.whiting2@exeter.ac.uk
#PBS -t 1-267

# Load bedops
module load BEDOPS/2.4.26

# Coverage directory
DIR=/gpfs/ts0/home/jw962/HP_LP_trials/deeptools/outputs/coverage/josie_poolseq/pop_avgs_filtered
mkdir $DIR
cd $DIR

# Read the input and output from command line
WINDOW=1000

# Set Pop_ARRAY
pop_array=(APHP APLP GH GL LT UT LMD UMD LO UQ)

# Loop over chrs
CHR=$(cut -f1 ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | sed "${MOAB_JOBARRAYINDEX}q;d")

# Make reference
awk '{print $1,0,$2}' ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | bedops --chop $WINDOW - | grep -w $CHR > winds_${CHR}.bed

# Loop over pops
for POP in "${pop_array[@]}"
do

# Set output
OUTPUT=${POP}_${CHR}_${WINDOW}_coverage

# Turn input into bed
awk -v CHR=$chr -v OFS="\t" 'NR>1 {print $1, $2, $2+50, 0, $3}' ${POP}_${CHR}_mashed2.bed > ${OUTPUT}_tmp.bed

# Do weighted mean within windows
bedmap --echo --wmean --skip-unmapped winds_${CHR}.bed ${OUTPUT}_tmp.bed | sed "s/|/\t/g" > ${OUTPUT}.bed

# Tidy
rm -f ${OUTPUT}_tmp.bed
done

# Tidy
rm -f winds_${CHR}.bed
