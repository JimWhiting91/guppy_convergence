#!/bin/bash
#PBS -d . 
#PBS -q pq  
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N bamCoverage_scafs
#PBS -e logs/bamCoverage_scafs.err.txt 
#PBS -o logs/bamCoverage_scafs.out.txt 
#PBS -V
##PBS -m e -M j.whiting2@exeter.ac.uk

# The array script runs jobs over chr

# Set the environment
POP=$1
BAM_DIR=/gpfs/ts0/home/jw962/GUPPYCon/STAR_bams/
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/deeptools
DATASET=five_aside_STAR

# Load what we need, including the deeptools anaconda library
module load Anaconda3/5.2.0
source activate deeptools-env

# Loop over scaffolds that are at least 10 kb in size
for SCF in $(cat $MASTER/data/STAR.chromosomes.release.fasta.fai | awk '$2 > 100000' | cut -f1 | grep -v 'chr')
do

# Make output dir
mkdir $MASTER/outputs/coverage
mkdir $MASTER/outputs/coverage/$DATASET
mkdir $MASTER/outputs/coverage/$DATASET/$POP
OUT=$MASTER/outputs/coverage/$DATASET/$POP

# Make the 'blacklisted' regions of 5kb at the beginning and end of chr
CHR_END=$(grep -w "$SCF" $MASTER/data/BUCK.phase-0.fasta.fai | cut -f2)
CHR_END2=$(bc <<< "scale = 5; $CHR_END - 5000")
echo -e "${SCF}\t0\t5000\tNA\n${SCF}\t$CHR_END2\t$CHR_END\tNA" > $OUT/${SCF}_black.bed

# Loop over bams
# Note - The Effective Genome Size has been reduced by 20,000 to account for Blacklist
for bam in $(ls ${BAM_DIR}/*.bam | grep $POP)
do
bamCoverage -b $bam -bl $OUT/${SCF}_black.bed -o ${bam}.${SCF}.coverage.bed -p 16 --outFileFormat=bedgraph --binSize=50 --smoothLength=75 --effectiveGenomeSize=528541853 --region=${SCF} --normalizeUsing=RPGC
mv ${bam}.${SCF}.coverage.bed $OUT
done

# Tidy up
rm -f $OUT/${SCF}_black.bed

# End scaffold loop
done

# Exit environment
source deactivate
