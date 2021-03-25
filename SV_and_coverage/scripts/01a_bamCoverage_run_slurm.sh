#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=bamCoverage_run
#SBATCH --error=logs/bamCoverage_run.err.txt
#SBATCH --output=logs/bamCoverage_run.out.txt
#SBATCH --export=All
#SBATCH --array=1-23
##SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

# The array script runs jobs over chr

# Set the environment
POP=$1
BAM_DIR=/gpfs/ts0/home/jw962/GUPPYCon/STAR_bams
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/deeptools
DATASET=five_aside_STAR

# Load what we need, including the deeptools anaconda library
module load Anaconda3/5.2.0
source activate deeptools-env

# Make output dir
mkdir $MASTER/outputs/coverage
mkdir $MASTER/outputs/coverage/$DATASET
mkdir $MASTER/outputs/coverage/$DATASET/$POP
OUT=$MASTER/outputs/coverage/$DATASET/$POP

# Make the 'blacklisted' regions of 100kb at the beginning and end of chr
CHR_END=$(grep -w "chr${SLURM_ARRAY_TASK_ID}" $MASTER/data/STAR.chromosomes.release.fasta.fai | cut -f2)
CHR_END2=$(bc <<< "scale = 5; $CHR_END - 100000")
echo -e "chr${SLURM_ARRAY_TASK_ID}\t0\t100000\tNA\nchr${SLURM_ARRAY_TASK_ID}\t$CHR_END2\t$CHR_END\tNA" > $OUT/chr${SLURM_ARRAY_TASK_ID}_black.bed

# Make a list of scaffolds to exclude for normalisation calculations
scafs=(`cat ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | grep -v "chr"`)
scaf_blacks=""i
for i in "${scafs[@]}"
do
    scaf_blacks+="$i "
done

# Loop over bams
# Note - The Effective Genome Size has been reduced by 100,000 to account for Blacklist
for bam in $(ls ${BAM_DIR}/*.bam | grep $POP)
do
bamCoverage -b $bam -bl $OUT/chr${SLURM_ARRAY_TASK_ID}_black.bed -o ${bam}.chr${SLURM_ARRAY_TASK_ID}.coverage.bed -p 4 --outFileFormat=bedgraph --binSize=50 --smoothLength=75 --effectiveGenomeSize=528351853 --region=chr${SLURM_ARRAY_TASK_ID} --normalizeUsing=RPGC
mv ${bam}.chr${SLURM_ARRAY_TASK_ID}.coverage.bed $OUT
done

# Tidy up
rm -f $OUT/chr${SLURM_ARRAY_TASK_ID}_black.bed

# Exit environment
source deactivate
