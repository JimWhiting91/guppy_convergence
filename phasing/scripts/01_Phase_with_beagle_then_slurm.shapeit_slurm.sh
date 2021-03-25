#!/bin/bash
#SBATCH -p pq
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=phase_VCF
#SBATCH --error=logs/phase_VCF.err.txt
#SBATCH --output=logs/phase_VCF.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk
#SBATCH --array=1-267%23

# Set up the environment
module load Java/1.8.0_74
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 BCFtools/1.9-foss-2018b
module load HTSlib/1.8-foss-2018a

# Which chr?
CHR=$(cat ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1)

# Set variables
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/phasing
POPMAP_DIR=/gpfs/ts0/home/jw962/HP_LP_trials/popmaps
DATASET=five_aside_STAR_v2_noMAF_lowNe
VCF=/gpfs/ts0/home/jw962/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_v2_noMAF_6526339_final.vcf.gz
BAM_DIR=~/GUPPYCon/STAR_bams
POP=five_aside_STAR

# Make the output
mkdir $MASTER/phased_vcfs/$DATASET
OUT=$MASTER/phased_vcfs/$DATASET

# Extract
bcftools view -r $CHR $VCF > $OUT/${POP}_${CHR}_tmp.recode.vcf

#-----------------------------
# Phase with beagle
java -Xss20m -jar ~/bin/beagle5.jar \
gt=$OUT/${POP}_${CHR}_tmp.recode.vcf \
out=$OUT/${POP}_${CHR}_beagle_phased.gt \
ne=10000 \
nthreads=12

# Tabix index
tabix $OUT/${POP}_${CHR}_beagle_phased.gt.vcf.gz

#-----------------------------
# Phase with shapeit
# First we need to list the BAMs
rm -f $MASTER/outputs/${POP}_${CHR}_bamlist
for ind in $(cat $POPMAP_DIR/${POP}.popmap)
do
ind2=$(echo -e ${ind}_)
bam_in=$(ls $BAM_DIR | grep "${ind2}" | grep ".bam" | grep -v ".bai" | grep -v "table")
echo -e "${ind}\t${bam_in}\t${CHR}" >> $MASTER/outputs/${POP}_${CHR}_bamlist
done

BAM_LIST=$MASTER/outputs/${POP}_${CHR}_bamlist

cd $BAM_DIR

# Extract phsae-informative reads
extractPIRs --bam $BAM_LIST \
            --vcf $OUT/${POP}_${CHR}_beagle_phased.gt.vcf.gz \
            --out $OUT/${POP}_${CHR}_beagle_PIRlist

# Now phase
shapeit -assemble \
	--force \
        --input-vcf $OUT/${POP}_${CHR}_beagle_phased.gt.vcf.gz \
        --input-pir $OUT/${POP}_${CHR}_beagle_PIRlist \
        -O $OUT/${POP}_${CHR}_shapeit_beagle_haps

# And convert to VCF
shapeit -convert \
        --input-haps $OUT/${POP}_${CHR}_shapeit_beagle_haps \
        --output-vcf $OUT/${POP}_${CHR}_shapeit_beagle.vcf

bgziptabix.sh $OUT/${POP}_${CHR}_shapeit_beagle.vcf

cd $MASTER

# Tidy
rm -f $OUT/${POP}_${CHR}_b* $BAM_LIST
