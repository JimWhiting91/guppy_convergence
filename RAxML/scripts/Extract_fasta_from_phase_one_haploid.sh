#!/bin/bash

# Set up the environment
module load Java/1.8.0_74
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 BCFtools/1.9-foss-2018b
module load HTSlib/1.8-foss-2018a
module load SAMtools/1.9-foss-2018b

# Set up vars
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/raxml
DATASET=five_aside_STAR
POPMAP_DIR=$MASTER/data
PHASE_DIR=/gpfs/ts0/home/jw962/HP_LP_trials/phasing/phased_vcfs/$DATASET
CHR=000094F_0
START=500000
END=1100002
VCF=$PHASE_DIR/${DATASET}_${CHR}_shapeit_beagle.vcf.gz
GENOME=~/guppy_research/STAR/STAR.chromosomes.release.fasta

#POP_ARRAY=(five_aside_STAR)
#POP_ARRAY=(scaf94_haplogroup1_v2 scaf94_haplogroup3_v2)
POP_ARRAY=(haplogroup1_CL_AP_region haplogroup3_CL_AP_region)

# Make outputs dir
mkdir $MASTER/data/haplotypes/$DATASET $MASTER/data/haplotypes/$DATASET/${CHR}_${START}_${END}
OUT=$MASTER/data/haplotypes/$DATASET/${CHR}_${START}_${END}

# Tabix if not done already...
tabix -f -p vcf $VCF

# Phase per chr over population in loop
for pop in "${POP_ARRAY[@]}";
do

# Convert the phase to fasta
for ind in $(cat $POPMAP_DIR/${pop}.popmap)
do
rand=$(shuf -i 1-2 -n 1)
samtools faidx $GENOME $CHR:${START}-${END} | bcftools consensus $VCF -H ${rand}pIu -s $ind -o $OUT/${ind}_phase${rand}.fa
sed -i "s/>.*/>${ind}_phase${rand}/" $OUT/${ind}_phase${rand}.fa
done

# Cat them all together and remove bad symbols, taking a random allele from each individual
cat $OUT/*_phase*.fa > $OUT/${pop}_${DATASET}_${CHR}_${START}_${END}_haps.fa

# Tidy up
rm -f $OUT/*phase* $MASTER/outputs/${pop}_${DATASET}_${CHR}_${START}_${END}_tmp*
done

# Merge final
cat $OUT/*_${DATASET}_${CHR}_${START}_${END}_haps.fa > $OUT/ALL_HAPLOID_${DATASET}_${CHR}_${START}_${END}.fa
