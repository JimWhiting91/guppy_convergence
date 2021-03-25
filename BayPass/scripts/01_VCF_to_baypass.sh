#!/bin/bash
#PBS -d .
#PBS -q sq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -A Research_Project-T110748
#PBS -N baypass_input_maker
#PBS -e logs/baypass_input_maker.err.txt
#PBS -o logs/baypass_input_maker.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk


# Script makes input files for baypass based on VCFs, sorted by chromosome

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass
DATASET=five_aside_STAR
#VCF_DIR=/gpfs/ts0/home/jw962/HP_LP_trials/phasing/phased_vcfs/$DATASET
VCF=/gpfs/ts0/home/jw962/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_3033083_final.vcf.gz

# Define the population
pop_array=(GH GL OH OL MADHP MADLP TACHP TACLP APHP APLP)

# Loop over popmaps
for POP in "${pop_array[@]}"
do

# Outputs go here
mkdir $MASTER/data/$POP

# Get counts from VCF files, adding the derived flag means that sites with an ancestral allele will be correctly ordered
vcftools --gzvcf $VCF --counts2 --keep ~/HP_LP_trials/popmaps/$DATASET/${POP}.popmap --out $MASTER/outputs/${DATASET}_${POP}_tmp

# Now we need to edit the format of the outputted file to Baypass standard - 2 cols per SNP, 1 row per SNP, Allele counts
tail -n+2 $MASTER/outputs/${DATASET}_${POP}_tmp.frq.count | cut -f 5,6 > $MASTER/data/${DATASET}_${POP}.geno

rm -f $MASTER/outputs/${DATASET}_${POP}_tmp*
done
