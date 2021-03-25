#!/bin/bash

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass/

VCF=~/guppy_research/five_aside_vcfs/five_aside_v3_pop_SNP.gatk.bi.depth4.miss.maf.final.filtered.AA.recode.vcf.gz
POPMAP_DIR=~/HP_LP_trials/popmaps/five_aside_v3/
pop_array=(GH GL OH OL MADHP MADLP APHP APLP TACHP TACLP)

for POP in "${pop_array[@]}"
do

# Get counts from VCF files, adding the derived flag means that sites with an ancestral allele will be correctly ordered
vcftools --gzvcf $VCF --counts2 --keep $POPMAP_DIR/${POP}.popmap --out $MASTER/outputs/${POP}_tmp

# Edit the output file to get final counts
tail -n+2 $MASTER/outputs/${POP}_tmp.frq.count | cut -f5,6 > $MASTER/data/${POP}_five_aside_v3.geno

# Tidy
rm -f $MASTER/outputs/${POP}_tmp*

done
