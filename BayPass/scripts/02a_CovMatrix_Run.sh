#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N baypass_core_run
#PBS -e logs/baypass_core_run.err.txt
#PBS -o logs/baypass_core_run.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

######################################
# This script filters the whole genome snp set for linkage and then runs the core model over this filtered dataset to output a covariance matrix
######################################

# We need these...
module load PLINK/2.00-alpha1-x86_64 ifort/2017.4.196-GCC-6.4.0-2.28 VCFtools/0.1.15-foss-2018a-Perl-5.26.1

# Environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/baypass
VCF=/gpfs/ts0/home/jw962/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_3033083_final.vcf.gz
DATASET=five_aside_STAR
POP_N=10

# Filter for linkage with plink
plink2 --vcf ${VCF} \
--out $MASTER/outputs/${DATASET}_plink_out_pruned --indep-pairwise 50 5 0.2 --allow-extra-chr

# Merge .geno with SNPs
gunzip -c $VCF | grep -v "#" | cut -f3 | paste - $MASTER/data/${DATASET}.geno > $MASTER/outputs/snp_labelled_geno.txt

# Filter .geno with pruned SNPs
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' $MASTER/outputs/${DATASET}_plink_out_pruned.prune.in $MASTER/outputs/snp_labelled_geno.txt | cut -f2- > $MASTER/data/${DATASET}_LD_pruned.geno

# Dry run of core model
baypass -npop $POP_N -gfile $MASTER/data/${DATASET}_LD_pruned.geno -outprefix $MASTER/outputs/${DATASET}_core -nthreads 16
