#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=plink_PCA_noCaroniPops
#SBATCH --error=logs/plink_PCA_noCaroniPops.err.txt
#SBATCH --output=logs/plink_PCA_noCaroniPops.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END 
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

module load PLINK/2.00-alpha1-x86_64
module load BCFtools/1.6-intel-2017b
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

echo "Start Plink"

MASTER=~/HP_LP_trials/fineSTRUCTURE/
VCF=~/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_3033083_final.vcf.gz
POPMAP_DIR=~/HP_LP_trials/popmaps/five_aside_STAR

cd $MASTER

# Make a vcf with IDs to use in plink
#bcftools annotate --set-id +'%CHROM\_%POS' $VCF > ${VCF}_IDs.vcf
CARONI_POPS=(Guanapo Aripo Tacarigua)
for POP in "${CARONI_POPS[@]}"
do

OUT=five_aside_STAR_No${POP}

# Filter for individuals
vcftools --gzvcf $VCF --remove $POPMAP_DIR/$POP.popmap --recode --out outputs/$OUT

# Make new
POP_VCF=outputs/${OUT}.recode.vcf

#Prune SNPs for LD
plink2 --vcf ${POP_VCF} \
--out $MASTER/outputs/${OUT}_plink_out_pruned --indep-pairwise 50 5 0.2 --allow-extra-chr

#Make Pruned Data
plink2 --allow-extra-chr --vcf ${POP_VCF} \
--extract $MASTER/outputs/${OUT}_plink_out_pruned.prune.in \
--make-bed --out $MASTER/outputs/${OUT}_plink_out_NoLD

# Prune with BCFtools
# Discard records with r2 bigger than 0.4, first removing records with more than 5% of genotypes missing
#bcftools +prune -l 0.4 -e'F_MISSING>=0.05' $VCF -Ob -o $VCF.pruned.vcf

#Run PCA
plink2 --bfile $MASTER/outputs/${OUT}_plink_out_NoLD --pca --out $MASTER/outputs/${OUT}_plink_out_NoLD_PCA_No${POP} --allow-extra-chr

#Convert plink pruned to vcf NOTE: ASSUMES REFERENCE ALLELE SO NOT COMPARABLE TO ORIGINAL VCF
#plink --ped five_aside_plink.ped \
#--map five_aside_plink.map \
#--extract five_aside_plink_out_pruned.prune.in \
#--out five_aside_plinkPRUNED_vcf --recode vcf --noweb --allow-extra-chr

# Remove the pop specific VCF
rm -f $POP_VCF

done

echo "End Plink"
