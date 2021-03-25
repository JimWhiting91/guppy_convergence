#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=prep_five_aside_beast_poec
#SBATCH --error=logs/prep_five_aside_beast_poec.err.txt
#SBATCH --output=logs/prep_five_aside_beast_poec.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

# Set environment
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 GATK/3.8-0-Java-1.8.0_144 BCFtools/1.9-foss-2018b
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/SNAPP
cd $MASTER

# First Prepare the input VCF
VCF=~/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_v3_picta_wingei_maf01_noINV_4829912.no_picta.merged.vcf.gz
REFERENCE_GENOME=~/guppy_research/STAR/STAR.chromosomes.release.fasta

RUN_NAME=five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2

# Retain top 3 individuals with lowest missing data
vcftools --gzvcf $VCF --missing-indv --out outputs/$RUN_NAME.missing

# Sort output and fetch 3 best
sort -nk5 outputs/$RUN_NAME.missing.imiss > outputs/$RUN_NAME.sorted.imiss
pops=(LT UT GH GL APHP APLP LO UQ LMD UMD)
rm -f outputs/$RUN_NAME.missing_keep.txt
cat ~/HP_LP_trials/popmaps/five_aside_STAR/WINGEI.popmap > outputs/$RUN_NAME.missing_keep.txt
for pop in "${pops[@]}"
do
  tail -n+2 outputs/$RUN_NAME.sorted.imiss | grep $pop | head -3 | cut -f1 >> outputs/$RUN_NAME.missing_keep.txt
done

# Subset the VCF + invariant filter again
bcftools view -S outputs/$RUN_NAME.missing_keep.txt -c 1 -c 1:nonmajor $VCF > outputs/$RUN_NAME.imiss_filtered.vcf

# Remove missing data
vcftools --vcf outputs/$RUN_NAME.imiss_filtered.vcf --max-missing 1 --recode --recode-INFO-all --out outputs/pic_wing_${RUN_NAME}.nomiss

# Thin, don't linkage...
vcftools --vcf outputs/pic_wing_${RUN_NAME}.nomiss.recode.vcf --thin 50000 --recode --out outputs/pic_wing_${RUN_NAME}.nomiss.thin50k

# Now just take a random sample of the thinned VCF for 1000 SNPs
grep "#" outputs/pic_wing_${RUN_NAME}.nomiss.thin50k.recode.vcf > header
grep -v "#" outputs/pic_wing_${RUN_NAME}.nomiss.thin50k.recode.vcf | grep -v "chr12" | shuf | head -n1000 | cat header - | bcftools sort > data/pic_wing_${RUN_NAME}.nomiss.thin50k.rand1000.recode.vcf
rm -f header

# Linkage filter
#bcftools +prune -w 50 -l 0.2 -e'F_MISSING>=0.02' outputs/pic_wing_${RUN_NAME}.nomiss.recode.vcf -Ov > outputs/pic_wing_${RUN_NAME}.nomiss.noLD.vcf

# Get per site Fst from R
#Rscript $MASTER/R/find_top_snps_for_vcf.R

# Final filter to get high Fst as well
bcftools annotate --set-id +'%CHROM\_%POS' outputs/pic_wing_${RUN_NAME}.nomiss.thin50k.recode.vcf | bcftools view -i"ID=@data/pic_wing_${RUN_NAME}_nochr12_highFST_RIVER.txt" - > data/pic_wing_${RUN_NAME}.nomiss.thin50k.highFST.nochr12.RIVER.vcf
bcftools annotate --set-id +'%CHROM\_%POS' outputs/pic_wing_${RUN_NAME}.nomiss.thin50k.recode.vcf | bcftools view -i"ID=@data/pic_wing_${RUN_NAME}_nochr12_highFST_HPLP.txt" - > data/pic_wing_${RUN_NAME}.nomiss.thin50k.highFST.nochr12.HPLP.vcf

# Clean
rm -f outputs/pic_wing_${RUN_NAME}*
