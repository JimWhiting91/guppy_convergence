#!/bin/bash
module load BCFtools/1.9-foss-2018b

REFERENCE=~/guppy_research/STAR/STAR.chromosomes.release.fasta
DATASET=five_aside_STAR
PHASED_VCFS=~/HP_LP_trials/phasing/phased_vcfs/five_aside_STAR_v3_wingei_noMAF_lowNe

cd $PHASED_VCFS

# Read in an array for chr intervals based on the ref genome fasta index
chr_array=($(awk '{print $1}' ${REFERENCE}.fai))

# Create an array of all gvcf names
for i in "${chr_array[@]}"
do
echo "$PHASED_VCFS/${DATASET}_${i}_shapeit_beagle.vcf.gz" >> batch_inputs.txt
done

# Remove empties
find $PHASED_VCFS/ -type f -empty -delete

# Remove empty VCFs
for file in $(cat batch_inputs.txt | grep -v "chr")
do

SNP_N=$(zcat $file | grep -v "#" | wc -l)
if [ $SNP_N -eq 0 ]
then
rm -f $file
fi

done

# Filter list
ls $PHASED_VCFS/${DATASET}_*_shapeit_beagle.vcf.gz >> batch_filter.txt
grep -Fxf batch_filter.txt batch_inputs.txt > batch_inputs_2.txt

rm -f batch_filter.txt batch_inputs.txt

# We have to add the contig ID back in...
for i in "${chr_array[@]}"
do
length=$(grep -w $i ${REFERENCE}.fai | cut -f2)
zcat ${DATASET}_${i}_shapeit_beagle.vcf.gz | sed "4i\##contig=<ID=${i},length=${length}>" > ${DATASET}_${i}_shapeit_beagle_contigID.vcf
done

# Amend
sed -i 's/.vcf.gz/_contigID.vcf/g' batch_inputs_2.txt

# Concat together
bcftools concat -o $PHASED_VCFS/${DATASET}_allchr_shapeit_beagle.vcf -f batch_inputs_2.txt
rm -f batch_inputs_2.txt

# Tidy up
~/bin/bgziptabix.sh $PHASED_VCFS/${DATASET}_allchr_shapeit_beagle.vcf
rm -f $PHASED_VCFS/*.vcf
