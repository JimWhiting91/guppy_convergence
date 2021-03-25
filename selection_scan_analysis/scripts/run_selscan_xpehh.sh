#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=16 # nodes=number of nodes required. ppn=number of processors per node
#PBS -A Research_Project-T110748
#PBS -q pq
#PBS -N run_selscan
#PBS -e logs/run_selscan.err.txt
#PBS -o logs/run_selscan.out.txt
#PBS -S /bin/bash
#PBS -V
#PBS -d .
#PBS -m e -M j.whiting2@exeter.ac.uk

# Script runs Selscan for iHS

# We need these...
module load BCFtools/1.9-foss-2018b

# Feed population from the command line

# Set up the environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/selscan
DATASET=five_aside_STAR
HAP_DIR=~/HP_LP_trials/phasing/phased_vcfs/$DATASET/
POP=$1
CPU=16
OUT=05_13_v3

# Go to working
cd $MASTER

# Source popmap and make input VCF
POPMAP_DIR=~/HP_LP_trials/popmaps/$DATASET

# Loop over CHR
for CHR in $(cat ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | awk '$2 > 1000000' | cut -f1)
do

  echo "STARTING $CHR"

  # Just in case...
  tabix -f -p vcf $HAP_DIR/${DATASET}_${CHR}_shapeit_beagle.vcf.gz

    # For each population we need to fetch the popmap
  for POPMAP in $(ls $POPMAP_DIR/ | grep $POP)
  do

    pop=$(echo $POPMAP | sed 's/.popmap//g')

    bcftools view -r $CHR -S $POPMAP_DIR/$POPMAP $HAP_DIR/${DATASET}_${CHR}_shapeit_beagle.vcf.gz > data/${pop}_${CHR}_phased.vcf
    bcftools annotate --set-id +'%CHROM\-%POS' data/${pop}_${CHR}_phased.vcf > data/${pop}_${CHR}_phased_ID.vcf
    rm -f data/${pop}_${CHR}_phased.vcf

  done

  # Make sure we have HP and LP
  HP=$( ls data/${POP}*_${CHR}_phased_ID.vcf | grep "H")
  LP=$( ls data/${POP}*_${CHR}_phased_ID.vcf | grep "L")

  # Make the map file
  grep -v "#" $HP | awk -v OFS='\t' '{print $1,$3,$2,$2}' > data/${POP}_${CHR}_phased_ID.vcf.mapfile

  # Run the file - can consider adding a maf?
  selscan --xpehh --threads $CPU --vcf $LP --vcf-ref $HP --map data/${POP}_${CHR}_phased_ID.vcf.mapfile --out outputs/${OUT}_${DATASET}_${POP}_${CHR}

done

# And normalise
norm --xpehh --crit-percent 0.05 --winsize 50000 --files outputs/${OUT}_${DATASET}_${POP}_*.xpehh.out
# And by window
norm --xpehh --crit-percent 0.05 --winsize 50000 --bp-win --files outputs/${OUT}_${DATASET}_${POP}_*.xpehh.out
#norm --pi --crit-percent 0.05 --winsize $win --bp-win --files outputs/${DATASET}_${POP}_*.pi.out
