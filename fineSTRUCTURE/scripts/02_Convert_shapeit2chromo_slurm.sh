#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -p highmem
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=conversion_for_fineSTRUC
#SBATCH --error=logs/conversion_for_fineSTRUCTURE.err.txt
#SBATCH --output=logs/conversion_for_fineSTRUCTURE.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END 
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

# This assumes that all separate chromosomes have been phased with shapeit and are in a directory somewhere...

module load Perl/5.28.0-GCCcore-7.3.0

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/fineSTRUCTURE/
DATASET=scaf94_distal
PHASE_DIR=$MASTER/data/${DATASET}_phase

# Phases can be cat together to give full chromosome
rm -f $MASTER/data/${DATASET}_allchr_shapeit_beagle_phase.haps
for chr in $(cut -f1 ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | grep -v "alt")
do
cat $PHASE_DIR/*_${chr}_*.haps >> $MASTER/data/${DATASET}_allchr_shapeit_beagle_phase.haps
done

# Now convert to chromopainter
perl ~/software/fs_4.0.1/impute2chromopainter.pl -J $MASTER/data/${DATASET}_allchr_shapeit_beagle_phase.haps $MASTER/data/${DATASET}_chromoIN

# Make rec file
perl /gpfs/ts0/home/jw962/software/fs_4.0.1/makeuniformrecfile.pl $MASTER/data/${DATASET}_chromoIN.phase $MASTER/data/${DATASET}_chromoIN.rec
