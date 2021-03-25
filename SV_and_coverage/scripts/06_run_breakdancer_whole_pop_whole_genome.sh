#!/bin/bash
#PBS -q pq
#PBS -l walltime=00:20:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N run_breakdancer
#PBS -e logs/run_breakdancer.err.txt
#PBS -o logs/run_breakdancer.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

# We need these
module load GLib/2.60.1-GCCcore-8.2.0 SAMtools/1.9-foss-2018b
module load GDGraph/1.54-GCCcore-8.2.0-Perl-5.28.1

# Set variables
MASTER=~/HP_LP_trials/breakdancer/
BAM_DIR=~/GUPPYCon/STAR_bams/
CHR=chr20

# Make all the config files
pops=(Aripo Guanapo Tacarigua Oropouche Madamas)
for POP in "${pops[@]}"
do
#POP=Guanapo
#REGION=1-4000000
# For whatever bizarre reason, isca won't accept this so I currenetly just run it in-line
cd $MASTER

# This needs running else we can't use the perl script
cpan GD::Graph::histogram
cpanm GD::Graph::histogram
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

# Set the popmap
POPMAP=/gpfs/ts0/home/jw962/HP_LP_trials/popmaps/five_aside_STAR/${POP}.popmap
OUT=$MASTER/outputs/$POP

# Build our bam file list
rm -f outputs/${POP}_bams.txt
for ind in $(cat $POPMAP)
do
  ls $BAM_DIR | grep "${ind}_" | grep -v ".bam." >> outputs/${POP}_bams.txt
done

# Make the config file
cmd=""
for file in $(cat outputs/${POP}_bams.txt)
do
cmd+="${BAM_DIR}/$file "
done

perl ~/software/breakdancer/perl/bam2cfg.pl -g $cmd > data/${POP}_merged.config

done

# Run breakdancer, the r flag here is being treated as a fairly stringent maf filter of ~ 0.1, ie. support must be at least ~ 15%, 6 RGs. Although removing this in favour of later filtering.
pops=(Guanapo Tacarigua Oropouche Madamas)
for POP in "${pops[@]}"
do
mkdir outputs/$POP
awk '$2 > 1000000' ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | parallel -j 8 "breakdancer-max -h -o {} -g outputs/$POP/${POP}_{}_breakdancer.bed data/${POP}_merged.config > outputs/$POP/${POP}_{}_breakdancer.out"
done

# Tidy up
rm -f data/${POP}*_${CHR}*.bam
