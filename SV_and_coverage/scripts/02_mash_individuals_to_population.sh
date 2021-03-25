#!/bin/bash
#PBS -q pq  
#PBS -l walltime=00:60:00
#PBS -l nodes=1:ppn=12
#PBS -A Research_Project-T110748
#PBS -N coverage_mash
#PBS -e logs/coverage_mash.err.txt 
#PBS -o logs/coverage_mash.out.txt 
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

POP=$1
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/deeptools/outputs/coverage/
DATASET=five_aside_STAR

mkdir $MASTER/$DATASET/$POP/pop_avg_filtered

cd $MASTER/$DATASET/$POP/

parallel -j 12 "cat *.chr{}.* | datamash -H -s -g 2 median 4 | sort -nk1 > ${POP}_chr{}_mashed.bed" ::: {1..23}

for chrN in {1..23}
do
CHR=chr${chrN}
awk -v chr="$CHR" '{print chr,"\t",$0}' $MASTER/$DATASET/$POP/${POP}_${CHR}_mashed.bed > $MASTER/$POP/${POP}_chr${chrN}_mashed2.bed
done

rm -f $MASTER/$POP/*_mashed.bed

mv $MASTER/$POP/*mashed* $MASTER/$DATASET/$POP/pop_avg_filtered

# Simpler loop for scaffolds

for scaf in $(cat /gpfs/ts0/home/jw962/HP_LP_trials/deeptools/data/STAR.chromosomes.release.fasta.fai | cut -f1 | grep "000")
do
cat *$scaf* | datamash -H -s -g 2 median 4 | sort -nk1 > ${POP}_${scaf}_mashed.bed
awk -v chr="$scaf" '{print chr,"\t",$0}' $MASTER/$DATASET/$POP/${POP}_${scaf}_mashed.bed > $MASTER/$DATASET/$POP/${POP}_${scaf}_mashed2.bed
done
rm -f $MASTER/$DATASET/$POP/*_mashed.bed
mv $MASTER/$DATASET/$POP/*mashed* $MASTER/$DATASET/$POP/pop_avg_filtered
done

