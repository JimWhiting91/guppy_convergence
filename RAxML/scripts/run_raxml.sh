#!/bin/bash
#PBS -q pq
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N run_raxml
#PBS -e logs/run_raxml.err.txt
#PBS -o logs/run_raxml.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

# We need these
module load flex/2.6.4-GCCcore-8.3.0 Bison/3.3.2

# Set up variables
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/raxml
REGION=000094F_0_500000_1100002
DATASET=five_aside_STAR
INPUT_FASTA=$MASTER/data/haplotypes/$DATASET/${REGION}/ALL_HAPLOID_${DATASET}_${REGION}.fa
SEED=2
OUTPUT=$MASTER/outputs/20_07_01_scaf94_CLAP_region_haplogroups
MODEL=GTR+G # Determine this first by running modeltest-ng

######## BEFORE RUNNING, RUN THIS LINE ON COMMAND-LINE TO CHECK INPUT, MAY HAVE TO SWITCH TO A REDUCED ALIGNMENT INPUT
# Parse input and get thread count
raxml-ng --parse --msa $INPUT_FASTA --model $MODEL --prefix ${OUTPUT}_T1
THREAD_N=$(grep "threads /" ${OUTPUT}_T1.raxml.log | cut -d ":" -f 2 | sed 's/ //g')

# Do we need to change the input file?
#if [ $(grep "Reduced" ${OUTPUT}_T1.raxml.log | wc -l) -eq 1 ]
#then
#  INPUT_FASTA=${OUTPUT}_T1.raxml.rba
#fi

# Infer the tree
raxml-ng --msa $INPUT_FASTA --model $MODEL --prefix ${OUTPUT}_T2 --seed $SEED --threads $THREAD_N

# Bootstrap and Convergence - Write as a loop and keep adding trees until we converge...
CONVERGE=NO
BOOTSTRAP_N=1
rm -f ${OUTPUT}_allbootstraps
while [ $CONVERGE = "NO" ]
do

# Change the seed
SEED=$(( $SEED + 1 ))
BOOTSTRAP_N=$(( $BOOTSTRAP_N + 1 ))

# Bootstrap
raxml-ng --bootstrap --msa $INPUT_FASTA --model $MODEL --prefix ${OUTPUT}_boot$BOOTSTRAP_N --seed $SEED --bs-trees 200 --threads $THREAD_N


# Merge
cat ${OUTPUT}_allbootstraps ${OUTPUT}_boot${BOOTSTRAP_N}.raxml.bootstraps >> ${OUTPUT}_allbootstraps

# Check for Bootstrap Convergence with 3% cut-off
raxml-ng --bsconverge --bs-trees ${OUTPUT}_allbootstraps --prefix ${OUTPUT}_converge --seed $SEED --bs-cutoff 0.03 --threads $THREAD_N

# Re-assign failsafe
FAILSAFE=$(cat ${OUTPUT}_converge.raxml.log | grep "did not" | wc -l)
if [ $FAILSAFE -ne 1 ]
then
CONVERGE=YES
fi

done

# Calculate bootstrap supports for nodes
raxml-ng --support --tree ${OUTPUT}_T2.raxml.bestTree --bs-trees ${OUTPUT}_allbootstraps --prefix ${OUTPUT}_TFINAL --threads $THREAD_N
