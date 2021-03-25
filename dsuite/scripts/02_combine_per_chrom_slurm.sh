#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=dsuite_run
#SBATCH --error=logs/dsuite_run.err.txt
#SBATCH --output=logs/dsuite_run.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

# We need these...
module load BCFtools/1.9-foss-2018b

# Set up environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/dsuite
DATASET=five_aside_STAR_v3_wingei_Guan_Tac_sisters_maf01_invariant_filtered_swap_TAC
TREE=$MASTER/data/five_aside_STAR_individuals_tree_outgroup_rooted_swap_TAC.nwk

cd $MASTER
mkdir outputs/$DATASET

# Move to output dir
cd outputs/$DATASET
inputs=""
for i in {1..23}
do
  CHR=chr$i
  echo $CHR
  if [[ "$CHR" != "chr12" ]]
then
  input_tmp=$(ls *Dmin.txt | grep -v "quartet" | grep "_${CHR}_" | sed 's/_Dmin.txt//g')
  inputs+="$input_tmp "
fi
done

# Loop to make input
Dsuite DtriosCombine $inputs -n ${DATASET}_combined_trios -t $TREE -o ${DATASET}_combined_trios
