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
#SBATCH --array=1-33
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk

# We need these...
module load BCFtools/1.9-foss-2018b

# Set up environment
MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/dsuite
DATASET=five_aside_STAR_v3_wingei_Guan_Tac_sisters_maf01_invariant_filtered_swap_TAC

cd $MASTER

mkdir outputs/$DATASET

CHR=$(awk '$2 > 1000000 {print $0}' ~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai | cut -f1 | sed "${SLURM_ARRAY_TASK_ID}q;d")
VCF_outgroup=~/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_v3_picta_wingei_maf01_noINV_4829912.no_picta.merged.vcf.gz
POPMAP_outgroup=$MASTER/data/five_aside_STAR_wingei.popmap
TREE=$MASTER/data/five_aside_STAR_individuals_tree_outgroup_rooted_swap_TAC.nwk

# Make temporary vcfs
bcftools view -r $CHR $VCF_outgroup > outputs/${DATASET}_${CHR}_outgroup.vcf

# Run first pass - with tree
Dsuite Dtrios outputs/${DATASET}_${CHR}_outgroup.vcf $POPMAP_outgroup -n ${DATASET}_${CHR} -t $TREE
# Run first pass - without tree
# #Dsuite Dtrios outputs/${DATASET}_${CHR}.vcf $POPMAP -n ${DATASET}_${CHR}
mv data/*${DATASET}_${CHR}_*txt outputs/$DATASET

# Run first pass - without tree and quartets
#Dsuite Dquartets outputs/${DATASET}_${CHR}.vcf $POPMAP -n ${DATASET}_${CHR}_quartets -t $TREE
#mv data/*${DATASET}_${CHR}_*txt outputs/$DATASET

# Loop over some chrs to test...
# #CHRS=(chr15 chr20 chr23 000119F_0 000094F_0)
# #CHRS=(chr19)
#
# for CHR in "${CHRS[@]}"
# do

# pops1=(TULP GH TUHP)
# pops2=(TUHP TULP TULP)
# pops3=(OH OH OH)
#
# WINDSIZE=50
# STEPSIZE=10
#
# for i in {0..2}
# do
#
# mkdir outputs/$DATASET/${CHR}_tmp
# cd outputs/$DATASET/${CHR}_tmp
#
# # Run the Test trio
# P1=${pops1[i]}
# P2=${pops2[i]}
# P3=${pops3[i]}
#
# echo -e "${P1}\t${P2}\t${P3}" > test_trios_${i}.txt
#
# # Run in sliding windows along genome, test_trios here is "APHP OL OH"
# Dsuite Dinvestigate -w $WINDSIZE,$STEPSIZE $MASTER/outputs/${DATASET}_${CHR}_outgroup.vcf $POPMAP_outgroup test_trios_${i}.txt
# mv ${P1}_${P2}_${P3}_localFstats__${WINDSIZE}_${STEPSIZE}.txt $MASTER/outputs/$DATASET/${P1}_${P2}_${P3}_trio_${CHR}_${DATASET}_localFstats_${WINDSIZE},${STEPSIZE}.txt
#
# done

cd $MASTER
rm -rf outputs/$DATASET/${CHR}_tmp

# Tidy up
rm -f $MASTER/outputs/${DATASET}_${CHR}.vcf $MASTER/outputs/${DATASET}_${CHR}_outgroup.vcf
