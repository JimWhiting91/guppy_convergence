#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=168:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=addon_five_aside_beast_poec
#SBATCH --error=addon_five_aside_beast_poec.err.txt
#SBATCH --output=addon_five_aside_beast_poec.out.txt
#SBATCH --export=All
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk
#SBATCH --array=1-3

# First make the input xml...
module load Ruby/2.5.0-intel-2017a

MASTER=/gpfs/ts0/home/jw962/HP_LP_trials/SNAPP
cd $MASTER

RUN_NAME=pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2_HPLP
VCF=$MASTER/data/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2.nomiss.thin50k.highFST.nochr12.HPLP.vcf
SAMPLES=$MASTER/data/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2_SAMPLES_HPLP.txt
CONSTRAINTS=$MASTER/data/five_aside_STAR_SNAPP_wingei_HPLP_constraints.txt

# Load the modules
module purge
module load Beast/2.6.3-GCC-9.3.0

# We pull the seed from the previous log file

# Run Beast
cd outputs/$RUN_NAME
beast -resume -threads 24 $MASTER/data/${RUN_NAME}_input_iter${SLURM_ARRAY_TASK_ID}.xml
