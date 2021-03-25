# fineSTRUCTURE

## Installation
fineSTRUCTURE and chromopainter are available here: https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure_info.html

Input data needs to be phased.

We also do PCA here using plink

## Summary of scripts
#### 01_plink_PCA			
Takes VCF as input, prunes for linkage, and performs PCA over all individuals. The variant `01a_plink_PCA_drop_Caroni_pops_slurm.sh` is the same but sequentially removes individuals associated with the three Caroni rivers to investigate sampling artefacts.

#### 02_Convert_shapeit2chromo_slurm.sh	
Converts shapeit2 formatted phased outputs to finestructure format.

#### 03_run_fineSTRUCTURE.sh			
This script is formatted as an HPC submission script, but in fact it should be run line by line on the command-line. Submission to the HPC is handled by the `qsub_run*.sh` scripts, this script merely formats the pipeline and makes/moves/edits inputs and outputs. For example, you can run up to stage 1, then wait for those jobs to finish, before running the next set of lines for stage 2, etc etc

#### qsub_run.sh and qsub_run_highmem.sh
Modified versions of the packaged qsub_run.sh for Exeter PBS system
