# BayPass

## Installation
Baypass software is available at: http://www1.montpellier.inra.fr/CBGP/software/baypass/

## Summary of scripts
#### 01_VCF_to_baypass.sh		
This script takes the VCF and loops over populations to produce population specific Baypass-formatted inputs. These can be combined to make the full input using the shell command `paste`. This can be run without submitting to an HPC with `vcf2baypass.sh`.

#### 02a_CovMatrix_Run.sh
Performs a linkage prune on the VCF, and retains rows from the full SNP input that are pruned for linkage. Then performs a core run of BayPass to produce a covariance matrix. This should be performed 10 times and the average covariance matrix used. `02b_CovMatrix_Run_per_pop.sh` is the same, but is run per population rather than over the whole dataset.

#### 03a_Run_baypass_subsets.sh
Runs BayPass over subsets of the full input. Assumes the data has already been subsetted by `R/subset_baypass_genos.R`. There are options in here to run either the xtx model or the aux model for HP-LP association.

#### 04_Analyse_baypass_XTX.sh
Submits R/get_xtx_significance.R to the HPC to calculate significance of XtX based on simulations. 
