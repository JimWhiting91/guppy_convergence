# SNAPP

## Installation
SNAPP is available as an add-on package for BEAST2. For producing the XML and implementing the clock analysis, it is necessary to follow this workflow: https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md and from here download the prep software `snapp_prep.rb`

The constraints file used is `data/five_aside_STAR_SNAPP_wingei_HPLP_constraints.txt`

## Summary of scripts
#### 01_Prepare_VCF_from_outgroup_batch_call.sh
Takes the VCF input and filters based on % missing data per individual, thins for linkage, and calls `R/find_top_snps_for_vcf.R` to assess which SNPs are most informative for populations. Results in a VCF of 1,000 SNPs with 6 individuals per river/wingei

#### 02_prep_snapp_xml_and_run.sh
Implements `snapp_prep.rb` and runs SNAPP. Set up to run 3 separate iterations with different starting seeds for 500,000 MCMC.

#### 03_run_MCMC_futher.sh
If it is necessary to run further MCMC, this script takes the same inputs as above and appends results.
