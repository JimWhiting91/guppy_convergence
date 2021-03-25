# Local PCA

## Installation
The lostruct R package is available here: https://github.com/petrelharp/local_pca

## Summary of scripts
#### R/run_localPCA_set_region.R
This R script loops over all 10 populations and runs local PCA over a user-defined chromosome. Can redefine window size (set at 100 SNPs) if desired. Requires BCFtools in order to subset the main VCF by chromosome prior to reading into R, else the VCF can be subsetted outside of R and read in on a per chromosome basis. The analysis shouldn't be run across chromosomes, only within.