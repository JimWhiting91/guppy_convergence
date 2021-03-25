# Dsuite

## Installation
Dsuite software is available at: https://github.com/millanek/Dsuite . There is also a link here to an excellent tutorial

## Summary of scripts
#### 01_Run_Dsuite_per_chrom_slurm.sh	
Subsets the VCF by chromosome and runs Dsuite over specific chromosomes. Commented out are options to explore sliding window analyses as well over certain trios.

#### 02_combine_per_chrom_slurm.sh
Combines per-chromosome files using `Dtrioscombine`

#### 03_calculate_fbranch_and_plot.sh
Takes final files and calculates Fbranch statistics and Z-scores. Produces Figure 2.