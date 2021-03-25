# Phasing
This repo takes the VCF and BAMs and phases using Beagle v5 (https://faculty.washington.edu/browning/beagle/beagle.html) and shapeit v2 (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html).

## Summary of Scripts
#### 01_Phase_with_beagle_then_slurm.shapeit_slurm.sh
Takes the VCF with all individuals and first performs imputation and phasing with beagle. Beagle-phased VCF is then parsed to shapeit, where phase-informative read information derived from the originla BAM files is used to re-phase the data. Script runs as an HPC array, and runs over individual chromosomes/scaffolds.

#### 02_concat_phased_vcfs.sh
Concatenates all of the per chromosome phased VCFs into a single VCF.