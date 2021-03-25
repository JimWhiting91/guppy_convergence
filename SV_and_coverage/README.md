# SV (smoove + breakdancer) and Coverage (deeptools)

## Installation
smoove is available here: https://github.com/brentp/smoove

deeptools is available here: https://deeptools.readthedocs.io/en/develop/index.html

breakdancer is available here: http://breakdancer.sourceforge.net/

## Summary of scripts
#### 01a_bamCoverage_run_slurm.sh
Calculates coverage for individuals from the same population based on BAM files. Runs as an array where each jobID is a chromosome, and can parse population identifier as a trailing variable. For example, `sbatch 01a_bamCoverage_run_slurm.sh GH` will calculate coverage for all bam files where `ls path/to/bams/ | grep "GH"`. Script also makes blacklist regions at the starts and ends of chromosomes where coverage can be difficult to estimate due to repeats. Removing these regions improves normalisation. The variant `01b_bamCoverage_scf.sh` is used over scaffolds, and has a less stringest trimming of data from the start/end of the scaffold

#### 02_mash_individuals_to_population.sh
Uses datamash to average across individuals within a population

#### 03_coverage_to_bed.sh			
Converts coverage outputs (50bp calculated) to bed files of any window size

#### 04_Run_R_analysis_for_ratios.sh
Submits `R/01_compare_HP_LP_coverage.R` to the HPC to calculate coverage ratios between HP and LP populations from the same river.

#### 05_run_smoove_lumpy.sh
Calls structural variants based on BAM files and excluding repeat regions (masked with bed file). Calls SVs within river pairs (N=5).

#### 06_run_breakdancer_whole_pop_whole_genome.sh
Calls structural variants based on BAM files with Breakdancer. Loops over populations (based on river popmaps, as above) and runs over all chr/scaf > 1000000 bp per pop with GNU parallel.

