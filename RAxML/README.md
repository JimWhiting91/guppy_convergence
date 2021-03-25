# RAxML-NG

## Installation
RAxML-NG is available here: https://github.com/amkozlov/raxml-ng
A tutorial is available here: https://github.com/amkozlov/raxml-ng/wiki/Tutorial

## Summary of scripts
#### Extract_fasta_from_phase_one_haploid.sh
Extracts haploid phases from individuals in the phased VCF and formats for RAxML. For the analyses here, individuals were from homozygote haplogroups, so a random haplotype is chosen for each homozygote individual at the CL-AP region.

#### run_raxml.sh
Runs RAxML over the input fasta's and performs bootstrapping. Bootstrapping takes place in a while loop that periodically checks the log files for convergence. Once convergence has taken place, node support is calculated.