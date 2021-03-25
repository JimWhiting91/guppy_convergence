# Analysis of scaffold 94 candidate region

Pulls together local PCA results, raxml results, and performs branch length analysis to assess evolutionary history of scaffold 94 candidate.

Makes use of GenotypePlot package to visualise haplotype structure, available here: https://github.com/JimWhiting91/genotype_plot

## Summary of Rscripts
#### scaf94_haplotype_branch_length_analysis.R
Performs branch length analysis over chromosome 20.

#### scaf94_figure.R
Produces majority of Figure 4. Visualises haplotype structure over all individuals with GenotypePlot, plots localPCA results (`data/scaf94_figs/five_aside_STAR_*_chr20_scaf94_merge.txt`), plots the raxml tree (`data/scaf94_figs/20_07_01_scaf94_CLAP_region_haplogroups_TFINAL.raxml.support`), and plots results of branch length analysis.
