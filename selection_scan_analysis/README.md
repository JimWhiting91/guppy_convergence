# Selection Scanning and Analysis
This repo includes scripts to calculate absolute allele frequency difference (AFD) and XP-EHH with selscan. Selscan is available here:

Rscripts include processing of outputs, overlapping outlier analysis, summary of selection, candidate gene identification, and functional enrichment analyses.

## Summary of Scripts
#### run_selscan_xpehh.sh
Calculates XP-EHH from phased VCF data. Scripts runs over a river, but produces inputs for HP and LP pops within the river for xp-ehh to be calculated between them. Calculates on a per chromosome basis, so each chromosome can be submitted to HPC separately.

## Summary of Rscripts
#### R/calculate_AFD_five_aside.R 
Uses `calc_AF()` from the global `R/five_aside_functions.R` to calculate the allele frequency changes among HP and LP populations. The inputs for this need to be allele frequency files (2 columns per pop, row per SNP), which is the same format as the input for BayPass, so the `BayPass/scripts/vcf2baypass.sh` script can be used to make this input. For e.g. the allele counts for all 10 pops looks like:
```
3	19	4	16	24	8	36	0	20	16	17	17	35	1	35	5	40 0	32	0
22	0	24	0	19	13	14	22	36	2	30	4	36	0	40	0	40 0	32	0
4	18	5	19	28	8	34	0	18	18	17	15	37	1	34	6	40 0	32	0
7	15	6	22	31	3	36	0	25	13	19	15	37	1	38	2	40 0	32	0
18	4	26	0	35	1	36	0	35	3	32	6	27	11	24	16	20 20	16	16
4	18	0	26	33	5	34	0	17	19	8	22	4	34	4	36	0 40	0	32
1	21	1	21	30	0	32	0	24	0	26	0	38	0	38	0	38 0	30	0
1	19	0	20	27	3	17	17	18	4	21	7	7	27	0	40	0 38	0	30
1	19	0	22	27	9	15	15	8	20	6	20	4	32	1	37	3 35	1	29
20	0	24	0	34	0	30	0	24	4	22	4	23	13	25	11	15 23	21	9
```

#### R/process_five_aside_STAR_XTX.R
Processes the XtX outputs and calculates 10kb windows

#### R/process_baypassAUX_perSNP.R
Processes the BayPass auxmodel HP-LP association outputs, for e.g. combines back together the runs from the separate subsets and makes windows and baypass result figures (Figure 3B-D). The variant of this `process_baypassAUX_perSNP_noMerge_chr20_scaf94.R` does the same, but does not merge scaffold 94 and chromosome 20 as described in the manuscript. Also produces supp figure for scaffold 94 region with coverage comparisons.

#### R/process_xpehh_perSNP.R
Processes the XP-EHH outputs from selscan (as described above) and turns to windows.

#### R/compare_all_selection_results.R
Compares and contrasts the processed outputs from each of the selection scan analyses and BayPass outputs. Looks for overlap among within-river selection, and also semi-overlap (where windows can 'almost' overlap). Produces Figure 3A and various summaries of outliers.

#### R/outlier_enrichment_analysis.R
Performs functional enrichment analysis within each river, and also assesses significance of cadherin signaling pathway enrichment in all rivers by permutations. Functional enrichment is based on zebrafish orthologues, which are identified here by pulling candidate windows out of the male guppy genome, aligning them to the female guppy genome (which has Ensembl annotations), and using those aligned, annotated regions to call zebrafish orthologues. Permutations make use of Panther outputs, which are provided in `outputs/pather_results_*` 

#### R/10kb_window_SNP_count_permutations.R
Uses counts of SNPs per window (counted during `calculate_AFD_five_aside.R`) to assess whether outlier windows have low SNP counts. Performs permutations against the genome-wide distribution of SNP count per 10kb window and compares observed median to permuted median for supp figs.
