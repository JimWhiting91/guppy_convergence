# Population genomics of convergent evolution in high- and low-predation guppies
### Associated data and scripts for "Drainage-structuring of ancestral variation and a common functional pathway shape limited genomic convergence in natural high- and low-predation guppies"
#### Contact jwhiting2315@gmail.com with any issues

#### Raw data and SNP-calling
Raw data for this project is available at ENA: PRJEB43917 (Aripo, Madamas, Tacarigua) and PRJEB36704 (Guanapo and Oropouche). SNP calling was performed using the pipeline described here: https://github.com/josieparis/gatk-snp-calling

#### Five Aside STAR
Many scripts reference `five_aside_STAR`, this is simply the shorthand name given to the main dataset (five HP and five LP populations, STAR = shorthand name for the genome assembly used).

#### VCF Data
VCF data is available from FigShare, doi:10.6084/m9.figshare.14315771, store in `VCFs/`. The main VCF here is used in the majority of analyses, and this version of the VCF has been phased by beagle/shapeit. For introgression analyses, we included additional available data for Poecilia wingei (see manuscript). Smoove-called SV VCFs are located in `SV_and_coverage/data/`.

#### General data and functions
General functions can be found in `R/`, these are sourced in various R scripts and the path may need updating to reflect their current path. These include important functions such as the liftover of scaffold 000094F_0 positions to the updated chr20.

`general_data/` includes the fasta index for the male guppy genome, which is used by various scripts.

`popmaps/` includes a popmap for all 10 populations + wingei individuals, popmaps for each of the five rivers, and a total popmap of all individuals in the main VCF.

#### General Usage
Each subdirectory includes an environment within which to perform specific analyses. For example, `BayPass/` includes scripts required to run BayPass aux model and calculate XtX within rivers starting with the VCF. Each analysis subdirectory includes a `data/`, `outputs/`, `figs/`, `R/`, `scripts/`, `scripts/logs/`, `tables/` subdirectory with scripts based on relative paths assuming this structure. If any of these are not present, they should be made prior to running with `mkdir`. Some paths may need to be updated before scripts will run, although in most cases this should be case of resetting the `MASTER` path directory at the top of scripts as the majority of paths are relative.

The following directories include individual analyses, i.e. they are not dependent on each other to run:
```
dsuite/
fineSTRUCTURE/
BayPass/
RAxML/
SNAPP/
SV_and_coverage/
localPCA/
phasing/
```

Note: The scripts associated with fastsimcoal and SFS production were provided by Vitor Sousa and so are not presented here.

The three directories with `*_analysis` names should be run at the end, as these include scripts related to the three main sections of the manuscript and require as inputs the outputs from various analyses.

