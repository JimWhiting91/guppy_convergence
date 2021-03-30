# Analysis of population structure outputs and introgression

## Summary of Rscripts
#### LP_LP_SFS_plot.R	
Plots 2dsfs for each LP-LP sfs (available in `data/LP_LP_comparisons`)

#### manuscript_structure_figures.R	
Main script for producing Figure 1. Map is produced based on shapefiles available in `data/trinidad_tobago` for Figure 1A. The results from SNAPP (`data/SNAPP_res`) are visualised for Figure 1B and summarised for divergence time estimates. FineSTRUCTURE outputs (`data/fineSTRUCTURE_res`) are read in and plotted in a heatmap. Plink-derived PCA are read in for all populations (`data/*.eigenval`, `data/*.eigenvec`) and plotted. Within-river SFS are plotted (`data/*sfs.txt`)

#### plot_map_with_rivers.R
Uses same inputs as above but with additional elevation data to plot rivers for supp figure. Elevation data downloaded from SRTM at https://earthexplorer.usgs.gov/

#### remove_caroni_pops_pca.R
Plots PCA with 1 of the 3 caroni rivers dropped for supp figure.