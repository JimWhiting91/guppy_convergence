#!/bin/bash

# Take our output, calc fbranch and plot supp fig
MASTER=~/Exeter/five_aside/dsuite
cd $MASTER

# Set
DATASET=five_aside_STAR_v3_wingei_Guan_Tac_sisters_maf01_invariant_filtered_swap_TAC
TREE=data/five_aside_STAR_individuals_tree_outgroup_rooted_swap_TAC.nwk
TREE_RES=data/$DATASET/${DATASET}_combined_trios_${DATASET}_combined_trios_combined_tree.txt

# Get fbranch
Dsuite Fbranch -p 0.001 $TREE $TREE_RES > outputs/${DATASET}_fbranch_res
Dsuite Fbranch -p 0.001 -Z test $TREE $TREE_RES > outputs/${DATASET}_fbranch_res_with_Z

# Plot
python3 ~/software/Dsuite/utils/dtools.py -n figs/${DATASET}_fbranch --outgroup Outgroup outputs/${DATASET}_fbranch_res $TREE
