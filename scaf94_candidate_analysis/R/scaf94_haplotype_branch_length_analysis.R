#####################################################################
# Analysis of local structure and branch lengths along chr20
lib<-c("ggridges","cowplot","ggplot2","data.table","viridis","vcfR","parallel","ape","ggtree","phytools","parallel","poppr","phyloch","adephylo")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in vcf
vcf <- read.vcfR("~/Exeter/VCFs/five_aside_STAR_chr20_scaf94_merged_shapeit_beagle.vcf.gz")

# Get windows
wind_size = 100000
winds <- seq(0,max(as.integer(vcf@fix[,2])),by=wind_size)
winds2 <- winds+wind_size
winds2[length(winds2)] <- max(as.integer(vcf@fix[,2]))

# Build the popmap
pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
pop_names <- c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
popmap <- data.frame(ind=colnames(vcf@gt)[2:ncol(vcf@gt)],
                     pop=NA)
for(i in 1:length(pops)){
  popmap[grep(pops[i],popmap$ind),"pop"] <- pop_names[i]
}

# Loop over windows and build NJ trees
window_trees <- lapply(1:length(winds),function(x){
  
  gen_dd <- vcfR2genind(vcf[as.integer(vcf@fix[,2]) > winds[x] &
                              as.integer(vcf@fix[,2]) <= winds2[x],])
  
  ######################################################################
  # Build tree and bootstrap
  pop_tree <- aboot(gen_dd, tree = "bionj", distance = nei.dist, sample = 1, showtree = F, cutoff = 0, quiet = F,threads = 6,root = FALSE)
  branch_lengths <- reshape2::melt(as.matrix(distTips(pop_tree)))
  
  # Fill the matrix
  comps_to_make <- combn(1:10,2)
  dist_mat <- matrix(ncol=10,nrow=10)
  colnames(dist_mat) <- pop_names
  rownames(dist_mat) <- pop_names
  for(i in 1:ncol(comps_to_make)){
    tmp1 <- popmap[grep(pop_names[comps_to_make[1,i]],popmap[,"pop"]),"ind"]
    tmp2 <- popmap[grep(pop_names[comps_to_make[2,i]],popmap[,"pop"]),"ind"]
    
    # Fetch distances
    dist_mat[comps_to_make[1,i],comps_to_make[2,i]] <- mean(branch_lengths[branch_lengths$Var1 %in% tmp1 & branch_lengths$Var2 %in% tmp2,"value"])
  }
  
  # Melt output
  melted_out <- na.omit(reshape2::melt(dist_mat))
  melted_out$start <- winds[x]
  melted_out$end <- winds2[x]
  melted_out$mid <- mean(c(winds[x],winds2[x]))
  
  return(melted_out)
  
})

# Merge these and visualise
merged_trees <- data.frame(rbindlist(window_trees))
merged_trees$comp <- paste0(merged_trees$Var1,"_",merged_trees$Var2)
comps_to_keep <- c("APHP_APLP","GHP_APHP","TACHP_APHP","APHP_OHP","APHP_MADHP")
merged_trees2 <- merged_trees[merged_trees$comp %in% comps_to_keep,]

# Re-name for neatness
new_names <- c("APHP-APLP","APHP-GHP","APHP-TACHP","APHP-OHP","APHP-MADHP")
for(i in 1:length(new_names)){
  merged_trees2$comp <- gsub(comps_to_keep[i],new_names[i],merged_trees2$comp)
}
merged_trees2$comp_F <- factor(merged_trees2$comp,levels=(new_names))

# Make ridge colours
ridge_colours <- data.frame(comp_F2 = c("APHP-APLP","APHP-GHP","APHP-TACHP","APHP-OHP","APHP-MADHP"),
                            fill_col = five_aside_colour_rivers$colour[c(3,2,1,4,5)])

# Visualise distance 
chr20_plot <- ggplot(merged_trees2,aes(x=mid,y=value,colour=comp_F))+
  annotate("rect", xmin = 1633448, xmax = 2133448, ymin = -Inf, ymax = Inf,
           alpha = .75,fill = "grey75")+
  geom_line()+
  facet_wrap(~comp_F,ncol=1,strip.position = "top")+
  labs(y="Mean Branch Distance",x="Chr 20 Position (Mb)")+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=14),
        legend.position = "none")+
  scale_colour_manual(breaks = ridge_colours$comp_F2,
                    values = ridge_colours$fill_col)+
  scale_x_continuous(breaks=c(0,1e7,2e7),labels=c(0,10,20))+
  scale_y_continuous(breaks=c(0.0,0.4,0.8),labels=c(0,0.4,0.8))

# Find the mean distance within the focal region
mean_CL <- mean(merged_trees2[merged_trees2$mid > 1633448 &
                                merged_trees2$ mid < 2133448 &
                                merged_trees2$comp == "APHP-APLP","value"])

# And visualise as the ridge
merged_trees2$comp_F2 <- factor(merged_trees2$comp,levels=new_names)

ridges <- ggplot(merged_trees2,aes(x=value,y=comp_F2,fill=comp_F2))+
  geom_density_ridges2(quantile_lines = TRUE, quantiles = 2,alpha=0.5)+
  geom_vline(xintercept = mean_CL,linetype="dashed")+
  labs(x="Mean Branch Distance",y="Comparison")+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none")+
  scale_fill_manual(breaks = ridge_colours$comp_F2,
                    values = ridge_colours$fill_col)+
  facet_wrap(~comp_F2,ncol=1,scales="free_y",strip.position = "top")
  

# Also make a newick summary to accompany the branch distances
fileConn<-file("data/tree.nwk")
writeLines("((((APLP,APHP),GHP),TACHP),(OHP,MADHP));", fileConn)
close(fileConn)
tree <- read.tree("data/tree.nwk",)

# Simple for nodes
ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

# Main
tree_fig <- ggtree(tree)+
  #geom_text(aes(label=node))+
  geom_tiplab(size=6)+
  geom_hilight(node=10, fill=five_aside_colour_rivers$colour[3]) + 
  geom_hilight(node=3, fill=five_aside_colour_rivers$colour[2]) + 
  geom_hilight(node=4, fill=five_aside_colour_rivers$colour[1]) + 
  geom_hilight(node=6, fill=five_aside_colour_rivers$colour[5]) + 
  geom_hilight(node=5, fill=five_aside_colour_rivers$colour[4]) +
  ggplot2::xlim(0, 6) 

tree_fig <- flip(tree_fig, 5, 6) %>% flip(2, 1)

# plot togerther
full_branch_fig <- cowplot::plot_grid(tree_fig,
                   ridges,
                   chr20_plot,
                   ncol=3,rel_widths = c(2,2,3),align = "h",axis = "tblr",labels=c("D","E","F"),label_size=24)
saveRDS(full_branch_fig,"outputs/full_branch_fig.rds")


