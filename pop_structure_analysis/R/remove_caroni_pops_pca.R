# Compare PCA when only including 2 Caroni Populations
library(ggplot2)
source("~/Exeter/five_aside/five_aside_functions.R")

# Make metadata
river_meta <- data.frame(pop=c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD"),
                         river=rep(c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"),each=2),
                         pred=rep(c("HP","LP"),5))

# Read in each of the inputs and plot
pops <- c("Aripo","Guanapo","Tacarigua")
pop_dropped_pca <- lapply(pops,function(pop){
  
  # Fetch files
  vals <- read.table(paste0("data/five_aside_STAR_No",pop,"_plink_out_NoLD_PCA_No",pop,".eigenval"))
  vecs <- data.frame(read.table(paste0("data/five_aside_STAR_No",pop,"_plink_out_NoLD_PCA_No",pop,".eigenvec")))
  colnames(vecs) <- c("ind1","ind2",paste0("PC",1:(ncol(vecs)-2)))
  
  # Calculate percent
  vals <- round((vals/sum(vals)) * 100,1)
  
  # Add in metadata
  vecs$river <- NA
  vecs$pred <- NA 
  for(i in river_meta$pop){
    vecs[grep(i,vecs$ind1),"river"] <- river_meta[river_meta$pop==i,"river"]
    vecs[grep(i,vecs$ind1),"pred"] <- river_meta[river_meta$pop==i,"pred"]
  }
  vecs$river_F <- factor(vecs$river,levels = five_aside_colour_rivers$river[five_aside_colour_rivers$river != pop])
  
  # # And plot
  # ggplot(vecs,aes(PC1,PC2,colour=river,shape=pred))+
  #   geom_point()+
  #   stat_ellipse(level=0.95)+
  #   scale_colour_manual(breaks = five_aside_colour_rivers$river[five_aside_colour_rivers$river != pop],
  #                       values = five_aside_colour_rivers$colour[five_aside_colour_rivers$river != pop])

  
  # Make plot
  ggplot(vecs,aes(PC1,PC2,colour=river_F))+
    geom_jitter(size=4,aes(shape=pred),alpha=0.5)+
    stat_ellipse(show.legend = F)+
    scale_colour_manual(values=five_aside_colour_rivers$colour,
                        breaks=five_aside_colour_rivers$river_F)+
    scale_shape_manual(values=c(19,17))+
    labs(x=paste0("PC1 (",vals$V1[1],"%)"),y=paste0("PC2 (",vals$V1[2],"%)"),shape="",colour="")+
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=16),
          legend.position = "right",
          legend.text = element_text(size=12),
          panel.grid = element_blank(),
          title = element_text(size=16))+
    ggtitle(paste0(pop, " removed"))
})

# Plot together
library(cowplot)
pdf("figs/FigureSX_PCA_with_caroni_pops_dropped.pdf",height=12,width=5)
plot_grid(plotlist = pop_dropped_pca,
          ncol=1)
dev.off()
