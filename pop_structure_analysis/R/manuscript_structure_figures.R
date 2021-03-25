####################################################
# Assemble Plots for main five aside manuscript figure 1

lib<-c("cowplot","ggplot2","data.table","viridis","ggtree","ape","wesanderson","patchwork")
lapply(lib,library,character.only=T)

# Get our functions
source("~/Exeter/five_aside/five_aside_functions.R")

##### FineStructure Matrix ######
# Read in results
chunks<-read.table("data/fineSTRUCTURE_res/five_aside_STAR_fineSTRUCTURE_linked.chunklengths.out",header=T)

# Remove Second Column
chunks<-chunks[,c(1,3:ncol(chunks))]
melt_chunks<-melt(chunks)
colnames(melt_chunks)[2]<-"Donor"

# Plot Order
inds<-as.character(unique(melt_chunks$Recipient))
ordered<-c(grep("LT",inds,value = T),
           grep("UT",inds,value = T),
           grep("GH",inds,value = T),
           grep("GL",inds,value = T),
           grep("APHP",inds,value = T),
           grep("APLP",inds,value = T),
           grep("LO",inds,value = T),
           grep("UQ",inds,value = T),
           grep("LMD",inds,value = T),
           grep("UMD",inds,value = T))

melt_chunks$Recipient_F<-factor(melt_chunks$Recipient,levels=ordered)
melt_chunks$Donor_F<-factor(melt_chunks$Donor,levels=ordered)

# Also get axis labels
labels<-data.frame(pop=c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP"),
                   true=c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD"))

# Find length of each
for(i in 1:10){
  labels$N[i]<-length(grep(labels$true[i],inds,value = T))
}

# Sum
labels$sum<-labels$N
for(i in 2:10){
  labels$sum[i]<-labels$N[i]+labels$sum[i-1]
}

# Mids
labels$mid<-labels$N/2
for(i in 2:10){
  labels$mid[i]<-mean(c(labels$sum[i],labels$sum[i-1]))
}

# Plot
FS_fig<-ggplot(melt_chunks,aes(x=Donor_F,y=Recipient_F,fill=log10(value)))+
  geom_tile()+
  scale_fill_viridis(option="magma")+
  theme_bw()+
  scale_x_discrete(breaks=ordered[round(labels$mid)],
                   labels=labels$pop)+
  scale_y_discrete(breaks=ordered[round(labels$mid)],
                   labels=labels$pop)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=3),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=30,hjust=1),
        axis.ticks = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.title = element_text(size=20),
        legend.text = element_text(size=14,angle=45,hjust=1))+
  labs(x="Donors",y="Recipients",fill=expression(log[10](Chunk~Length)))

##### PLINK PCA #######
# Read in Plink Results
pca_dd<-read.table("data/five_aside_STAR_plink_out_NoLD_PCA.eigenvec")
pca_vals<-read.table("data/five_aside_STAR_plink_out_NoLD_PCA.eigenval")
pca_vals$percent<-round((pca_vals$V1*100)/sum(pca_vals$V1),digits = 1)

# Plot PC1 and PC2
pca_dd$pop<-NA
for(i in 1:nrow(labels)){
  pca_dd[grep(labels$true[i],pca_dd$V1,"pop"),"pop"]<-as.character(labels$pop[i])
}

# Set rivers
pca_dd$river<-"Aripo"
pca_dd[pca_dd$pop %in% c("GHP","GLP"),"river"]<-"Guanapo"
pca_dd[pca_dd$pop %in% c("MADHP","MADLP"),"river"]<-"Madamas"
pca_dd[pca_dd$pop %in% c("OLP","OHP"),"river"]<-"Oropouche"
pca_dd[pca_dd$pop %in% c("TACHP","TACLP"),"river"]<-"Tacarigua"

# Set predation
pca_dd$pred<-"HP"
pca_dd[grep("LP",pca_dd$pop),"pred"]<-"LP"

pca_dd$river_F<-factor(pca_dd$river,levels = c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"))

# Make labels
to_label<-data.frame(labels=c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"),
                     x_pos=c(0.035,0.03,0.01,-0.03,-0.075),
                     y_pos=c(-0.125,0.11,0.02,-0.02,0.02))
to_label$river_F<-factor(to_label$labels,levels=to_label$labels)

# Plot
pca_fig<-ggplot(pca_dd,aes(V3,V4,colour=river_F))+
  geom_jitter(size=4,aes(shape=pred),alpha=0.5)+
  stat_ellipse(show.legend = F)+
  scale_colour_manual(values=five_aside_colour_rivers$colour,
                      breaks=five_aside_colour_rivers$river_F,guide=FALSE)+
  scale_shape_manual(values=c(19,17))+
  labs(x="PC1 (33.5%)",y="PC2 (19.8%)",shape="",colour="")+
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = c(0.05,0.95),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.text = element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 8)))+
  geom_label(data=to_label,aes(x=x_pos,y=y_pos,label=labels,fill=river_F),colour="white",fontface="bold",size=6)+
  scale_fill_manual(values=five_aside_colour_rivers$colour,
                    breaks=five_aside_colour_rivers$river_F,guide=FALSE)

# Also get subsets of PC1-PC2 for each individual river...
river_PCA_fig <- ggplot(pca_dd,aes(x=V3,y=V4,shape=pred,colour=river_F))+
    geom_point(size=4,alpha=0.75)+
    stat_ellipse(show.legend = F)+
    scale_colour_manual(values=five_aside_colour_rivers$colour,
                        breaks=five_aside_colour_rivers$river_F,guide=FALSE)+
    scale_shape_manual(values=c(19,17))+
    labs(x="PC1 (33.5%)",y="PC2 (19.8%)",shape="",colour="")+
    theme_bw()+
    theme(axis.text = element_blank(),
          axis.title = element_text(size=16),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.text = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid = element_blank(),
          axis.ticks = element_blank())+
    scale_fill_manual(values=five_aside_colour_rivers$colour,
                      breaks=five_aside_colour_rivers$river_F,guide=FALSE)+
    facet_wrap(~river_F,ncol=1,scales = "free",strip.position = "right")

pca_fig2<-ggplot(pca_dd,aes(V3,V5,colour=river_F))+
  geom_jitter(size=4,aes(shape=pred),alpha=0.5)+
  stat_ellipse(show.legend = F)+
  scale_colour_manual(values=five_aside_colour_rivers$colour,
                      breaks=five_aside_colour_rivers$river_F)+
  scale_shape_manual(values=c(19,17))+
  labs(x="PC1",y="PC3 (17.9%)",shape="",colour="")+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_text(size=14),
        legend.position = "none",
        legend.text = element_text(size=16),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

# PCA plot inset 
PCA_fig_inset <-
  ggdraw() +
  draw_plot(pca_fig) +
  draw_plot(pca_fig2, x = 0.2, y = .135, width = .4, height = .36)

# And combine with the per rivers...
full_pca_fig <- plot_grid(PCA_fig_inset,river_PCA_fig,
                          ncol=2,rel_widths = c(5,2),axis = "tblr",align = "h")

##### Plot SNAPP Tree ######
library(treeio)
bs_trees <- read.beast("data/SNAPP_res/downsample.tree")
snapp_tree <- read.beast("data/SNAPP_res/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2_HPLP_input_final_MCC.tree")
class(bs_trees) <- "multiPhylo"

# Plot bootstraps
p <- ggdensitree(bs_trees, layout="slanted", color="black", alpha=0.1)

# Add colours
grp <- list("Guanapo" = c('GHP', 'GLP'),
            "Tacarigua" = c('TACHP', 'TACLP'),
            "Aripo" = c('APHP', 'APLP'),
            "Oropouche" = c('OHP', 'OLP'),
            "Madamas" = c('MADHP', 'MADLP'),
            "P wingei" = "WINGEI")



# Plot Consensus without wingei
snapp_tree_no_wingei <- drop.tip(snapp_tree,"WINGEI")
snapp_tree_no_wingei_grp <- groupOTU(snapp_tree_no_wingei, grp)
consensus_tree <- ggtree(snapp_tree_no_wingei_grp,aes(colour=group),layout="slanted",size=2) +
  theme_tree2(legend.position="none") +
  theme(panel.grid.major.x = element_line(size=1))+
  scale_x_continuous(labels = abs)+
  geom_tiplab(size=8)+
  scale_colour_manual(breaks = c(five_aside_colour_rivers$river,"0"),
                      values = c(five_aside_colour_rivers$colour,"black"))+
 xlim_tree(xlim=c(NA,0.12))+
  geom_range("height_0.95_HPD",center="height", color='firebrick', size=2, alpha=.3)+
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior) < 0.9, 
                 x=branch), hjust=1,colour="black",size=7)

consensus_tree <- consensus_tree+xlab("Time (mya)") + theme(axis.text = element_text(size=20),
                                                            axis.title = element_text(size=22))
consensus_tree <- revts(consensus_tree)
consensus_tree

# Make a table of all of the tree info
ggtree(snapp_tree) + geom_text(aes(label=node), hjust=-.3)

# Get node info
node_info <- data.frame(node_lab=c("Madamas","Oropouche","Aripo","Guanapo","Tacarigua","Northern/Oropouche","Caroni/Northern/Oropouche","All Guppy","Guppy + P.wingei"),
                        node=c(21,19,18,17,14,19,15,13,12),
                        age=NA,
                        age_range=NA)

tree_meta <- data.frame(snapp_tree@data)

# Fill node_info
for(node in node_info$node){
  node_info[node_info$node == node,"age"] <- round(tree_meta[tree_meta$node == node,"height"]*1000000)
  node_info[node_info$node == node,"age_range"] <- paste(as.character(round(tree_meta[tree_meta$node == node,"height_0.95_HPD"][[1]]*1000000)),collapse = "-")
}

# Rename columns
colnames(node_info) <- c("Node Label","Node","Age (Years)","95% HPD Interval (Years)")
write.table(node_info,
            "tables/TableSX_SNAPP_node_ages.txt",
            row.names = F,quote = F,sep = "\t")

# Make even larger plot with introgression and sampling info
map<-readRDS("outputs/five_aside_rivers.rds")


# SFS Plots ---------------------------------------------------------------
rivers <- c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
pops <- list(c("TACHP","TACLP"),
             c("GHP","GLP"),
             c("APHP","APLP"),
             c("OHP","OLP"),
             c("MADHP","MADLP"))
names(pops) <- rivers
sfs_figs <- lapply(rivers,function(river){
  
  # Read in SFS
  sfs <- read.table(paste0("data/",river,"_sfs.txt"))
  
  # Melt
  freq <- melt(as.matrix(sfs), value.name=c("Frq"))
  freq$Var1 <- as.integer(gsub("d1_","",freq$Var1))
  freq$Var2 <- as.integer(gsub("d0_","",freq$Var2))
  
  freq$Frq<- log(freq$Frq)
  freq<-freq[freq$Var2 < nrow(sfs)-freq$Var1,]
  
  
  # Plot
  sfs_plot <- ggplot(freq, aes(Var2, Var1, fill = Frq))+
    geom_tile(size=0.6)+
    scale_fill_gradient2(low = "#00008b", high = "#ffff00", mid = "#4eee94", 
                         midpoint = 5.0, 
                         name="Freq\n(log)") +
    theme_void()+ 
    theme(axis.title.x = element_text(size=24),
          axis.title.y = element_text(size=24,angle=90),
          legend.justification = c(1, 1), 
          legend.position = c(1, 1),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20),
          legend.direction = "vertical")+
    # theme(panel.grid = element_blank(), 
    #       axis.title.x = element_text(face = "bold",vjust = -1.2),
    #       axis.title.y = element_text(face = "bold"), 
    #       legend.position = "none",
    #       axis.text.x = element_text(size = 10, vjust=0.3, margin=margin(-12,0,0,0)),axis.text.y = element_text(size = 10, hjust=0.5, margin=margin(0,-12,0,0)))+
    coord_fixed(ratio=1)  +
    scale_x_continuous(breaks=seq(0,16,1))+
    scale_y_continuous(breaks=seq(0,16,1))+
    labs(x = pops[[river]][2], y= pops[[river]][1])
  
  # Return
  return(sfs_plot)
})

# Plot together
sfs_row <- plot_grid(plotlist = sfs_figs,
                     nrow=1,ncol=5)


# Bring them all together -------------------------------------------------
pdf("figs/five_aside_figure1_final_revised_new_tree.pdf",width=24,height=14)
plot_grid(
  plot_grid(plot_grid(map,consensus_tree,ncol=1,labels=c("A","B"),label_size=30),
            plot_grid(FS_fig,labels="C",label_size=30),
            plot_grid(full_pca_fig,ncol=1,labels=c("D"),label_size=30),
            rel_widths=c(1,1.4,1.3),ncol = 3),
  sfs_row,nrow=2,ncol=1,labels = c("","E"),label_size=30,rel_heights = c(3,1))
dev.off()







