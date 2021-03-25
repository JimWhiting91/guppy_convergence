###########################
# Figures for scaffold 94/chr20 region
lib<-c("ggplot2","data.table","vcfR","cowplot","ggtree","GenotypePlot")
lapply(lib,library,character.only=T)
source("~/software/genotype_plot.R")
source("~/Exeter/five_aside/five_aside_functions.R")
devtools::load_all("~/Exeter/genotype_plot/")

###########################
# A) Genotypes
vcf="~/Exeter/VCFs/five_aside_STAR_chr20_scaf94_merged_shapeit_beagle.vcf.gz"

popmap=data.frame(ind=unlist(list(read.table("~/Exeter/VCFs/LT.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/UT.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/GH.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/GL.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/APHP.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/APLP.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/LO.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/UQ.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/LMD.popmap",stringsAsFactors = F)[,1],
            read.table("~/Exeter/VCFs/UMD.popmap",stringsAsFactors = F)[,1])))

pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
pop_labs <-c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
names(pop_labs) <- pops
for(pop in pops){
  popmap[grep(pop,popmap$ind),"pop"] <- pop_labs[pop]
}
CL_genotype_plot<-genotype_plot(vcf=vcf,
                             chr="chr20",
                             start=0,
                             end=2633448,
                             popmap=popmap,
                             cluster=TRUE,
                             snp_label_size=500000)
                           #  colour_scheme = c("#FFCC99","#FF8000","#660000"))

# Assemble the genotype plot
# Dendrogram
dendrogram_metadata <- data.frame(ind = CL_genotype_plot$dendro_labels)
dendrogram_metadata$tip <- 1:nrow(dendrogram_metadata)
rivers <- c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
rivers2 <- rep(rivers,each=2)
pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
pred <- rep(c("HP","LP"),5)

# Build metadata
dendrogram_metadata$river <- NA
dendrogram_metadata$pred <- NA
for(i in 1:length(pops)){
  dendrogram_metadata[grep(pops[i],dendrogram_metadata$ind),"river"]<-rivers2[i]
  dendrogram_metadata[grep(pops[i],dendrogram_metadata$ind),"pred"]<-pred[i]
}

# Add tip labels
CL_genotype_plot$dendrogram <- CL_genotype_plot$dendrogram + 
  geom_jitter(aes(y=-2.5,x=dendrogram_metadata$tip,colour=dendrogram_metadata$river,shape=dendrogram_metadata$pred),
              size = 2.2,height = 2.2,width=0,alpha=0.75)+
  scale_colour_manual(breaks=five_aside_colour_rivers$river,
                      values=five_aside_colour_rivers$colour)+
  scale_shape_manual(breaks=c("HP","LP"),
                     values=c(19,17))+
  theme(legend.position="left",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  labs(colour="River",shape="Predation")+
  guides(colour = guide_legend(override.aes = list(size=8)),
         shape = guide_legend(override.aes = list(size=8)))



# # Find the SNPs that fall within the CL-AP region
# SNPs between (836423 + 647025) & 836423+1297025)
vcf_in<-read.vcfR(vcf)
vcf_in<-vcf_in[as.integer(vcf_in@fix[,2]) < 2633448,]
bp<-as.integer(vcf_in@fix[,2])
plotting_pos<-seq(min(bp),max(bp),by=(max(bp)-min(bp))/length(bp))[1:length(bp)]
CLAP_start<-plotting_pos[min(which(as.integer(vcf_in@fix[,2]) > 1483448))]
CLAP_end<-plotting_pos[max(which(as.integer(vcf_in@fix[,2]) < 2133448))]

# Merge with genotype plot and annotate
CL_genotype_plot$genotypes <- CL_genotype_plot$genotypes+
  geom_vline(xintercept=c(CLAP_start,CLAP_end),size=1.5,colour="black")

# Make full
CL_genotype_plot$positions <- CL_genotype_plot$positions +
  theme(axis.text = element_text(size=20))
CL_combined <- combine_genotype_plot(CL_genotype_plot,heights = c(1,12))

###########################
# B) LocalPCAs
pops<-c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
new_pops<-c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
names(new_pops)<-pops

local_pca_res<-data.frame(rbindlist(lapply(pops,function(pop){
  tmp<-read.table(paste0("data/scaf94_figs/five_aside_STAR_",pop,"_chr20_scaf94_merge.txt"),header=T)
  tmp$pop<-new_pops[pop]
  return(tmp)
})))

# Set river
local_pca_res$river<-"Aripo"
local_pca_res[grep("G",local_pca_res$pop),"river"]<-"Guanapo"
local_pca_res[grep("TAC",local_pca_res$pop),"river"]<-"Tacarigua"
local_pca_res[grep("O",local_pca_res$pop),"river"]<-"Oropouche"
local_pca_res[grep("MAD",local_pca_res$pop),"river"]<-"Madamas"
local_pca_res$river_F<-factor(local_pca_res$river,levels=c("Tacarigua",
                                                           "Guanapo",
                                                           "Aripo",
                                                           "Oropouche",
                                                           "Madamas"))
# Set up correct order
local_pca_res$pop_F<-factor(local_pca_res$pop,levels=new_pops)

# Add High LP
local_pca_res$pred <- "HP"
local_pca_res[grep("LP",local_pca_res$pop),"pred"] <- "LP"

# Subset for HP
#local_pca_res<-local_pca_res[local_pca_res$pop %in% c("APHP","TACHP","GHP","OHP","MADHP"),]

# Separate start and end for windows
windows <- sapply(strsplit(local_pca_res$window_id,":"),'[[',2)
local_pca_res$window_start <- as.numeric(sapply(strsplit(windows,"-"),'[[',1))

mds_plot <-ggplot(local_pca_res,aes(x=window_start,y=mds1,colour=river_F,shape=pred))+
  annotate("rect", xmin = 0, xmax = 836423+1797025, ymin = -Inf, ymax = Inf,
           alpha = .5,fill = "grey75")+
  geom_step()+
  facet_grid(river_F~pred,scales="free_y")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size=15),
        legend.position = "none",
       strip.background = element_rect(fill="white"),
        title = element_text(size = 22),
        axis.title=element_text(size=20))+
  scale_x_continuous(breaks = seq(0,max(local_pca_res$window_pos),5000000),
                     labels = seq(0,max(local_pca_res$window_pos),5000000)/1000000)+
  labs(y="MDS1",x="Chromosome 20 Pos (Mb)")+
  geom_hline(yintercept=0,linetype="dotted")+
  scale_colour_manual(breaks=five_aside_colour_rivers$river_F,
                      values = five_aside_colour_rivers$colour)+
  scale_shape_manual(breaks=c("HP","LP"),
                     values=c(19,17))+
  geom_vline(xintercept = c((836423 + 647025),(836423+1297025)))

# Test
# pdf("figs/test.pdf",width=16,height=8)
# plot_grid(genos_with_lines,
#           mds_plot,ncol=2,rel_widths = c(7,3))
# dev.off()

###########################
# C) Tree
scaf94_tree<-read.tree("data/scaf94_figs/20_07_01_scaf94_CLAP_region_haplogroups_TFINAL.raxml.support")

# Inspect
ggtree(scaf94_tree)+geom_text(aes(label=node), hjust=-.3,size=2)
ggtree(scaf94_tree)+geom_tiplab()

# ggtree(scaf94_tree,layout = "equal_angle")

# Get the tip points
tree_metadata <- data.frame(taxa=scaf94_tree$tip.label)
tree_metadata$river<-"Aripo"
tree_metadata$pred<-"HP"

# Set rivers
tree_metadata[grep("GH",tree_metadata$taxa),"river"]<-"Guanapo"
tree_metadata[grep("GL",tree_metadata$taxa),"river"]<-"Guanapo"
tree_metadata[grep("UT",tree_metadata$taxa),"river"]<-"Tacarigua"
tree_metadata[grep("LT",tree_metadata$taxa),"river"]<-"Tacarigua"
tree_metadata[grep("LO",tree_metadata$taxa),"river"]<-"Oropouche"
tree_metadata[grep("UQ",tree_metadata$taxa),"river"]<-"Oropouche"
tree_metadata[grep("UMD",tree_metadata$taxa),"river"]<-"Madamas"
tree_metadata[grep("LMD",tree_metadata$taxa),"river"]<-"Madamas"

# Set predation regime
# Set rivers
tree_metadata[grep("APLP",tree_metadata$taxa),"pred"]<-"LP"
tree_metadata[grep("GL",tree_metadata$taxa),"pred"]<-"LP"
tree_metadata[grep("UT",tree_metadata$taxa),"pred"]<-"LP"
tree_metadata[grep("UQ",tree_metadata$taxa),"pred"]<-"LP"
tree_metadata[grep("UMD",tree_metadata$taxa),"pred"]<-"LP"

# And final label
tree_metadata$pop<-paste0(tree_metadata$river,"_",tree_metadata$pred)

# Set factor
tree_metadata$river_F<-factor(tree_metadata$river,levels=c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"))

# Plot tree
p <- ggtree(scaf94_tree,layout = "daylight",colour="gray60")
p2 <- p %<+% tree_metadata + 
  geom_tippoint(aes(color=river,shape=pred),size=4,alpha=0.75)+
  scale_colour_manual(breaks = five_aside_colour_rivers$river_F,
                      values = five_aside_colour_rivers$colour)+
  scale_shape_manual(breaks = c("HP","LP"),
                     values = c(19,17))+
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = "none",
        legend.direction = "vertical")+
  labs(shape="Predation",colour="River")+
  coord_flip()
plot(p2)

#### Final plot together with branch length analysis #####
scaf94_branch_fig <- readRDS("outputs/full_branch_fig.rds")
# pdf("figs/Figure3_scaf94_analyses.pdf",width=20,height=13)
# plot_grid(
#   plot_grid(CL_combined,
#           plot_grid(mds_plot,p2,labels=c("B","C"),label_size = 24,ncol=1,nrow=2,rel_heights = c(1.2,1),axis = "tblr"),
#           labels = c("A",""),label_size = 24,rel_heights = c(1,1),ncol = 2,nrow = 1,axis = "tblr",align="h"),
#   scaf94_branch_fig,
#   ncol = 1,nrow = 2,rel_heights = c(2,1))
# 
# dev.off()


pdf("figs/Figure3_scaf94_analyses.pdf",width=22,height=17)
plot_grid(
  CL_combined,
  plot_grid(mds_plot,
  plot_grid(p2,scaf94_branch_fig,labels = c("C",""),rel_widths = c(1,2),label_size = 24),
  ncol=1,nrow=2,
  labels=c("B",""),label_size=24),
  ncol = 2,nrow = 1,rel_widths = c(1.5,2),
  labels = c("A",""),label_size=24,align = "h",axis = "v")

dev.off()

