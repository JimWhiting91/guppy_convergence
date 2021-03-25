# Script for processing baypass outputs into windowed outliers...

# Load packages
lib<-c("gggenes","ggExtra","ggplot2","data.table","reshape2","parallel","dplyr","tidyr","cowplot")
lapply(lib,library,character.only=TRUE)
source("~/Exeter/five_aside/five_aside_functions.R")

# Which data are we working with
DATASET="five_aside_STAR"
WINDOW=10000

# Get our SNP locations
snps<-data.frame(fread("data/AF/five_aside_STAR_snps",header = F))
snp_list<-lapply(1:16,function(x){
  SNPs<-data.frame(fread(paste0("data/xtx/",DATASET,"_sub",x,".snps"),header=F))
  SNP_dd <- SNPs %>% separate(V1, c("chr", "BP","BP2"))
  SNP_dd[!(SNP_dd$chr %in% paste0("chr",1:23)), "BP" ]<- SNP_dd[!(SNP_dd$chr %in% paste0("chr",1:23)), "BP2" ]
  SNP_dd<-SNP_dd[,1:2]
  SNP_dd$chr<-gsub("F","F_0",SNP_dd$chr)
  return(SNP_dd)
})

# Read in results
dd<-data.frame(rbindlist(mclapply(1:16,function(x){
  tmp<-data.frame(fread(paste0("data/baypass/five_aside_STAR_auxmodel_sub",x,"_summary_betai.out")))
  
  tmp2<-cbind(tmp,snp_list[[x]])
  return(tmp2)
},mc.cores=4)))

# Liftover chr20
dd$BP <- as.integer(dd$BP)
dd[dd$chr == "chr20","BP"]<-merge_scaf94_chr20(dd[dd$chr == "chr20","BP"],"chr20")
dd[dd$chr == "000094F_0","BP"]<-merge_scaf94_chr20(dd[dd$chr == "000094F_0","BP"],"94")
dd[dd$chr == "000094F_0","chr"]<-"chr20"

# Save these
saveRDS(dd,
        "outputs/baypass_results_perSNP.rds")

##################################################################################################
# Build windowed averages
chrs<-unique(dd$chr)[2:length(unique(dd$chr))]

if(!(file.exists(paste0("outputs/baypass_AUXmodel_windows_",WINDOW,".rds")))){
winds<-data.frame(rbindlist(lapply(chrs,function(x){
  tmp<-dd[dd$chr == x,]
  tmp<-tmp[order(as.integer(tmp$BP)),]
  tmp$BP<-as.integer(tmp$BP)
  tmp_wind<-make_windows(tmp,"BF.dB.",0,max(tmp$BP),WINDOW,avg="mean")
  tmp_wind2<-make_windows(tmp,"BF.dB.",0,max(tmp$BP),WINDOW,avg="count")
  tmp_wind$chr<-x
  tmp_wind$count<-tmp_wind2$BF.dB.
return(na.omit(tmp_wind))
})))

colnames(winds)[3]<-"BayP"

# Save the baypass windows to an output
winds$window_id<-paste0(winds$chr,":",as.integer(winds$BP1),"-",as.integer(winds$BP2))
saveRDS(winds,
        paste0("outputs/baypass_AUXmodel_windows_",WINDOW,".rds"))
}

winds<-readRDS(paste0("outputs/baypass_AUXmodel_windows_",WINDOW,".rds"))

##################################################################################################
# Filter windows with low support
winds<-winds[winds$count > 5 & winds$count < 500,]

# Fetch outliers 
upper95<-quantile(winds$BayP,0.95)
upper99<-quantile(winds$BayP,0.99)
upper999<-quantile(winds$BayP,0.999)

# Set them
winds$outlier<-"No"
winds[winds$BayP > upper95,"outlier"]<-"95"
winds[winds$BayP > upper99,"outlier"]<-"99"
winds[winds$BayP > upper999,"outlier"]<-"999"

# Check each set of outliers
best<-winds[winds$outlier=="999",]
best[order(-best$BayP),]
better<-winds[winds$outlier=="99",]
good<-winds[winds$outlier=="95",]


# And plot to check
ggplot(winds[winds$chr == "chr20",],aes(x=rowMeans(winds[winds$chr == "chr20",c("BP1","BP2")]),y=BayP))+
  geom_line()

# Plot the whole genome?
winds$chrN<-gsub("chr","",winds$chr)
winds$chr_F<-factor(winds$chrN,levels = 1:23)
winds$chrN<-gsub("chr","",winds$chr_F)
winds<-na.omit(winds)

# Fix this
winds$window_id<-paste0(winds$chr,":",as.integer(winds$BP1),"-",as.integer(winds$BP2))

# Plot whole genome
baypass_genome<-ggplot(na.omit(winds),aes(x=rowMeans(winds[,c("BP1","BP2")]),y=BayP))+
        geom_point(alpha=0.5,size=2,aes(colour=outlier),show.legend = F)+
        facet_wrap(~chr_F,ncol=23,scales = "free_x",strip.position = "bottom")+
        theme_minimal()+
        theme(axis.text.x = element_blank(),
              panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        axis.title = element_text(size=20),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=18))+
  scale_colour_manual(breaks = c("No","95","99","999"),
                      values = c("gray50","gold2","orange2","red4"))+
  xlab("Genome Position (Chr)")+
  ylab("HP/LP Association (BF)")+
 # geom_hline(yintercept = upper95,colour="red2",linetype="dotted")+
 # geom_hline(yintercept = upper99,colour="red2",linetype="dashed")+
  geom_hline(yintercept = upper999,colour="red2",linetype="solid")+
  geom_hline(yintercept = median(winds$BayP),colour="black",linetype="dashed")

# and plot specific chrs
baypass_chr_plot<-function(x){
  ggplot(na.omit(winds[winds$chr== paste0("chr",x),]),aes(x=rowMeans(winds[winds$chr==paste0("chr",x),c("BP1","BP2")]),y=BayP))+
    geom_point(alpha=0.5,size=2,aes(colour=outlier),show.legend = F)+
    #facet_wrap(~chr_F,ncol=23,scales = "free_x")+
    theme_minimal()+
    theme(panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.title = element_text(size=20),
          strip.text = element_text(size=16),
          axis.text.y = element_text(size=18),
          axis.text.x = element_text(size=14))+
    scale_colour_manual(breaks = c("No","95","99","999"),
                        values = c("gray50","gold2","orange2","red4"))+
  xlab(paste0("Chr ",x," Position (Mb)"))+
  ylab("HP/LP Association (BF)")+
 # geom_hline(yintercept = upper95,colour="red2",linetype="dotted")+
  geom_hline(yintercept = upper999,colour="red2",linetype="solid")+
  geom_hline(yintercept = median(winds$BayP),colour="black",linetype="dashed")+
  scale_x_continuous(breaks = seq(0,max(winds[winds$chr==paste0("chr",x),"BP2"]),2000000),
                     labels = seq(0,max(winds[winds$chr==paste0("chr",x),"BP2"]),2000000)/1000000)
}

# Plot them
baypass_chr_plot(20)
baypass_chr_plot(8)
baypass_chr_plot(4)
baypass_chr_plot(1)
baypass_chr_plot(9)
baypass_chr_plot(13)

# Plot Just chr8 and 20
winds$chr_F2<-factor(winds$chr,levels = paste0("chr",1:23))
baypass_8_and_20<-ggplot(na.omit(winds[winds$chr %in% c("chr8","chr20"),]),aes(x=rowMeans(winds[winds$chr %in% c("chr8","chr20"),c("BP1","BP2")]),y=BayP))+
  geom_point(alpha=0.5,size=2,aes(colour=outlier),show.legend = F)+
  facet_wrap(~chr_F2,ncol=2,scales = "free_x",strip.position = "right")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        axis.title = element_text(size=20),
        strip.text = element_text(size=16),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=14))+
  scale_colour_manual(breaks = c("No","95","99","999"),
                      values = c("gray50","gold2","orange2","red4"))+
  xlab("Chromosome Position (Mb)")+
  ylab("HP/LP Association (BF)")+
  # geom_hline(yintercept = upper95,colour="red2",linetype="dotted")+
  # geom_hline(yintercept = upper99,colour="red2",linetype="dashed")+
  geom_hline(yintercept = upper999,colour="red2",linetype="solid")+
  geom_hline(yintercept = median(winds$BayP),colour="black",linetype="dashed")+
  scale_x_continuous(breaks = seq(0,max(winds[winds$chr %in% c("chr8","chr20"),"BP2"]),2000000),
                     labels = seq(0,max(winds[winds$chr %in% c("chr8","chr20"),"BP2"]),2000000)/1000000)


# View the best
winds[winds$outlier == 999,]

# Save
saveRDS(winds,
        paste0("outputs/baypass_AUXmodel_windows_",WINDOW,"_processed_Scaf94_merged.rds"))

# Also save figures
saveRDS(list(baypass_genome,
             baypass_8_and_20),
        "outputs/baypass_figures_genome_chr20,8.rds")


# Plot SNP baypass results ------------------------------------------------
# First merge scaf94 and chr20
dd2<-dd[2:nrow(dd),]
  
snp_baypass_fig <- function(data,chr,start,end){
  # Filter data
  tmp <- data[data$chr == chr &
                data$BP < end &
                data$BP > start,]
  
  # Plot
  ggplot(tmp,aes(y=BF.dB.,x=BP))+
    geom_point()
}

snp_baypass_fig(dd2,"chr20",0,4000000)
snp_baypass_fig(dd2,"chr15",0,6000000)
snp_baypass_fig(dd2,"chr1",11000000,13000000)
snp_baypass_fig(dd2,"chr1",12200000,12300000)
snp_baypass_fig(dd2,"chr8",17000000,18000000)

# Plot CL-AP region SNPs with coverage diffs...
CL_baypass <- snp_baypass_fig(dd2,"chr20",836423+1797025-1100000,836423+1797025-500000)
CL_baypass <- CL_baypass + 
  theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=10))+
  xlab("Chromosome 20 Position (Mb)")+
  ylab("HP/LP Association (BF)")+
  scale_x_continuous(breaks=seq(1550000,2100000,100000),
                     labels=seq(1550000,2100000,100000)/1000000)+
  geom_hline(yintercept = upper999,colour="red2",linetype="solid")

  
  
# Fetch coverage data
coverage <- data.frame(fread("~/Exeter/five_aside/five_aside_coverage/data/five_aside_STAR_coverage_ratios_1000.txt"))
coverage_CL <- coverage[coverage$chr == "000094F_0" &
                          coverage$BP1 >= 500000 &
                          coverage$BP2 <= 1100000,]
coverage_CL$BP1_STAR <- merge_scaf94_chr20(coverage_CL$BP1,scaf = "94")
coverage_CL$BP2_STAR <- merge_scaf94_chr20(coverage_CL$BP2,scaf = "94")


# Make coverage tracks...
coverage_CL$river_F <- factor(coverage_CL$river,levels=five_aside_colour_rivers$river)
coverage_fig <- ggplot(coverage_CL,aes(x=rowMeans(coverage_CL[,c("BP1_STAR","BP2_STAR")]),y=log2(cov_ratio),colour=river_F))+
  geom_line(show.legend = F,size=1)+
  facet_wrap(~river_F,ncol=1,strip.position = "right",scales = "free_y")+
  scale_colour_manual(breaks=five_aside_colour_rivers$river_F,
                      values=five_aside_colour_rivers$colour)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size=10))+
  labs(y=expression(HP/LP~Coverage~(log[2])))+
  scale_x_continuous(breaks=seq(1550000,2100000,100000),
                     labels=seq(1550000,2100000,100000)/1000000)

# Stack
plot_grid(coverage_fig,CL_baypass,
          ncol=1,axis = "lr",align="v",rel_heights = c(2,1))

# Make geneplot...
CL_genes <- read.csv("data/scaf94_annotated_genes_BUCK.csv",header=F)
CL_genes <- CL_genes[CL_genes$V2 < 1100000 &
                       CL_genes$V3 > 500000,]
CL_genes$V3 <- merge_scaf94_chr20(CL_genes$V3,scaf="94")
CL_genes$V2 <- merge_scaf94_chr20(CL_genes$V2,scaf="94")
CL_genes$gene <- c("palmdelphin","plppr4","plppr5","Novel","Novel","Novel","Novel","snx7","LOC103482394","LOC108165682")

# Make gene plot...
clap_genes <- ggplot(CL_genes,aes(xmin=V3,xmax=V2,y=gene,fill=gene))+
  geom_gene_arrow(arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_x_continuous(breaks=seq(1550000,2100000,100000),
                     labels=seq(1550000,2100000,100000)/1000000)+
  xlim(1533448,2133448)

# Stack
full_clap_region_fig <- plot_grid(clap_genes,coverage_fig,CL_baypass,
              ncol=1,axis = "lr",align="v",rel_heights = c(0.5,2,0.75))

# Repeat the plot but limit exclusively the plppr5 regions.
plppr5_exons <- read.csv("data/plppr5_exons.csv",header=F)
plppr5_exons$V2 <- merge_scaf94_chr20(plppr5_exons$V2,scaf = "94")
plppr5_exons$V3 <- merge_scaf94_chr20(plppr5_exons$V3,scaf = "94")
plppr5_exons$V4 <- "plppr5 exon"
plppr5_exons_fig <- ggplot(plppr5_exons,aes(xmin=V3,xmax=V2,y=V4,fill=V4))+
  geom_gene_arrow(arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank())+
  xlim(1899738,2023476)

plppr5_region_fig <- plot_grid(plppr5_exons_fig,
                                  coverage_fig+xlim(1899738,2023476),
                                  CL_baypass+xlim(1899738,2023476),
                                  ncol=1,axis = "lr",align="v",rel_heights = c(0.5,2,0.75))

# Plot together...
pdf("figs/FigureSX_CLAP_region_BF_coverage_genes.pdf",height=10,width=8)
plot_grid(full_clap_region_fig,
          plppr5_region_fig,
          ncol=2,labels="AUTO",label_size=20)
dev.off()

# plot_grid(clap_genes+xlim(1859738,2023476),
#           coverage_fig+xlim(1859738,2023476),
#           CL_baypass+xlim(1859738,2023476),
#           ncol=1,axis = "lr",align="v",rel_heights = c(0.5,2,1))


##################################################################################################
# Frequency of highest Baypass SNP - 000094F_0:556282


# ##################################################################################################
# # Set the cut-off
# upperQ<-quantile(dd$BF.dB.,0.999)
# 
# # Build windows
# window_size<-WINDOW
# chrs<-unique(dd$chr)
# chrs<-chrs[chrs!="x"]
# 
# # Set
# dd$BP<-as.integer(dd$BP)
# 
# baypass_winds<-data.frame(rbindlist(lapply(chrs,function(x){
#   
#   tmp<-dd[dd$chr == x,]
#   
#   winds<-seq(0,max(tmp$BP),window_size)
#   winds2<-winds+window_size
#   
#   # Go through windows
#   window_out<-data.frame(rbindlist(lapply(1:length(winds),function(y){
#     
#     tmp2<-tmp[tmp$BP < winds2[y] &
#                 tmp$BP > winds[y],]
#     
#     if(nrow(tmp2) > 0){
#       SNP_N<-nrow(tmp2)
#       outlier_N<-nrow(tmp2[tmp2$BF.dB. >= upperQ,])
#     } else {
#       SNP_N<-NA
#       outlier_N<-NA
#     }
#     
#     out<-data.frame(chr=x,
#                     start=winds[y],
#                     end=winds2[y],
#                     SNP_N=SNP_N,
#                     outlier_N=outlier_N)
#     
#     return(out)
#   })))
#   
#   return(window_out)
#   
# })))
# 
# # Filter duds
# baypass_winds<-na.omit(baypass_winds)
# 
# ##### FIND OUTLIERS #######
# # Calculate cutoff
# SNP_sim<-as.data.frame(seq(1,max(baypass_winds$SNP_N),by=1))
# colnames(SNP_sim)<-"SNP_N"
# 
# # We need to know what the probability for expected number of SNPs is which = N of outlier SNPs over total SNPs
# p<-sum(baypass_winds$outlier_N)/sum(baypass_winds$SNP_N)
# 
# # Here we calculate binomial expectation for 0.99 quantile
# for (k in 1:max(baypass_winds$SNP_N)){
#   SNP_sim$exp[k]<-qbinom (0.999, k, p)
# }
# 
# #Calculate outliers
# baypass_winds$exp<-qbinom(0.999, baypass_winds$SNP_N, p)
# 
# # Write intermediate output to file
# top_outliers<-baypass_winds[baypass_winds$outlier_N > baypass_winds$exp,]
# 
# # Reorder my strength of effect
# top_outliers$residual_N<-top_outliers$outlier_N-top_outliers$exp
# top_outliers<-top_outliers[order(-top_outliers$residual_N),]
# 
# write.table(top_outliers,
#             paste0("outputs/",DATASET,"_auxmodel_baypass_BF_windows_",window_size,"_OUTLIERS.txt"),
#             row.names=F,quote=F,sep="\t")
