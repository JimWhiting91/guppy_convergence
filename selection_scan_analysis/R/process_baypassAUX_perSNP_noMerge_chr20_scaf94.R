# Script for processing baypass outputs into windowed outliers...

# Load packages
lib<-c("cowplot","ggExtra","ggplot2","data.table","reshape2","parallel","dplyr","tidyr","gggenes")
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

# Save these
saveRDS(dd,
        "outputs/baypass_results_perSNP_noMerge.rds")

##################################################################################################
# Build windowed averages
chrs<-unique(dd$chr)[2:length(unique(dd$chr))]

if(!(file.exists(paste0("outputs/baypass_AUXmodel_windows_",WINDOW,"_noMerge.rds")))){
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
          paste0("outputs/baypass_AUXmodel_windows_",WINDOW,"_noMerge.rds"))
}

winds<-readRDS(paste0("outputs/baypass_AUXmodel_windows_",WINDOW,"_noMerge.rds"))

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

# View the best
winds[winds$outlier == 999,]

# Save
saveRDS(winds,
        paste0("outputs/baypass_AUXmodel_windows_",WINDOW,"_processed_noMerge.rds"))


# Plot chr20 candidate region ---------------------------------------------
outlier_regions <- readRDS("outputs/outlier_candidate_regions_for_all_rivers_10000.rds")
rivers <- names(outlier_regions)
outlier_dd <- data.frame(rbindlist(lapply(rivers,function(x){
  return(data.frame(window_id=outlier_regions[[x]],
                    river=x))
})))

# Separate out
outlier_dd <- outlier_dd %>% separate(window_id,sep=":",into = c("chr","window"))
outlier_dd <- outlier_dd %>% separate(window,sep="-",into = c("start","end"))
outlier_dd$start <- as.integer(outlier_dd$start)
outlier_dd$end <- as.integer(outlier_dd$end)

# Transform for scaf94
outlier_dd[outlier_dd$chr == "chr20","start"] <- merge_scaf94_chr20(outlier_dd[outlier_dd$chr == "chr20","start"],scaf = "chr20")
outlier_dd[outlier_dd$chr == "chr20","end"] <- merge_scaf94_chr20(outlier_dd[outlier_dd$chr == "chr20","end"],scaf = "chr20")
outlier_dd[outlier_dd$chr == "000094F_0","start"] <- merge_scaf94_chr20(outlier_dd[outlier_dd$chr == "000094F_0","start"],scaf = "94")
outlier_dd[outlier_dd$chr == "000094F_0","end"] <- merge_scaf94_chr20(outlier_dd[outlier_dd$chr == "000094F_0","end"],scaf = "94")
outlier_dd[outlier_dd$chr == "000094F_0","chr"] <- "chr20"

# Filter for region
outlier_chr20 <- outlier_dd[outlier_dd$chr == "chr20" & outlier_dd$start < 2633448,]

# Plot SNP only Baypass
snp_dd <- dd[dd$chr %in% c("chr20","000094F_0") & as.integer(dd$BP) < 4000000,]
snp_dd$BP <- as.integer(snp_dd$BP)
snp_dd[snp_dd$chr == "chr20","BP"] <- merge_scaf94_chr20(snp_dd[snp_dd$chr == "chr20","BP"] ,"chr20")
snp_dd[snp_dd$chr == "000094F_0","BP"] <- merge_scaf94_chr20(snp_dd[snp_dd$chr == "000094F_0","BP"] ,"94")
snp_dd <- snp_dd[snp_dd$BP < 2633448,]

# Plot SNP baypass
baypass_snps <- ggplot(snp_dd,aes(x=BP,y=BF.dB.))+
  geom_point(alpha=0.5)+
  labs(y="BF",x="Chromosome 20 Pos (Mb)")+
  theme_bw()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size = 16))+
  scale_x_continuous(breaks=seq(0,2633448,1000000),
                     labels=seq(0,2633448,1000000)/1000000)

# And outlier windows
outlier_chr20$river_F <- factor(outlier_chr20$river,levels = rivers)
outliers <- ggplot(outlier_chr20)+
  geom_segment(aes(x=start,xend=end,y=river_F,yend=river_F,colour=river_F),size=5,show.legend = F)+
  scale_colour_manual(breaks = five_aside_colour_rivers$river_F,
                      values = five_aside_colour_rivers$colour)+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=14),
        panel.grid = element_blank())

# Plot genes in region...
region_genes <- read.csv("tables/TableSX_CL_region_genes_chr20_transformed.csv")
region_genes <- region_genes[region_genes$start < 2133448 & region_genes$end > 1533448 &
                               region_genes$gene != "None",]

# Tidy up
region_genes$gene <- gsub("LOC103482405","plppr4",region_genes$gene)
region_genes$gene <- gsub("LOC103482406","palmdelphin",region_genes$gene)
region_genes$direction <- ifelse(region_genes$strand == "+",1,-1)

genes <- ggplot(region_genes, aes(xmin = start, y=strand,xmax = end, fill = gene,forward=direction,label=gene)) +
        geom_gene_arrow() +
        scale_fill_brewer(palette = "Set3") +
        theme_void()+
        theme(legend.position = "top") +
        labs(fill="")

subset_snps <- baypass_snps +
  xlim(c(1533448,2133448))+
  xlab("Chromosome 20 Pos (BP)")

# Combine
baypass_snp_genes <- plot_grid(genes,subset_snps,
                               ncol=1,nrow=2,rel_heights = c(2,5),
                               axis = "tblr",align="v")
pdf("figs/FigureSX_baypass_per_snp_candidate_genes.pdf",width=10,height=6)
baypass_snp_genes
dev.off()




<