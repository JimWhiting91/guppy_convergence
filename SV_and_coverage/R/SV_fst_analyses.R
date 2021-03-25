######### Quick analysis of SV FSTs
lib<-c("cowplot","vcfR","hierfstat","data.table","ggplot2","qqman","adegenet","gggenes")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

rivers<-c("Tacarigua",
          "Guanapo",
          "Aripo",
          "Oropouche",
          "Madamas")

ids<-list(c("LT","UT"),
          c("GH","GL"),
          c("APHP","APLP"),
          c("LO","UQ"),
          c("LMD","UMD"))

names(ids)<-rivers

#  Per river Fst ----------------------------------------------------------

fst_results<-lapply(rivers,function(river){
  
  # Read in the popmap
  popmap<-read.table(paste0("~/Exeter/VCFs/",river,".popmap"))[,1]
  
  # Divide into HP LP
  HP <- grep(ids[[river]][1],popmap,value = T)
  LP <- grep(ids[[river]][2],popmap,value = T)
  popmap2<-c(rep("HP",length(HP)),
             rep("LP",length(LP)))
  
  # Read in SVs
  sv_vcf<-read.vcfR(paste0("~/Exeter/VCFs/",river,"_SV_filtered_smaller_indels-smoove.genotyped.vcf.gz"))
  
  # Transform...
  dat<-vcfR2genind(sv_vcf)
  pop(dat)<-popmap2
  dat2<-genind2hierfstat(dat)
  
  # Calculate
  stats<-basic.stats(dat2,diploid = 2,digits = 2)
  
  # Make output
  fst_res<-na.omit(data.frame(chr = sv_vcf@fix[,1],
                      pos = as.integer(sv_vcf@fix[,2]),
                      fst = stats$perloc$Fst,
                      river = river,
                      snp = sv_vcf@fix[,3]))

  fst_res2<-fst_res[grep("chr",fst_res$chr),]
  # Set scaffolds to "chr24"
  fst_res2
  fst_res2$chr2<-as.integer(gsub("chr","",fst_res2$chr))
  
  # And plot
  manhat<-manhattan(fst_res2, chr="chr2", bp="pos", snp="snp", p="fst",logp = FALSE)
  
  # Return
  return(list(fst_res,manhat))
})

# View each set of top Fst SVs...
head(fst_results[[1]][[1]][order(-fst_results[[1]][[1]]$fst),],20)
head(fst_results[[2]][[1]][order(-fst_results[[2]][[1]]$fst),],20)
head(fst_results[[3]][[1]][order(-fst_results[[3]][[1]]$fst),],20)
head(fst_results[[4]][[1]][order(-fst_results[[4]][[1]]$fst),],20)
head(fst_results[[5]][[1]][order(-fst_results[[5]][[1]]$fst),],20)

####################################################################
# Plot Oropouche Fst/SV Fst around B-cadherin on chr15...
xtx <- data.frame(rbindlist(readRDS("~/Exeter/five_aside/selection_scanning/outputs/processed_baypass_xtx_all_pops_window10000.rds")))
xtx <- xtx[xtx$river == "Oropouche" & xtx$chr == "chr15",]
xtx$mid <- rowMeans(xtx[,c("BP1","BP2")])
fst_plot <- fst_results[[4]][[1]]
fst_plot <- fst_plot[fst_plot$chr == "chr15",]

# Plot them both
svfst <- ggplot(fst_plot,aes(x=pos,y=fst))+
  geom_point()+
  theme_minimal()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y=expression(SV~F[ST]))+
  scale_x_continuous(breaks=seq(0,max(xtx$mid),5e6),
                     labels=seq(0,max(xtx$mid),5e6)/1e6)

chr15_xtx <- ggplot(xtx,aes(x=mid,y=M_XtX))+
  geom_point()+
  theme_minimal()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  labs(y="SNP XTX",x="Chromosome 15 Position (Mb)")+
  scale_x_continuous(breaks=seq(0,max(xtx$mid),5e6),
                     labels=seq(0,max(xtx$mid),5e6)/1e6)

# Also add B-cadherin at chr15:5030641-5056311
Bcad_pos <- data.frame(molecule="B-Cadherin",
                       gene="B-Cadherin",
                       start=5030641,
                       end=5056311,
                       strand="forward",
                       direction=1)
Bcadherin <- ggplot(Bcad_pos, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left") +
  scale_fill_brewer(palette = "Set3") +
  theme_void()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title=element_blank())+
  xlim(5e6,5.1e6)

Bcad_fst <- svfst + xlim(5e6,5.1e6)
Bcad_xtx <- chr15_xtx + xlim(5e6,5.1e6) + theme(axis.text.x = element_text(size=12,angle=30,hjust=1)) + xlab("Chromosome 15 Position (BP)")

# Merge them all
pdf("figs/FigureSX_Oropouche_SV_FST_chr15.pdf",width=16,height=5)
plot_grid(plot_grid(svfst,chr15_xtx,ncol=1,axis = "tblr",align="v",rel_heights = c(4,6)),
          plot_grid(Bcadherin,Bcad_fst,Bcad_xtx,ncol=1,axis = "tblr",align="v",rel_heights = c(1,4,6)),
          ncol=2,labels="AUTO",rel_widths = c(3,2),label_size = 20)
dev.off()
