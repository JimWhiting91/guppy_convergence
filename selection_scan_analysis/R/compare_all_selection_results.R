#####################################################################################
# Combine all selection scan results and combine results to compare

# Load packages
lib<-c("Rfast","vcfR","regioneR","ggplotify","ggplot2","data.table","reshape2","parallel","dplyr","tidyr","UpSetR","cowplot")
lapply(lib,library,character.only=TRUE)
source("~/Exeter/five_aside/five_aside_functions.R")

# Window size
wind_size<-10000

# Read in selection results
xtx<-data.frame(rbindlist(readRDS(paste0("outputs/processed_baypass_xtx_all_pops_window",wind_size,".rds"))))
AFD<-readRDS(paste0("outputs/AFD_calculations_five_aside_STAR_",wind_size,".rds"))[[2]]
xpehh<-readRDS(paste0("outputs/xpehh_windows_window_",wind_size,".rds"))

# Which rivers
rivers<-unique(AFD$river)

#####################################################################################
# Set outliers
# AFD
AFD$window_id<-paste0(AFD$chr,":",as.integer(AFD$BP1),"-",as.integer(AFD$BP2))
AFD<-data.frame(rbindlist(lapply(rivers,function(x){
  tmp<-AFD[AFD$river == x,]
  
  # Filter out windows with few SNPs
  tmp<-tmp[tmp$SNP_count > 5,]
  
  # Define outliers
  tmp$outlier<-"No"
  # cutoff<-tmp[tmp$AFD > quantile(tmp$AFD,0.95),"outlier"]<-95
  # tmp[tmp$AFD > quantile(tmp$AFD,0.99),"outlier"]<-99
  # tmp[tmp$AFD > quantile(tmp$AFD,0.999),"outlier"]<-999
  
  cutoff<-quantile(tmp$AFD,0.95)
  if(cutoff > 0.5){
    cutoff<-0.5
  }
  
  # Set outliers
  tmp[tmp$AFD > cutoff,"outlier"]<-95
  
  # Return
  return(tmp)
})))

# XP-EHH
xpehh<-na.omit(xpehh)
xpehh$outlier<-"No"
xpehh$normxpehh<-as.numeric(xpehh$normxpehh)
xpehh[abs(xpehh$normxpehh) > 2,"outlier"]<-"95"
xpehh[abs(xpehh$normxpehh) > 3,"outlier"]<-"99"

# XTX
xtx$window_id<-paste0(xtx$chr,":",as.integer(xtx$BP1),"-",as.integer(xtx$BP2))

#####################################################################################
# For each measure, plot genome-wide distribution of outliers, fixing the scaffold 94 misplacement.
AFD$mid<-rowMeans(AFD[,c("BP1","BP2")])
AFD[AFD$chr == "chr20","mid"]<-merge_scaf94_chr20(AFD[AFD$chr == "chr20","mid"],"chr20")
AFD[AFD$chr == "000094F_0","mid"]<-merge_scaf94_chr20(AFD[AFD$chr == "000094F_0","mid"],"94")
AFD[AFD$chr == "000094F_0","chr"]<-"chr20"

xtx$mid<-rowMeans(xtx[,c("BP1","BP2")])
xtx[xtx$chr == "chr20","mid"]<-merge_scaf94_chr20(xtx[xtx$chr == "chr20","mid"],"chr20")
xtx[xtx$chr == "000094F_0","mid"]<-merge_scaf94_chr20(xtx[xtx$chr == "000094F_0","mid"],"94")
xtx[xtx$chr == "000094F_0","chr"]<-"chr20"

xpehh$mid<-rowMeans(xpehh[,c("BP1","BP2")])
xpehh[xpehh$chr == "chr20","mid"]<-merge_scaf94_chr20(xpehh[xpehh$chr == "chr20","mid"],"chr20")
xpehh[xpehh$chr == "000094F_0","mid"]<-merge_scaf94_chr20(xpehh[xpehh$chr == "000094F_0","mid"],"94")
xpehh[xpehh$chr == "000094F_0","chr"]<-"chr20"
xpehh$abs_xpehh<-abs(xpehh$normxpehh)

plot_measure_genome_wide<-function(x,measure){
  tmp<-x
  tmp$chrN<-gsub("chr","",tmp$chr)
  tmp$chr_F<-factor(tmp$chrN,levels = 1:23)
  tmp$river_F<-factor(tmp$river,levels=rivers)
  
  # Fetch the cutoffs
  cutoffs<-data.frame(river=rivers,
                      q95=NA,
                      q99=NA,
                      q999=NA)
  
  if(measure != "abs_xpehh"){
  for(i in rivers){
    cutoffs[cutoffs$river == i,"q95"]<-quantile(tmp[tmp$river == i,measure],0.95)
    cutoffs[cutoffs$river == i,"q99"]<-quantile(tmp[tmp$river == i,measure],0.99)
    cutoffs[cutoffs$river == i,"q999"]<-quantile(tmp[tmp$river == i,measure],0.999)
  }
  } else {
    cutoffs$q95<-2
    cutoffs$q99<-3
    cutoffs$q999<-0
  }
  
  cutoffs$river_F<-factor(cutoffs$river,levels=rivers)
  
  tmp<-na.omit(tmp)
  
 p1<- ggplot(tmp,aes(x=mid,y=tmp[,measure],colour=river_F))+
    geom_point(alpha=0.5,size=0.75)+
    facet_grid(river_F~chr_F,scales = "free")+
   theme_minimal()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          panel.spacing.y=unit(0.1, "lines"),
          axis.title = element_text(size=20),
          strip.text = element_text(size=16),
          axis.text.y = element_text(size=12),
          legend.position="none")+
   scale_colour_manual(values = five_aside_colour_rivers$colour,
                       breaks = five_aside_colour_rivers$river_F)+
    xlab("Genome Position (Chr)")+
    ylab(measure)
   #  geom_hline(data = cutoffs, aes(yintercept = q95),colour="red2",linetype="dotted")+
   # geom_hline(data = cutoffs,aes( yintercept = q99),colour="red2",linetype="dashed")+
   # geom_hline(data = cutoffs, aes(yintercept = q999),colour="red2",linetype="solid")
return(p1)  
}

# Plot them
# AFD
afd_cutoffs<-data.frame(rbindlist(lapply(rivers,function(x){
  tmp<-AFD[AFD$river == x,]
  cutoff<-quantile(tmp$AFD,0.95)
  if(cutoff > 0.5){
    cutoff<-0.5
  }
  data.frame(river=x,
             cutoff=cutoff)
})))
afd_cutoffs$river_F<-factor(afd_cutoffs$river,levels = afd_cutoffs$river)
afd_plot<-plot_measure_genome_wide(AFD,"AFD")+
  geom_hline(data=afd_cutoffs,aes(yintercept = cutoff))

# XPEHH
xpehh_plot<-plot_measure_genome_wide(xpehh,"abs_xpehh")+ylab("XP-EHH")+
  geom_hline(yintercept = 2)

# XtX
xtx_cutoffs<-data.frame(rbindlist(lapply(rivers,function(x){
  tmp<-xtx[xtx$river == x,]
  cutoff<-min(tmp[tmp$outlier != "No","M_XtX"])
  data.frame(river=x,
             cutoff=cutoff)
})))
xtx_cutoffs$river_F<-factor(xtx_cutoffs$river,levels = xtx_cutoffs$river)
xtx_plot<-plot_measure_genome_wide(xtx,"M_XtX")+ylab("XTX")+
  geom_hline(data=xtx_cutoffs,aes(yintercept = cutoff))

# Plot these supps
pdf("figs/Figure_SX_Genome_scan_results.pdf",width=12,height=10)
afd_plot
xpehh_plot
xtx_plot
dev.off()

# Also save them as a list to use elsewhere
saveRDS(list(afd_plot,xpehh_plot,xtx_plot),"outputs/genome_scan_figures.rds")

#####################################################################################
# Plot chr20 to highlight aripo peak...
chr20_xtx_data <- xtx[xtx$chr=="chr20",]
chr20_xtx_data$river_F <- factor(chr20_xtx_data$river,levels = five_aside_colour_rivers$river)
chr20_xtx <- ggplot(chr20_xtx_data,aes(mid,M_XtX,colour=river_F))+
  geom_point(alpha=0.5)+
  facet_wrap(~river_F,ncol=1,strip.position = "right")+
  scale_colour_manual(values=five_aside_colour_rivers$colour,
                      breaks = five_aside_colour_rivers$river_F)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        strip.text = element_text(size=12))+
  scale_x_continuous(breaks=seq(0,max(chr20_xtx_data$mid),2e6),
                     labels=seq(0,max(chr20_xtx_data$mid),2e6)/1e6)+
  labs(y="XTX",x="Chromosome 20 Position (Mb)")

# Plot to output
pdf("figs/FigureSX_chr20_xtx_zoom.pdf",width=9,height=6)
chr20_xtx
dev.off()

#####################################################################################
# For each measure, plot density distribution of each as well
density_dd<-data.frame(river=c(AFD$river,xpehh$river,xtx$river),
                       measure=c(AFD$AFD,abs(xpehh$normxpehh),xtx$M_XtX),
                       group=c(rep("AFD",nrow(AFD)),rep("XP-EHH",nrow(xpehh)),rep("XTX",nrow(xtx))))
density_dd$river_F<-factor(density_dd$river,levels=rivers)

density_cutoffs<-rbind(afd_cutoffs,xtx_cutoffs,
                       data.frame(river=rivers,
                                  cutoff=2,
                                  river_F=afd_cutoffs$river_F))
density_cutoffs$group<-rep(c("AFD","XTX","XP-EHH"),each=5)

scan_density<-lapply(unique(density_dd$group),function(x){
tmp <- ggplot(density_dd[density_dd$group == x,],aes(x=measure,fill=river_F))+
  geom_density(alpha=0.6)+
  facet_wrap(~river_F,scales = "free_y",ncol=1,strip.position = "right")+
  geom_vline(data=density_cutoffs[density_cutoffs$group == x,],aes(xintercept=cutoff))+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.position = "none")+
  scale_fill_manual(values=five_aside_colour_rivers$colour,
                    breaks=five_aside_colour_rivers$river_F)+
  ylab("Density")+
  xlab(x)
return(tmp)
})

# plot them all together
pdf("figs/Figure_SX_Distributions_of_genome_scans.pdf",width=12,height=8)
plot_grid(plotlist = scan_density,
          ncol=3,labels = "AUTO")
dev.off()

#####################################################################################
# For each river, pull the outliers that are relevant to that river
outlier_pool<-lapply(rivers,function(x){
  
  afd_tmp<-data.frame(window_id=AFD[AFD$river == x & AFD$outlier %in% c("95","99"),"window_id"],
                      measure="AFD")
  xtx_tmp<-data.frame(window_id=xtx[xtx$river == x & xtx$outlier %in% c("95","99","999"),"window_id"],
                      measure="XTX")
  xpehh_tmp<-data.frame(window_id=xpehh[xpehh$river == x & xpehh$outlier %in% c("95","99"),"window_id"],
                        measure="XPEHH")
  
  # Merge these all together
  outliers_tmp<-rbind(afd_tmp,xtx_tmp,xpehh_tmp)
  
  # Add river
  outliers_tmp$river<-x
  
  # Count
  counts<-table(outliers_tmp$window_id)
  to_keep<-outliers_tmp[outliers_tmp$window_id  %in% names(counts[counts >= 2]),]
  
return(list(to_keep,unique(to_keep$window_id)))
})

# Fetch the data.frames
outlier_dd<-data.frame(rbindlist(lapply(outlier_pool,function(x){x[[1]]})))
outlier_dd <- outlier_dd %>% separate(window_id,into = c("chr","window"),sep = ':')
outlier_dd <- outlier_dd %>% separate(window,into = c("BP1","BP2"),sep = '-')
outlier_dd$BP1<-as.integer(outlier_dd$BP1)
outlier_dd$BP2<-as.integer(outlier_dd$BP2)
outlier_dd$window_id<-paste0(outlier_dd$chr,":",outlier_dd$BP1,"-",outlier_dd$BP2)
saveRDS(outlier_dd,
        paste0("outputs/full_outlier_info_for_all_rivers_wind",wind_size,".rds"))

# Fetch all the outliers
river_outliers<-lapply(outlier_pool,function(x){x[[2]]})
names(river_outliers)<-rivers
saveRDS(river_outliers,
        paste0("outputs/outlier_candidate_regions_for_all_rivers_",wind_size,".rds"))

# Plot the overlap as an Upset
outlier_upset<-upset(fromList(river_outliers),
                     sets = rev(five_aside_colour_rivers$river),
                     keep.order = T,
                     sets.bar.color = rev(five_aside_colour_rivers$colour),
                     
                     text.scale = 1.8,
                     number.angles = -45,
                     set_size.angles = 45,
                     mb.ratio = c(0.5,0.5),
                     point.size = 3
                     )

saveRDS(outlier_upset,"outputs/outlier_upset_figure.rds")

# Count up overlap
all_outliers<-unlist(river_outliers)
count<-table(all_outliers)
max(count)

# Count overlapping outliers that occur in more than river
overlapping_windows<-data.frame(window_id=names(count[count >= 2]))
saveRDS(overlapping_windows,
        paste0("outputs/AFD,XTX,XPEHH_overlapping_outliers_",wind_size,".rds"))
overlapping_windows_chr <- overlapping_windows %>% separate(window_id,into = c("chr","window"),sep = ":")

#####################################################################################
# Summarise these windows into a supp table
overlapping_windows$rivers_count <- NA
overlapping_windows$rivers <- NA
overlapping_windows$rivers_scans <- NA
for(i in 1:nrow(overlapping_windows)){
  river_tmp <- na.omit(sapply(rivers,function(x){
    if(overlapping_windows$window_id[i] %in% river_outliers[[x]]){
      return(x)
    } else {
      return(NA)
    }
  }))
  
  # For each, get the line of evidence as well
  river_evidence <- sapply(river_tmp,function(x){
    outlier_sub <- outlier_dd[outlier_dd$river == x &
                                outlier_dd$window_id == overlapping_windows$window_id[i],]
    return(paste0(x,"-",paste(outlier_sub$measure,collapse = "/")))
  })
  
  overlapping_windows$rivers[i] <- paste(river_tmp,collapse = " & ")
  overlapping_windows$rivers_count[i] <- length(river_tmp)
  overlapping_windows$rivers_scans[i] <- paste(river_evidence,collapse = " & ")
}

overlapping_windows <- overlapping_windows[order(-overlapping_windows$rivers_count),]

# Write the table
write.table(overlapping_windows,
            "tables/TableSX_overlapping_outlier_scans_among_rivers.txt",
            row.names = F,sep = "\t",quote = F)

#####################################################################################
# Summarise numbers of outlier types in table...
outlier_counts <- matrix(ncol=8,nrow=5)
outlier_counts[,1] <- rivers
colnames(outlier_counts) <- c("River","AFD","XTX", "XPEHH","AFD & XTX","AFD & XPEHH","XTX & XPEHH", "AFD & XTX & XPEHH")
# outlier_overlaps <- list("AFD & XTX" = c("AFD","XTX"),
#                          "AFD & XPEHH" = c("AFD","XPEHH"),
#                          "XTX & XPEHH" = c("XTX","XPEHH"),
#                          "AFD & XTX & XPEHH" = c("AFD","XTX","XPEHH"))

# Fill matrix
for(i in 1:5){
  
  outlier_river <- outlier_dd[outlier_dd$river == rivers[i],]
  
  # Get grouped vars
  outlier_types <- sapply(unique(outlier_river$window_id),function(x){
    scans_tmp <- outlier_river[outlier_river$window_id == x,"measure"]
    return(paste(scans_tmp,collapse = " & "))
  })
  
  # Count types
  measure_overlap_counts <- table(outlier_types)
  
  # Fill
  for(overlap in names(measure_overlap_counts)){
    outlier_counts[i,overlap] <- measure_overlap_counts[overlap]
  }
  
  # Get single counts
  single_counts <- table(outlier_river$measure)
  for(measure in names(single_counts)){
    outlier_counts[i,measure] <- single_counts[measure]
  }
}

# Save this
write.table(outlier_counts,
            "tables/TableSX_outlier_and_overlaps.txt",
            row.names = F,quote = F,sep = "\t")

#####################################################################################
# Also look for "semi-overlapping" windows
all_outliers<-lapply(outlier_pool,function(x){x[[2]]})
names(all_outliers)<-rivers

# Visualise
upset(fromList(all_outliers))
semi_overlap<-semi_overlap_scan(all_outliers,100000,50000)
semi_overlap<-semi_overlap[semi_overlap$overlapping > 1,]
nrow(semi_overlap[semi_overlap$overlapping > 2,])

#####################################################################################
# Confirm there are no fully overlapping outliers for any one statistic...
# AFD
AFD_outliers<-AFD[AFD$outlier %in% c(95,99),]
AFD_count<-table(AFD_outliers$window_id)
max(AFD_count)
AFD_max<-AFD_outliers[AFD_outliers$window_id %in% names(AFD_count[AFD_count == max(AFD_count)]),]
AFD_max[order(AFD_max$window_id),]

# And semi...
AFD_outlier_winds<-lapply(rivers,function(x){return(AFD_outliers[AFD_outliers$river == x,"window_id"])})
names(AFD_outlier_winds)<-rivers
semi_AFD<-semi_overlap_scan(AFD_outlier_winds,100000,50000)
semi_AFD[semi_AFD$overlapping == max(semi_AFD$overlapping),]

# xtx
XTX_outliers<-xtx[xtx$outlier != "No",]
XTX_count<-table(XTX_outliers$window_id)
max(XTX_count)
XTX_max<-XTX_outliers[XTX_outliers$window_id %in% names(XTX_count[XTX_count == max(XTX_count)]),]
XTX_max[order(XTX_max$window_id),]

# And semi...
XTX_outlier_winds<-lapply(rivers,function(x){return(XTX_outliers[XTX_outliers$river == x,"window_id"])})
names(XTX_outlier_winds)<-rivers
semi_XTX<-semi_overlap_scan(XTX_outlier_winds,100000,50000)
semi_XTX[semi_XTX$overlapping == max(semi_XTX$overlapping),]

# xpehh
xpehh_outliers<-xpehh[xpehh$outlier != "No",]
xpehh_count<-table(xpehh_outliers$window_id)
max(xpehh_count)
xpehh_max<-xpehh_outliers[xpehh_outliers$window_id %in% names(xpehh_count[xpehh_count == max(xpehh_count)]),]
xpehh_max[order(xpehh_max$window_id),]

# And semi...
XPEHH_outlier_winds<-lapply(rivers,function(x){return(xpehh_outliers[xpehh_outliers$river == x,"window_id"])})
names(XPEHH_outlier_winds)<-rivers
semi_XPEHH<-semi_overlap_scan(XPEHH_outlier_winds,100000,50000)
semi_XPEHH[semi_XPEHH$overlapping == max(semi_XPEHH$overlapping),]

################################################
# Also get the BayPass results...
baypass_fig<-readRDS("outputs/baypass_figures_genome_chr20,8.rds")
pathway_fig<-readRDS("outputs/pathway_plot_figure.rds")


# Fix the upset
outlier_upset2 <- cowplot::plot_grid(NULL, outlier_upset$Main_bar, outlier_upset$Sizes, outlier_upset$Matrix,
                           nrow=2, align='hv', rel_heights = c(1.5,1),
                           rel_widths = c(2,3))

# Plot together
pdf("figs/Figure2_selection_scan_overlap_pathway_baypass.pdf",width=16,height=10)
plot_grid(plot_grid(outlier_upset2,pathway_fig,ncol=2,labels = c("A","B"),rel_widths = c(1.5,1),label_size = 20),
          #plot_grid(baypass_fig[[1]],label=c("C)"),label_size = 20),
          baypass_fig[[1]]+theme(axis.title.y = element_text(size=14)),
          baypass_fig[[2]]+theme(axis.title.y = element_text(size=14)),
          labels = c("","C","D"),label_size=20,
          nrow=3,ncol=1,rel_heights = c(1.2,0.8,0.8),
          align = "v",axis = "tblr",
          label_y=1.1)
dev.off()

################################################
# Associations between Baypass and selection scans
baypass_outliers <- readRDS("outputs/baypass_AUXmodel_windows_10000_processed_NoMerge.rds")
baypass_outliers <- baypass_outliers[baypass_outliers$outlier == 999,"window_id"]

# Make a dataframe of outliers from selection scanning
river_outliers_dd <- data.frame(window_id = unlist(river_outliers),
                              river = c(rep("Tacarigua",length(river_outliers[[1]])),
                                        rep("Guanapo",length(river_outliers[[2]])),
                                        rep("Aripo",length(river_outliers[[3]])),
                                        rep("Oropouche",length(river_outliers[[4]])),
                                        rep("Madamas",length(river_outliers[[5]]))))


# Intersect them...
baypass_intersect <- data.frame(rbindlist(lapply(baypass_outliers,function(window){
  
  # Get rivers
  tmp_rivers <- river_outliers_dd[river_outliers_dd$window_id == window,"river"]
  
  # Make out
  out <- data.frame(window_id = window)
  
  if(length(tmp_rivers) == 0){
    out$rivers <- "None"
  } else {
    out$rivers <- paste0(tmp_rivers,collapse = ", ")
  }

  # Now do individuals
  afd_tmp <- paste(na.omit(unlist(lapply(rivers,function(x){
    
    tmp2 <- AFD_outlier_winds[[x]][AFD_outlier_winds[[x]] == window]
    if(length(tmp2) == 0){
      return(NA)
    } else {return(x)}
  }))),collapse=", ")
  
  xtx_tmp <-   paste(na.omit(unlist(lapply(rivers,function(x){
    
    tmp2 <- XTX_outlier_winds[[x]][XTX_outlier_winds[[x]] == window]
    if(length(tmp2) == 0){
      return(NA)
    } else {return(x)}
  }))),collapse=", ")
  
  
  xpehh_tmp <-   paste(na.omit(unlist(lapply(rivers,function(x){
    
    tmp2 <- XPEHH_outlier_winds[[x]][XPEHH_outlier_winds[[x]] == window]
    if(length(tmp2) == 0){
      return(NA)
    } else {return(x)}
  }))),collapse=", ")

# Combine
out$afd_outliers <- afd_tmp
out$xtx_outliers <- xtx_tmp
out$xpehh_outliers <- xpehh_tmp
return(out)
})))


# Write to output
write.table(baypass_intersect,
            paste0("outputs/TableSX_baypass_selection_scan_intersection.txt"),
            quote = F,row.names = F,sep="\t")

