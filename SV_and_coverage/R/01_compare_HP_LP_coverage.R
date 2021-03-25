###########################################################################################################################################################################
# This script takes a directory of deeptools coverage bams merged over populations and averaged into 1kb windows and examines coverage and coverage ratios between HP and LP pops
###########################################################################################################################################################################

## Load packages
lib = c("parallel","data.table","tidyr","dplyr","ggplot2")
lapply(lib, library, character.only=TRUE)

# Get parameters from commandline to find outputs
setwd("/gpfs/ts0/projects/Research_Project-T110748/people/jimbo/HP_LP_trials/deeptools")
input<-"/gpfs/ts0/home/jw962/HP_LP_trials/deeptools/outputs/coverage/five_aside_STAR_FINAL"
  
# Define pops
pops<-c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
pairs<-lapply(seq(1,length(pops),2),function(x){return(pops[c(x,x+1)])})
names(pairs)<-c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
rivers<-rep(names(pairs),each=2)

# Define window size
WINDOW=10000

# List all the files
inputs<-list.files(input)

# Read in all the data to 1 file
dd<-data.frame(rbindlist(lapply(1:length(pops),function(x){

# Get files
to_read<-grep(pops[x],inputs,value=T)
to_read<-grep(paste0(WINDOW,"_"),to_read,value=T)

# Read in chr at a time
tmp_dd<-data.frame(rbindlist(lapply(to_read,function(y){

# Read file
tmp<-data.frame(fread(paste0(input,"/",y)))

# Rename
colnames(tmp)<-c("chr","BP1","BP2","coverage")

# Add cols
tmp$pop<-pops[x]

if(x %in% seq(1,10,2)){
tmp$pred<-"High"
} else {
tmp$pred<-"Low"
}

tmp$river<-rivers[x]

return(tmp)
})))

return(tmp_dd)
})))

######################################################################
# Do Plots and visualisations
chrs_to_plot<-read.table("~/guppy_research/STAR/STAR.chromosomes.release.fasta.fai")
chrs_to_plot<-chrs_to_plot[chrs_to_plot$V2 > 1000000,"V1"]
chrs<-unique(dd$chr)

# Plot a file per chromosome
ratio_out<-data.frame(rbindlist(mclapply(chrs,function(x){
	print(x)

tmp<-dd[dd$chr == x,]

# First do plot of raw data
raw_plot<-ggplot(tmp,aes(x=BP1,y=coverage))+
          geom_line()+
          facet_grid(river~pred)+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(tmp$BP2),2),by=2000000),
                       labels=seq(0,signif(max(tmp$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="none")+
          ylab("Coverage")

# Now we need to calculate coverage ratio per pair
ratios<-data.frame(rbindlist(lapply(1:length(pairs),function(y){

tmp2<-tmp[tmp$pop %in% pairs[[y]],]
BP_to_keep<-duplicated(tmp2$BP1)
tmp2<-tmp2[tmp2$BP1 %in% tmp2$BP1[BP_to_keep],]

high_cov<-tmp2[tmp2$pop == pairs[[y]][1],"coverage"]
low_cov<-tmp2[tmp2$pop == pairs[[y]][2],"coverage"]

out<-data.frame(chr=tmp2$chr,
                BP1=tmp2[tmp2$pop== pairs[[y]][1],"BP1"],
                BP2=tmp2[tmp2$pop== pairs[[y]][1],"BP2"],
                pop=tmp2[tmp2$pop== pairs[[y]][1],"pop"],
                river=tmp2[tmp2$pop== pairs[[y]][1],"river"],
                pred=tmp2[tmp2$pop== pairs[[y]][1],"pred"],
                cov_ratio=high_cov/low_cov,
                high_cov_est=high_cov,
                low_cov_est=low_cov)

return(out)
})))

if(nrow(ratios) > 0){

# Plot raw values as above as ln-transform
raw_ratio_plot<-ggplot(ratios,aes(x=BP1,y=log2(cov_ratio)))+
          geom_line()+
          facet_wrap(~river,ncol=1,strip.position="right")+
          xlab(x)+
          scale_x_continuous(breaks = seq(0,signif(max(ratios$BP2),2),by=2000000),
                       labels=seq(0,signif(max(ratios$BP2),2),by=2000000)/1000000)+
                          theme_bw()+
          theme(panel.grid = element_blank(),
          axis.title = element_text(size=24),
          axis.text = element_text(size=14),
          strip.text = element_text(size=12),
          legend.position="none")+
          ylab("Coverage Ratio (H/L)")

# Plot them both together if big enough
if(x %in% chrs_to_plot){
pdf(paste0("figs/five_aside_STAR_",x,"_",WINDOW,"_Coverage_Raw_figs.pdf"),width=16)
print(raw_plot)
print(raw_ratio_plot)
dev.off()
}

}

return(ratios)
},mc.cores=12)))

# Save the ratios
write.table(ratio_out,
            paste0("outputs/five_aside_STAR_coverage_ratios_",WINDOW,".txt"),row.names=F,quote=F,sep="\t")
            
