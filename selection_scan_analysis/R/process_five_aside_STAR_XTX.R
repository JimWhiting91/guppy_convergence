require(corrplot) ; require(ape) ; require(data.table) ; require(dplyr) ; require(tidyverse) ; require(parallel)
#source the baypass R functions (check PATH)

source("R/baypass_utils.R")
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in info for SNPsr
snps<-data.frame(fread("data/AF/five_aside_STAR_snps",header = F))
snps$index<-1:nrow(snps)
chr_snps<-snps[grep("chr",snps$V1),] %>% separate(V1,sep="_",into=c("chr","BP"))
scaf_snps<-snps[grep("chr",snps$V1,invert=T),]
scaf_snps$V1<-gsub("F_0","F0",scaf_snps$V1)
scaf_snps<- scaf_snps %>% separate(V1,sep="_",into=c("chr","BP"))
scaf_snps$chr<-gsub("F0","F_0",scaf_snps$chr)
snps<-rbind(chr_snps,scaf_snps)

# Read in genos to filter invariants
genos<-data.frame(fread("data/AF/five_aside_STAR.geno"))

# List rivers
rivers<-c("TAC","G","AP","O","MAD")
river_names<-c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
names(river_names)<-rivers
# First for each river we need to know which sites are and are not invariants
invariant_filter<-function(x){
  freqs<-cbind(x,snps)
  freqs2<-freqs[!(freqs[,1] == 0 & freqs[,3] == 0),]
  freqs2<-freqs2[!(freqs2[,2] == 0 & freqs2[,4] == 0),]
  return(freqs2[,5:6])
}

# Make list of variant sites for each river
river_variants<-mclapply(seq(1,17,4),function(i){invariant_filter(genos[,i:(i+3)])},mc.cores=5)
names(river_variants)<-rivers

# Set window size
window_size<-50000

# Run over all rivers
baypass_list<-lapply(rivers,function(river){
  
  #Estimates of the XtX differentiation measures
  dd<-data.frame(rbindlist(lapply(1:16,function(x){
    # Get the subset SNPs
    snp_sub<-data.frame(fread(paste0("data/xtx/five_aside_STAR_sub",x,".snps")))
    snp_sub$index<-1:nrow(snp_sub)
    chr_snp_sub<-snp_sub[grep("chr",snp_sub$x),] %>% separate(x,sep="_",into=c("chr","BP"))
    scaf_snp_sub<-snp_sub[grep("chr",snp_sub$x,invert=T),]
    scaf_snp_sub$x<-gsub("F_0","F0",scaf_snp_sub$x)
    scaf_snp_sub<- scaf_snp_sub %>% separate(x,sep="_",into=c("chr","BP"))
    scaf_snp_sub$chr<-gsub("F0","F_0",scaf_snp_sub$chr)
    snp_sub<-rbind(chr_snp_sub,scaf_snp_sub)
    
    tmp<-data.frame(fread(paste0("data/xtx/",river,"_five_aside_STAR_sub",x,"_summary_pi_xtx.out")))
    tmp<-cbind(tmp,snp_sub[,c("chr","BP")])
    return(tmp)
  })))

  # Reorder
  dd$BP<-as.integer(dd$BP)
  dd<-dd[order(dd$chr,dd$BP),]
  
  # Filter out the invariants
  variants<-river_variants[[river]]
  dd<-dd[paste0(dd$chr,"_",dd$BP) %in% paste0(variants$chr,"_",variants$BP),]

  
  #######################################################
  # Fetch the XtX calibrations
  #######################################################
  
  #get the pod XtX
  pod.xtx=read.table(paste0("data/xtx/",river,"_five_aside_STAR_10k_sim_summary_pi_xtx.out"),h=T)$M_XtX
  #compute the thresholds
  pod.thresh95=quantile(pod.xtx,probs=0.95)
  pod.thresh99=quantile(pod.xtx,probs=0.99)
  pod.thresh999=quantile(pod.xtx,probs=0.999)
  
  
  ####################################################
  # Highlight 'outlier' Windows
  ####################################################
  # Mark outliers in main file
  dd$outliers<-"No"
  dd[dd$M_XtX > pod.thresh95,"outliers"]<-"Yes"
  
  # Build windows
  chrs<-unique(dd$chr)
  dd_winds<-data.frame(rbindlist(lapply(chrs,function(chr){
    tmp<-dd[dd$chr == chr,]
    winds<-make_windows(tmp,"M_XtX",0,max(tmp$BP),window_size,avg = "mean")
    winds$chr<-chr
    return(winds)
  })))

  # Filter duds
  baypass_winds<-na.omit(dd_winds)
  
  # Mark outliers
  baypass_winds$outlier<-"No"
  baypass_winds[baypass_winds$M_XtX > pod.thresh95,"outlier"]<-"95"
  baypass_winds[baypass_winds$M_XtX > pod.thresh99,"outlier"]<-"99"
  baypass_winds[baypass_winds$M_XtX > pod.thresh999,"outlier"]<-"999"
  
  # Mark river
  baypass_winds$river<-river_names[river]
  
  # And return
  return(baypass_winds)
  
  # Practise plot
  # ggplot(baypass_winds[baypass_winds$chr == "chr20",],aes(x=BP1,y=M_XtX))+
  #   geom_line()
  
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
  #       outlier_N<-nrow(tmp2[tmp2$outliers == "Yes",])
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
  # baypass_winds<-na.omit(dd_winds)
  # 
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
  # #
  # 
  # write.table(top_outliers,
  #             paste0("outputs/five_aside_STAR_",river,"_XTX_",window_size,"_TOP_OUTLIERS.txt"),
  #             row.names=F,quote=F,sep="\t")
  # 
  # write.table(top_outliers,
  #             paste0("outputs/five_aside_STAR_",river,"_XTX_",window_size,"_ALL_WINDOWS.txt"),
  #             row.names=F,quote=F,sep="\t")
  
})

# Save as RDS
saveRDS(baypass_list,
        paste0("outputs/processed_baypass_xtx_all_pops_window",window_size,".rds"))


