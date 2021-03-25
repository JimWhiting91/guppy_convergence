# ------------------------------
# From a file of allele counts per pop, calculate and assess AFD
lib<-c("data.table","ggplot2","dplyr","tidyr","Rfast","parallel")
lapply(lib,library,character.only=T)

source("~/Exeter/five_aside/five_aside_functions.R")

# Set up pops
pops<-c("TACHP","TACLP","GH","GL","APHP","APLP","OHP","OLP","MADHP","MADLP")

# Read in info for SNPs
snps<-data.frame(fread("data/AF/five_aside_STAR_snps",header = F))
snps$index<-1:nrow(snps)
chr_snps<-snps[grep("chr",snps$V1),] %>% separate(V1,sep="_",into=c("chr","BP"))
scaf_snps<-snps[grep("chr",snps$V1,invert=T),]
scaf_snps$V1<-gsub("F_0","F0",scaf_snps$V1)
scaf_snps<- scaf_snps %>% separate(V1,sep="_",into=c("chr","BP"))
scaf_snps$chr<-gsub("F0","F_0",scaf_snps$chr)
snps<-rbind(chr_snps,scaf_snps)

# Correct this
snps$BP<-as.integer(snps$BP)

# Fetch AFs - these are ordered according to pops
AF<-data.frame(fread("data/AF/five_aside_STAR.geno",header = F))

# For each pop calculate the AFDs
rivers<-c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
AFD_list<-mclapply(seq(1,20,4),function(x){
 tmp <- calc_AF(AF[,x:(x+3)])
 return(tmp)
},mc.cores=5)
names(AFD_list)<-rivers

# Now window-ise AFD for each pop and output the result to a data.frame
wind_size<-10000
AFD_winds<-data.frame(rbindlist(lapply(rivers,function(x){
  
  # Run window calculations
  chrs<-unique(AFD_list[[x]]$chr)
  winds<-na.omit(data.frame(rbindlist(lapply(chrs,function(y){
    tmp<-AFD_list[[x]][AFD_list[[x]]$chr == y,]
    tmp_winds<-make_windows(tmp,"AFD",0,max(tmp$BP),wind_size,avg = "median")
    counts<-make_windows(tmp,"AFD",0,max(tmp$BP),wind_size,avg = "count")
    tmp_winds$chr<-y
    tmp_winds$river<-x
    tmp_winds$SNP_count<-counts$AFD
    return(tmp_winds)
  }))))
  
  # Return
  return(winds)
})))

# Save these files as an object
AFD_results<-list(AFD_list,AFD_winds)
saveRDS(AFD_results,paste0("outputs/AFD_calculations_five_aside_STAR_",wind_size,".rds"))

# Check some summaries
AFD_winds %>% group_by(river) %>%
  summarise(AFD=median(AFD))

