#################################################
# Get statistics on size of SVs from smoove

lib<-c("data.table","ggplot2","vcfR","dplyr","tidyverse")
lapply(lib,library,character.only=T)

# Fetch river VCFs
rivers <- c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
sv_vcf <- lapply(rivers,function(x){
  
  # Read it in
  vcf<-read.vcfR(paste0("~/Exeter/VCFs/",x,"_SV_filtered-smoove.genotyped.vcf.gz"))
  
  return(vcf)
})

# Counts of each kinds of SV
counts_matrix <- matrix(ncol=5,nrow=3)
rownames(counts_matrix)<-c("BND","DUP","INV")
colnames(counts_matrix)<-rivers
for(i in 1:5){
  fix<-data.frame(sv_vcf[[i]]@fix)
  svtype<-strsplit(fix$INFO,";")
  svtype<-sapply(svtype,function(x){return(x[1])})
  sv_counts<-table(svtype)
  counts_matrix[1,i]<-sv_counts["SVTYPE=BND"]
  counts_matrix[2,i]<-sv_counts["SVTYPE=DUP"]
  counts_matrix[3,i]<-sv_counts["SVTYPE=INV"]
}

# Because of low counts of inversions, lets look at those specifically...
inv_vcfs<-lapply(sv_vcf,function(x){
  fix<-data.frame(x@fix)
  svtype<-strsplit(fix$INFO,";")
  svtype<-sapply(svtype,function(x){return(x[1])})
  return(x[which(svtype=="SVTYPE=INV")])
})

# Fetch the location and size of each inversion in each population
inv_metadata<-lapply(inv_vcfs,function(x){
  fix<-data.frame(x@fix)
  svlen<-strsplit(fix$INFO,";")
  svlen<-sapply(svlen,function(x){return(gsub("SVLEN=","",x[2]))})
  
  inv_pos<-paste0(fix[,1],":",fix[,2])
  out<-data.frame(pos=inv_pos,
                  size=svlen)
  return(out)
})

# Now we'll check for SVs in the target regions around chr20 and scaf94...
candidate_vcf<-lapply(sv_vcf,function(x){
  fix<-data.frame(x@fix)
  to_keep<-which(fix[,1] %in% c("chr20","000094F_0") & as.integer(fix[,2]) < 4000000)
  return(x[to_keep,])
})
