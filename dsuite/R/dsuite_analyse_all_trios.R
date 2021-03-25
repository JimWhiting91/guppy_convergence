############################################################
# Analyse the Dsuite trios
lib<-c("data.table","ggplot2","viridis","tidyverse")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Dir name for outputs
dsuite_run <- "five_aside_STAR_v3_wingei_Guan_Tac_sisters_maf01_invariant_filtered_swap_TAC"

dd <- read.table(paste0("data/",dsuite_run,"/",dsuite_run,"_combined_trios_",dsuite_run,"_combined_trios_combined_tree.txt"),header=T)
dd_Dmin <- read.table(paste0("data/",dsuite_run,"/",dsuite_run,"_combined_trios_",dsuite_run,"_combined_trios_combined_Dmin.txt"),header=T)
dd_Dbbaa <- read.table(paste0("data/",dsuite_run,"/",dsuite_run,"_combined_trios_",dsuite_run,"_combined_trios_combined_BBAA.txt"),header=T)
dd<-dd[order(-dd$Dstatistic),]

dd$bonferroni<-p.adjust(dd$p.value,method = "fdr")

outliers<-dd[dd$bonferroni < 0.001,]
nrow(outliers)

pops<-list(c("GH","GL"),
           c("TACHP","TACLP"),
           c("APHP","APLP"),
           c("OH","OL"),
           c("MADHP","MADLP"))

river_res<-lapply(pops,function(x){
  
  rbind(outliers[outliers$P1 == x[[2]] & outliers$P2 == x[[1]],], outliers[outliers$P1 == x[[1]] & outliers$P2 == x[[2]],])
  
})

# Combine and save as a supp table...
within_river_outliers <- data.frame(rbindlist(river_res))
write.table(within_river_outliers,
            "tables/TableSX_within_river_D_outliers.tsv",
            col.names = T,row.names = F,quote = F,sep = "\t")

# Output all trios, ordered by f4-ratio
write.table(dd[order(-dd$f4.ratio),],
            "tables/TableSX_ordered_Dsuite_results.tsv",
            col.names = T,row.names = F,quote = F,sep = "\t")

                              
##### For significant trios, which chroms is introgression strongest on?? #####
chrs <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")
chrs <- chrs[chrs$V2 > 1000000,1]
inputs <- list.files(paste0("data/",dsuite_run))

# Take only outliers with an f4 ratio > 5%
within_river_outliers2 <- within_river_outliers[within_river_outliers$f4.ratio > 0.05,]

# Fetch the chrom_res
chrom_res <- data.frame(rbindlist(lapply(chrs,function(chr){
  print(chr)
  chr_inputs <- grep(paste0("_",chr,"_"),inputs,value = T)
  chr_inputs <- grep("scaf94",chr_inputs,invert = T,value = T)
  tmp <- read.table(paste0("data/",dsuite_run,"/",grep("tree.txt",chr_inputs,value=T)),header=T)
  tmp$chr = chr
  return(tmp)
})))

# Filter for trios of interest
trio_chrom_res <- lapply(1:nrow(within_river_outliers2),function(x){
  trio_focal <- chrom_res[chrom_res$P1 == within_river_outliers2$P1[x] & chrom_res$P2 == within_river_outliers2$P2[x] & chrom_res$P3 == within_river_outliers2$P3[x],]
  trio_focal[order(-trio_focal$f4.ratio),]
})

# Compare APLP and APHP introgression by chrom to see if they are associated
APHP_intro <- trio_chrom_res[[2]]
APHP_intro$pair <- "APHP-OHP"
APHP_intro <- APHP_intro[order(APHP_intro$chr),]
APLP_intro <- trio_chrom_res[[4]]
APLP_intro$pair <- "APLP-GHP"
APLP_intro <- APLP_intro[order(APLP_intro$chr),]

# Combine and plot
APHP_intro$APLP_pair <- APLP_intro$f4.ratio
ggplot(APHP_intro,aes(f4.ratio,APLP_pair))+
  geom_point()


