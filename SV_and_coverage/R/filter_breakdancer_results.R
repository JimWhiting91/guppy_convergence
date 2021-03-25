#################################################
# Filter breakdancer results

lib<-c("data.table","ggplot2","vcfR","dplyr","tidyverse","parallel")
lapply(lib,library,character.only=T)

# Fetch river VCFs
rivers <- c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
chrs <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")
chrs <- chrs[chrs$V2 > 1000000,1]

# Read in and filter
breakdancer_filtered <- lapply(rivers,function(x){
  
  # Run by chr
  chrs_out <- data.frame(rbindlist(mclapply(chrs,function(chr){
    
    # Read in results and info
    skipthese <- readLines(paste0("data/breakdancer_res/",x,"/",x,"_",chr,"_breakdancer.out"))
    skip_to <- max(grep("#",skipthese))
    res <- read.table(paste0("data/breakdancer_res/",x,"/",x,"_",chr,"_breakdancer.out"),skip = skip_to,fill=T)[,1:10]
    
    # Set colnames
    colnames(res)<-c("chr","pos","strand_rp1","chr2","pos2","strand_rp2","svtype","svlen","confidence","support")
    
    # Fetch filtering info
    out <- res[res$svlen > 1000 &
                 res$confidence == 99,]
    out <- out[out$support > median(out$support),]
    out$river<-x
    return(out)
  },mc.cores=detectCores()-1)))

return(chrs_out)
})

# Set names
names(breakdancer_filtered)<-rivers

# How many SVs per 
dd <- data.frame(rbindlist(breakdancer_filtered))
table(dd$river)

# How about SVs in our target region for each river...
candidate_svs <- lapply(rivers,function(x){
  
  tmp<-breakdancer_filtered[[x]]
  
  # Subset for region
  subset <- tmp[tmp$chr %in% c("chr20","000094F_0") &
                  tmp$pos < 4000000,]
  return(subset)
})

# Save the results as an R object..
saveRDS(breakdancer_filtered,
        "outputs/filtered_breakdancer_SVs_all_rivers.rds")

