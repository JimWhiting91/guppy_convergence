################################################################################################
# Based on HiC information, merge scaffold 94 with chr20 

merge_scaf94_chr20<-function(x,scaf=NULL){
  
  ## Example usage
  # merge_scaf94_chr20(dd[dd$chr == "chr20","BP"],scaf="chr20")
  # merge_scaf94_chr20(dd[dd$chr == "000094F_0","BP"],scaf="94")
  
  merged_out<-sapply(x,function(bp){
    
    # Is it on SCF_204?
    if(scaf=="chr20"){
      if(bp > 836423 & bp <= 3164071){
        bp_out<- 3164071 - bp + 836423 + 1797025
      } else if (bp < 836423) {
        bp_out<-bp
      } else {
        bp_out<-bp+1797025
      }
    } else if (scaf=="94"){
      bp_out<-1797025 - bp + 836423
    }
    return(bp_out)
  })
  return(as.integer(merged_out))
}

################################################################################################
# Transform BUCK co-ords to STAR positions for chr15 region
buck2star_chr15<-function(x){
  
  chr15_out <-sapply(x,function(bp){

      if(bp > 4913036 & bp <= 5120295){ # Main bit
        bp_out<- bp + 29701
      } else if (bp > 1610024 & bp <= 1642740) {
      	# Add on and invert
        bp_out<- 4932724 - (bp-1610024)
        # And invert       
      } else {
      	bp_out<-NA
      }

    return(bp_out)
  })
  return(chr15_out)
}


################################################################################################

make_windows<-function(data,colname,minBP,maxBP,windsize,avg){
	
  if(!is.character(avg) | !is.character(colname)){
  	stop("ERROR: Either colname or avg is not a character string")
  }

  winds1<-seq(minBP,maxBP,windsize)
  winds2<-winds1+windsize

  # Summarise for each
  wind_sum<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){

    tmp<-data

    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]

    if(nrow(tmp2) > 0){

    # Mean or median
    if(avg == "mean"){
    out<-data.frame(BP1=winds1[y],
                    BP2=winds2[y],
                     sum=mean(tmp2[,colname]))
    } else if (avg == "median"){
    out<-data.frame(BP1=winds1[y],
                    BP2=winds2[y],
                    sum=median(tmp2[,colname]))
    } else if (avg == "count"){
      out<-data.frame(BP1=winds1[y],
                      BP2=winds2[y],
                      sum=length(tmp2[,colname]))
    }

   colnames(out)[3]<-colname

    return(out)
    } else {
      out<-data.frame(BP1=winds1[y],
                      BP2=winds2[y],
                      sum=NA)
    colnames(out)[3]<-colname
      return(out)
    }
  },mc.cores=detectCores()-1)))
  return(wind_sum)
}

################################################################################################
# Calculate allele frequencies using count matrices

calc_AF<-function(freqs){
  
  # Function assumes an N x 4 matrix, 2 cols per pop, N SNPs
  freqs<-cbind(freqs,snps[,1:2])
  
  # First we remove any rows in which allele frequency is 1 in both pops or 0 as these are "invariant"
  freqs2<-freqs[!(freqs[,1] == 0 & freqs[,3] == 0),]
  freqs2<-freqs2[!(freqs2[,2] == 0 & freqs2[,4] == 0),]
  
  # Calculate AF for each SNP per pop
  freqs2$freq1<-freqs2[,1]/rowSums(freqs2[,1:2])
  freqs2$freq2<-freqs2[,3]/rowSums(freqs2[,3:4])
  
  # Calculate AFD
  freqs2$AFD<-abs(freqs2$freq1-freqs2$freq2)
  
  # Return the output
  return(freqs2[,5:ncol(freqs2)])
}

################################################################################################
# Function for sliding window assessment of close but not fully-overlapping outliers

semi_overlap_scan<-function(windows,window_size,sliding=window_size){
  # ------------------
  # Input info
  # Windows should be a named list of outliers in form of chrX:XXX-XXX
  # window_size = how large you want scan windows to be, eg. windows of 100kb = 100000
  # sliding = how windows move across the genome, eg. 50000 to move in 50kb increments. For non-overlapping, sliding = window_size. This is default
  # ------------------
  
  # Load these
  library(tidyr)
  library(parallel)
  
  # Merge windows to a data.frame
  outliers<-data.frame(window_id=unlist(windows))
  outliers$pop<-gsub('[[:digit:]]+', '', rownames(outliers))
  outliers<-outliers %>% separate(window_id,into = c("chr","window"),sep = ":")
  outliers<-outliers %>% separate(window,into = c("BP1","BP2"),sep = '-')
  outliers$BP1<-as.integer(outliers$BP1)
  outliers$BP2<-as.integer(outliers$BP2)
  outliers$mid<-rowMeans(outliers[,c("BP1","BP2")])
  
  # Make windowed input
  tmp<-data.frame(window=unlist(windows))
  tmp <- tmp %>% separate(window,into = c("chr","window"),sep = ":")
  tmp$BP2 <- sapply(tmp$window,function(x){return(strsplit(x,"-")[[1]][2])})
  
  # Make windows over each chr
  chrs<-unique(tmp$chr)
  scan_windows<-data.frame(rbindlist(lapply(chrs,function(x){
    maxBP<-max(as.integer(tmp[tmp$chr == x,"BP2"]))
    winds<-seq(0,maxBP,sliding)
    out<-data.frame(chr=x,
                    BP1=as.integer(winds),
                    BP2=as.integer(winds+window_size))
    return(out)
  })))
  
  # Now run through all the windows and harvest the number of outliers
  overlapping<-na.omit(data.frame(rbindlist(mclapply(1:nrow(scan_windows),function(x){
    print(x)
    
    # Subset
    outlier_sub<-outliers[outliers$chr == scan_windows$chr[x] &
                            outliers$mid > scan_windows$BP1[x] &
                            outliers$mid < scan_windows$BP2[x], ]
    
    # If empty return nothing
    if(nrow(outlier_sub) == 0){
      out<-data.frame(window_id=paste0(scan_windows$chr[x],":",as.integer(scan_windows$BP1[x]),"-",as.integer(scan_windows$BP2[x])),
                      overlapping=NA,
                      pops=NA)
      return(out)
    } else {
    
    # Which populations?
    count<-length(unique(outlier_sub$pop))
    pops<-paste(sort(unique(outlier_sub$pop)),collapse = ",")
    out<-data.frame(window_id=paste0(scan_windows$chr[x],":",as.integer(scan_windows$BP1[x]),"-",as.integer(scan_windows$BP2[x])),
                    overlapping=count,
                    pops=pops)
    return(out)
    }
  },mc.cores=detectCores()-1))))

return(overlapping)
}

######################################################################################################
# Function takes a set of genome regions, maps back to genome and returns genes with NCBI, Uniprot and KEGG IDs etc
# Fetch the Uniprot IDs...
ensembl_digger<-function(outliers=NULL,genome1=NULL,genome2=NULL,aln_block_length=1000,aln_qual=40,biomart_attributes=c("ensembl_gene_id","entrezgene_id","uniprot_gn_id","kegg_enzyme","chromosome_name","start_position","end_position"),n.cores=1,active_mart=guppy){
          
  # Firstly, take all the outlier regions and make a multi-fasta
  lapply(outliers,function(x){
  system(paste0("samtools faidx ",genome1," ",x," >> outputs/tmp_multi.fa"))
  })
  
  # Align to old genome
  system2('/Users/jimwhiting/bin/minimap2',
          args=c(paste0("-t ",n.cores),genome2,"outputs/tmp_multi.fa"),
          stdout ="outputs/tmp_multi.paf",wait=T)
            
  # Fetch regions
  aln<-read.table("outputs/tmp_multi.paf",fill=T)

  # Keep tidy!
  system(paste0("rm -f outputs/tmp_multi*"))

  # Only keep "high support alignment regions"
  aln<-aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]
  regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))

  # Pull uniprot genes from biomaRt for each region
		tmp_biomart<-getBM(attributes = biomart_attributes,
                   filters= "chromosomal_region",
                   values=regions,
                   mart=active_mart)
	
  # Return
return(tmp_biomart)
 }

#########################################################################################################
# Five Aside colour scheme...
library(RColorBrewer)
#five_aside_colour_pops<-data.frame(pop = c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","TACHP","TACLP"),
#								   colour = c("deepskyblue4","deepskyblue1"))

five_aside_colour_rivers<-data.frame(river = c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"),
									 colour = c(brewer.pal(n = 8, name = "Dark2")[1:3],
									 			brewer.pal(n = 8, name = "Dark2")[6],
									 			brewer.pal(n = 8, name = "Dark2")[8])) 
five_aside_colour_rivers$river_F<-factor(five_aside_colour_rivers$river,levels= five_aside_colour_rivers$river)