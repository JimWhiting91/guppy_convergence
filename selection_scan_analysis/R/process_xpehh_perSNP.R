###############################################
# Read in SelScan results for X-EHH and evaluate

# Get these
lib<-c("data.table","ggplot2","UpSetR","parallel")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Assign window size and variables
wind<-50000
pops<-c("TAC","G","AP","O","MAD")
rivers<-c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
chrs<-as.character(read.table("~/Exeter/VCFs/five_aside_STAR_chrs.txt")[,1])

# Read in results
dd<-data.frame(rbindlist(lapply(1:5,function(x){
  inputs<-list.files("data/xpehh/")
  inputs<-grep(pops[x],inputs,value=T)
  inputs<-grep("out.norm",inputs,value=T)
  inputs<-grep("windows",inputs,value=T,invert = T)
  
  # Read em in
  tmp<-data.frame(rbindlist(lapply(chrs,function(chr){
    tmp_input<-grep(paste0("_",chr,".xpehh"),inputs,value=T)
    #print(chr)
    if(length(tmp_input) > 0){
      if (file.size(paste0("data/xpehh/",tmp_input)) > 0){
        tmp<-data.frame(fread(paste0("data/xpehh/",tmp_input),header=T))
        if(nrow(tmp) > 0){
        tmp$chr<-chr
        
        # And make windows here
        tmp$BP<-tmp$pos
        winds<-make_windows(na.omit(tmp),"normxpehh",minBP = 0,maxBP = max(tmp$pos),avg="mean",windsize=wind)
        winds$chr<-chr
        winds$river=rivers[x]
        return(winds)
        }
      }
    }
  })))
  tmp$river<-rivers[x]
  return(tmp)
})))

# Make window IDS
dd$window_id<-paste0(dd$chr,":",as.integer(dd$BP1),"-",as.integer(dd$BP2))

# Look at xpehh
ggplot(dd,aes(x=normxpehh,fill=river))+
  #geom_density(alpha=0.5)+
  geom_histogram()+
  facet_wrap(~river,ncol=1,scales = "free_y")+
  geom_vline(xintercept = 2)

# Most significant xpehh windows per pop
outliers<-lapply(rivers,function(x){
  tmp<-dd[abs(dd$normxpehh) > 2  &  dd$river == x,]
  tmp<-na.omit(tmp[order(-tmp$normxpehh),"window_id"])
})

# Merge to dataframe
outlier_dd<-data.frame(window_id=unlist(outliers))
outlier_dd$river<-c(rep(names(outliers)[1],length(outliers[[1]])),
                    rep(names(outliers)[2],length(outliers[[2]])),
                    rep(names(outliers)[3],length(outliers[[3]])),
                    rep(names(outliers)[4],length(outliers[[4]])),
                    rep(names(outliers)[5],length(outliers[[5]])))

# Calculate overlap
overlap<-sort(table(unlist(outliers)),decreasing = T)
head(overlap,10)
names(outliers)<-rivers
upset(fromList(outliers),keep.order = F)

# Write the xtx results to Rdata
saveRDS(dd,paste0("outputs/xpehh_windows_window_",wind,".rds"))
