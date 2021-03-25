# Visualise LP-LP SFS
library(ggplot2)
library(data.table)

# Set our pops
pop_order <- c("TACLP","GLP","APLP","OLP","MADLP")

# Fetch the SFS
sfs_files <- list.files("data/LP_LP_comparisons/")
sfs_files <- grep("TULP",sfs_files,invert = T,value = T)
sfs_dd <- data.frame(rbindlist(lapply(sfs_files,function(file){
  
  # Read it in
  sfs_in <- read.table(paste0("data/LP_LP_comparisons/",file),skip = 1,fill=T,header=T)
  #sfs_in <- sfs_in[3:nrow(sfs_in),2:ncol(sfs_in)]

  # Melt
  freq <- melt(as.matrix(sfs_in), value.name=c("Frq"))
  freq$Var1 <- as.integer(gsub("d1_","",freq$Var1))
  freq$Var2 <- as.integer(gsub("d0_","",freq$Var2))
  freq$Frq<- log(freq$Frq)
  freq<-freq[freq$Var2 < nrow(sfs_in)-freq$Var1,]

  # And set up pops 1 and 2, these are organised alphabetically
  pops <- sort(unlist(strsplit(file,"_"))[1:2])
  pop1 <- pops[1]
  pop2 <- pops[2]
  
  # Flip it if we need to..
  if(which(pop_order == pop1) > which(pop_order == pop2)){
    tmp <- freq$Var1
    freq$Var1 <- freq$Var2
    freq$Var2 <- tmp
    freq$pop1 <- pop2
    freq$pop2 <- pop1
  } else {
    freq$pop1 <- pop1
    freq$pop2 <- pop2
  }
  
return(freq)
  
})))

# Set order
sfs_dd$pop1_F <- factor(sfs_dd$pop1,levels=pop_order)
sfs_dd$pop2_F <- factor(sfs_dd$pop2,levels=pop_order)

# Plot
sfs_plot <- ggplot(sfs_dd, aes(Var2, Var1, fill = Frq))+
  geom_tile(size=0.2)+
  #scale_fill_gradient2(low = "white", high = "red3",midpoint = 5) +
  scale_fill_gradient2(low = "#00008b", high = "#ffff00", mid = "#4eee94",midpoint = 5.0) +
  # scale_fill_viridis()+
  theme_void()+ 
  theme(legend.position = c(0.1,0.1),
        legend.justification = c(0,0),
        axis.title = element_blank(),
        strip.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))+
  # theme(panel.grid = element_blank(), 
  #       axis.title.x = element_text(face = "bold",vjust = -1.2),
  #       axis.title.y = element_text(face = "bold"), 
  #       legend.position = "none",
  #       axis.text.x = element_text(size = 10, vjust=0.3, margin=margin(-12,0,0,0)),axis.text.y = element_text(size = 10, hjust=0.5, margin=margin(0,-12,0,0)))+
  # coord_fixed(ratio=1) +
  facet_grid(pop1_F~pop2_F,scales = "free")+
  labs(fill="Freq\n(log)")

# Save
pdf("figs/FigureSX_LP-LP_SFS_figs.pdf",width=8,height=6)
sfs_plot
dev.off()

