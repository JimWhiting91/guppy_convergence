##########################################################################################
# Script runs local PCA across a chromosomal PCA to infer windows of shared ancestry

# Version 2 of this script runs off a whole VCF and subsets in-line
#########################################################################################

# Set working directory
setwd("~/Exeter/guppy_convergence/localPCA")

## Load packages
lib = c("parallel","data.table","tidyr","dplyr","ggplot2","lostruct","ggpubr","vcfR")
lapply(lib, library, character.only=TRUE)

# This is the name of the original vcf
vcf<-"~/Exeter/VCFs/five_aside_STAR_3033083_final.vcf.gz"
window_N<-100
chr<-"chr20"

#####################################
# LOOP OVER POPS
all_pops<-c("GH","GL","LO","UQ","APHP","APLP","LT","UT","UMD","LMD")

lapply(all_pops,function(i){
#####################################

pops<-i

#output<-"6nations_GH_chr12"
output<-paste0("five_aside_STAR_",i,"_",chr)

# Make the input
system(paste0("bcftools view -r ",chr," ",vcf," > data/localPCA_tmp.vcf"),wait=TRUE)

# Read in VCF in popGenome for info on BP locations for unphased vcf
vcf_PG<-read.vcfR("data/localPCA_tmp.vcf")
BP1<-as.integer(vcf_PG@fix[,2])
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]
inds_vec<-colnames(vcf_PG@gt)[2:length(colnames(vcf_PG@gt))]

# Filter for inds we want
inds_to_keep<-unlist(lapply(pops,function(x){return(grep(x,inds_vec))}))
vcf_PG@gt<-vcf_PG@gt[,c(1,inds_to_keep+1)]
# Get males and female columns
#females<-match(female_vec,inds_vec)[is.na(match(female_vec,inds_vec))==FALSE]
#males<-setdiff(1:length(inds_vec),females)

# Return logical vector for polymorphic sites
invariant_filter<-is.polymorphic(vcf_PG,na.omit=T)

# Read in VCF for localPCA
vcf_in<-read_vcf("data/localPCA_tmp.vcf")

# And tidy up...
system("rm -f data/localPCA_tmp.vcf")

# Filter for invariants
vcf_in<-vcf_in[invariant_filter,inds_to_keep]
BP1<-BP1[invariant_filter]
window_BP<-BP1[seq((window_N/2),length(BP1),by=window_N)]

# Read in eigen_windows, win=number of rows of matrix, k = number of eigenvector/value pairs to return, mc.cores compatible but does seem to crash sometimes with mc.cores
eigenstuff <- eigen_windows(vcf_in, win=window_N, k=2)

# Calculate the distances, npc=N of PCs computed, mc.core compatiblke
windist <- pc_dist(eigenstuff, npc=2,mc.cores=detectCores()-1)

# Remove any windows for which we are getting NA values
to_keep_logical<-rowSums(is.na(windist)) != ncol(windist)
to_keep<-(1:ncol(windist))[to_keep_logical]
windist2<-windist[to_keep,to_keep]

# Principal Coordinates Analysis, k = dimension of space for data to represented in, eig = return eigenvalues
fit2d <- cmdscale( windist2, eig=TRUE, k=2 )
eigenvals <- fit2d$eig/sum(fit2d$eig)
eigenvals_10 <- eigenvals[eigenvals > 0.1]

# Plot along a chromosome
plot_dd<-data.frame(mds1=fit2d$points[,1],
                    mds2=fit2d$points[,2],
                    window_N=1:length(fit2d$points[,1]),
                    window_pos=window_BP[to_keep])

# Calculate cut-offs for outliers by trimming distributions and taking 3*SD of trimmed distribution
trim_mds1<-plot_dd[plot_dd$mds1 > quantile(plot_dd$mds1,probs=0.05) & plot_dd$mds1 < quantile(plot_dd$mds1,probs=0.95),"mds1"]
mds1_cutoffs<-c(mean(trim_mds1)+3*sd(trim_mds1),
                mean(trim_mds1)-3*sd(trim_mds1))

trim_mds2<-plot_dd[plot_dd$mds2 > quantile(plot_dd$mds2,probs=0.05) & plot_dd$mds2 < quantile(plot_dd$mds2,probs=0.95),"mds2"]
mds2_cutoffs<-c(mean(trim_mds2)+3*sd(trim_mds2),
                mean(trim_mds2)-3*sd(trim_mds2))

##################
# We also want to return the windows at the extremes of the distribution
##################

# Give window IDs
window_ids<-rep(0)
for(i in 2:length(to_keep)){
  window_ids[i]<-paste0(chr,":",((BP1[(i-1)*window_N])+1),"-",BP1[(i)*window_N])
}
window_ids[1]<-paste0(chr,":1-",BP1[window_N])
plot_dd$window_id <- window_ids[to_keep]


windows_to_keep1<-unique(c(plot_dd[plot_dd$mds1 < mds1_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds1 > mds1_cutoffs[1],"window_id"]))

windows_to_keep2<-unique(c(plot_dd[plot_dd$mds2 < mds2_cutoffs[2],"window_id"],
                           plot_dd[plot_dd$mds2 > mds2_cutoffs[1],"window_id"]))


# Reassign to plot_dd
plot_dd$outlier<-rep("N")
plot_dd[plot_dd$window_id %in% windows_to_keep2,"outlier"]<-"mds2"
plot_dd[plot_dd$window_id %in% windows_to_keep1,"outlier"]<-"mds1"

# Vals
cols1 <- c("mds1" = "red", "mds2" = "black", "N" = "black")
cols2 <- c("mds2" = "red", "mds1" = "black", "N" = "black")

plot_dd<-na.omit(plot_dd)

p1<-ggplot(plot_dd,aes(x=mds1,y=mds2, colour=window_N))+
  geom_point()+
  xlab("Coordinate 1")+
  ylab("Coordinate 2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.position="none",
        strip.text = element_text(size=24))+
  geom_hline(yintercept=mds2_cutoffs,linetype="dashed")+
  geom_vline(xintercept=mds1_cutoffs,linetype="dashed")    				

p2<-ggplot(plot_dd,aes(x=window_pos,y=mds1,colour=outlier))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,max(plot_dd$window_pos),by=1000000),
                     labels=seq(0,max(plot_dd$window_pos),by=1000000)/1000000)+
  ylab("MDS 1 Score")+
  xlab(paste0("Chr Pos (Mb)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=24),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=16,angle=90,hjust=1),
        strip.text = element_text(size=24),
        legend.position="none")+
  scale_colour_manual(values=cols1)

p2.2<-ggplot(plot_dd,aes(x=window_pos,y=mds2,colour=outlier))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,max(plot_dd$window_pos),by=1000000),
                     labels=seq(0,max(plot_dd$window_pos),by=1000000)/1000000)+
  ylab("MDS 2 Score")+
  xlab(paste0("Chr Pos (Mb)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=24),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=16,angle=90,hjust=1),
        strip.text = element_text(size=24),
        legend.position="none")+
  scale_colour_manual(values=cols2)


# Save figure
pdf(paste0("figs/",output,".pdf"),width=15)
print(ggarrange(p1,ggarrange(p2,p2.2,nrow=2),ncol=2,widths=c(2,3)))
dev.off()

# Output outlier windows
write.table(plot_dd,
            paste0("outputs/",output,".txt"),quote=F,row.names=F,sep="\t")

})
