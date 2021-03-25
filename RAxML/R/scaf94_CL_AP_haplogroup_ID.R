##### This script designates haplogroups in the Caroni drainage for the scaffold 94 inversion to determine age.
# For now, we assume the aripo-specific inversion is 0.5-1.1, and the whole scaffold is part of a larger inversion in Guanapo and Tac

# We need these
lib<-c("data.table","ggplot2","adegenet","vcfR","seqinr","ape","wesanderson","cowplot")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Fetch the VCF for 1 chr and
system("bcftools view -r 000094F_0 ~/Exeter/VCFs/five_aside_STAR_3033083_final.vcf.gz > outputs/scaf94_tmp.vcf",wait=TRUE)

# Read in VCFs
scaf94_vcf <- read.vcfR("outputs/scaf94_tmp.vcf")

# Filter
CL_AP_vcf <- scaf94_vcf[as.integer(scaf94_vcf@fix[,2]) > 500000 &
                          as.integer(scaf94_vcf@fix[,2]) < 1100000,]

############
# Define haplogroups based on DAPC of the inversion
############ 
# Convert
aripo_inversion.gen<-vcfR2genlight(CL_AP_vcf,n.cores=6)

# Filter for Caroni
inds<-row.names(aripo_inversion.gen$ind.names)

# Confirm 'inversion' structuring with glPCA
inversion.pca<-glPca(aripo_inversion.gen)
scores<-data.frame(inversion.pca$scores)
eigs<-data.frame(inversion.pca$eig)
eigs$vals<-round(eigs$inversion.pca.eig/sum(eigs$inversion.pca.eig),digits=2)

pops<-c("APHP","APLP","GH","GL","LMD","UMD","UQ","LO","LT","UT")
rivers<-c("Aripo","Guanapo","Madamas","Oropouche","Tacarigua")
rivers2<-rep(rivers,each=2)
HP<-c("GH","LO","APHP","LT","LMD")

scores$pop<-NA
scores$river<-NA
scores$pred<-NA
scores$ID<-rownames(scores)

# Add them in
for(i in 1:length(pops)){
  scores[grep(pops[i],scores$ID),"pop"]<-pops[i]
  scores[grep(pops[i],scores$ID),"river"]<-rivers2[i]
  
}

scores[scores$pop %in% HP,"pred"]<-"HP"
scores[!(scores$pop %in% HP),"pred"]<-"LP"

# Plot the main PCs
PC_plot<-ggplot(na.omit(scores),aes(PC1,PC2,colour=river,shape=pred))+
  geom_point(size=5,alpha=0.8)+
  scale_colour_manual(breaks = five_aside_colour_rivers$river,
                      values = five_aside_colour_rivers$colour)+
  scale_shape_manual(values=c(19,4))+
  theme_bw()+
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16))+
  labs(x=paste0("PC1 (",eigs$vals[1]*100,"%)"),
       y=paste0("PC2 (",eigs$vals[2]*100,"%)"),
       colour="River",shape="Predation")

# Save this
pdf("figs/CL_AP_region_PCA.pdf")
PC_plot+
  geom_vline(xintercept = c(-3,6),linetype="dashed")
dev.off()

# What SNPs are responsible for this
loadings<-data.frame(inversion.pca$loadings)
loadings$bp<-as.integer(CL_AP_vcf@fix[,2])
ggplot(loadings,aes(x=bp,y=Axis1))+geom_point()
ggplot(loadings,aes(x=bp,y=Axis2))+geom_point()

# Show histogram
ggplot(scores,aes(x=PC1,fill=pred))+
  geom_density(alpha=0.5)  +
  theme_bw()+
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16))+
  labs(fill="Predation",y="Density")

# Define Haplogroups. These cutoffs are based on the genotype plot below
haplogroup1<-scores[scores$PC1 < -3,"ID"]
haplogroup3<-scores[scores$PC1 > 6,"ID"]

# Run a genotype plot to check
source("~/software/genotype_plot.R")
popmap<-list(haplogroup1,haplogroup3)
names(popmap)<-c("haplogroup_REF","haplogroup_ALT")
genos<-genotype_plot(vcf="~/Exeter/VCFs/five_aside_STAR_3033083_final.vcf.gz",
                     chr = "000094F_0",
                     start=500000,
                     end=1100000,
                     popmap=popmap,
                     cluster=FALSE)
pdf("figs/CL_AP_region_haplogroup_genotypes_scaf94.pdf",width=12,height=8)
plot_grid(genos$positions,
          genos$genotypes,rel_heights = c(0.4,4),ncol=1,axis="tblr",align = "v")
dev.off()

# Haplogroups are all good

# These look ok so run with these...
out1<-data.frame(haplo1=haplogroup1)
write.table(out1,
            "outputs/haplogroup1_CL_AP_region.popmap",
            row.names = F,col.names = F,quote = F,sep="\t")
out3<-data.frame(haplo3=haplogroup3)
write.table(out3,
            "outputs/haplogroup3_CL_AP_region.popmap",
            row.names = F,col.names = F,quote = F,sep="\t")
