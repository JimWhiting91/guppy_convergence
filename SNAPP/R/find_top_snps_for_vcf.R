library(vcfR)
library(hierfstat)
library(adegenet)

vcf = read.vcfR("~/Exeter/tmp/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2.nomiss.thin50k.recode.vcf")

# Get pop strata
inds <- colnames(vcf@gt)[2:ncol(vcf@gt)]

pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
river_labs <- c("TAC","TAC","G","G","AP","AP","OROP","OROP","MAD","MAD")
names(river_labs) <- pops

metadata <- data.frame(ind=inds,
                       pop="WINGEI",
                       pop2="WINGEI")
for(pop in pops){
  metadata[grep(pop,metadata$ind),"pop2"] <- river_labs[pop]
  metadata[grep(pop,metadata$ind),"pop"] <- pop
}

# First get the pred level...

# Get Fst by converting and calculating
genind_dd <- vcfR2genind(vcf)
pop(genind_dd) <- metadata$pop
dat2 <- genind2hierfstat(genind_dd)
f_stats <- basic.stats(dat2)

# And keep top 1000 SNPs based on perloc Fst and that are on high quality autosomes...
f_stats_chr <- f_stats$perloc
f_stats_chr <- f_stats_chr[grep("chr",rownames(f_stats_chr)),]
f_stats_chr <- f_stats_chr[grep("chr12",rownames(f_stats_chr),invert = T),]
top_snps1000 <- rownames(head(f_stats_chr[order(-f_stats_chr$Fst),],1000))
out1000 <- data.frame(snp=rownames(f_stats_chr)[rownames(f_stats_chr) %in% top_snps1000])

write.table(out1000,"~/Exeter/tmp/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2_nochr12_highFST_HPLP.txt",
            row.names = F,quote = F,sep = "\t",col.names=F)

# Then the river level...

# Get Fst by converting and calculating
pop(genind_dd) <- metadata$pop2
dat2 <- genind2hierfstat(genind_dd)
f_stats <- basic.stats(dat2)

# And keep top 1000 SNPs based on perloc Fst and that are on high quality autosomes...
f_stats_chr <- f_stats$perloc
f_stats_chr <- f_stats_chr[grep("chr",rownames(f_stats_chr)),]
f_stats_chr <- f_stats_chr[grep("chr12",rownames(f_stats_chr),invert = T),]
top_snps1000 <- rownames(head(f_stats_chr[order(-f_stats_chr$Fst),],1000))
out1000 <- data.frame(snp=rownames(f_stats_chr)[rownames(f_stats_chr) %in% top_snps1000])

write.table(out1000,"~/Exeter/tmp/pic_wing_five_aside_STAR_v3_RIVERS_WINGEI_SNAPP_v2_nochr12_highFST_RIVER.txt",
            row.names = F,quote = F,sep = "\t",col.names=F)

