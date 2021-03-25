library("data.table")

snps<-data.frame(fread("data/five_aside_STAR_snps"))

snp_max<-round(nrow(snps)/16)

genos<-data.frame(fread("data/five_aside_STAR.geno"))

# Sample every snp_interval
for (i in 1:16){

rows<-seq(i,nrow(snps),16)

datasub<-genos[rows,]
datasub2<-snps[rows,]

write.table(datasub,
paste0("/gpfs/ts0/home/jw962/HP_LP_trials/baypass/data/subsets/five_aside_STAR_sub",i,".geno"),
row.names=F,quote=F,sep="\t")

write.table(datasub2,
paste0("/gpfs/ts0/home/jw962/HP_LP_trials/baypass/data/subsets/five_aside_STAR_sub",i,".snps"),
row.names=F,quote=F,sep="\t")

}
