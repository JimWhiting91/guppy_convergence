# Permutation analysis of outlier window size relative to general distribution...
lib <- c("ggplot2","parallel","data.table")
sapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Get the outlier regions
outlier_info <- readRDS("outputs/full_outlier_info_for_all_rivers_wind10000.rds")

# Get the SNP positions for all SNPs in the VCF
vcf_snps <- data.frame(fread("data/five_aside_STAR.snps",header=F))
colnames(vcf_snps) <- c("chr","BP")
chrs <- unique(vcf_snps$chr)
chr_sizes <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")

# Convert to windows and return the number of SNPs
window_snp_counts <- data.frame(rbindlist(lapply(chrs,function(chr){
  print(paste0("Running ",chr))
  counts <- make_windows(data = vcf_snps[vcf_snps$chr == chr,],
                         minBP = 0,
                         maxBP = max(vcf_snps[vcf_snps$chr==chr,2]),
                         colname = "BP",
                         windsize = 10000,
                         avg="count")
  counts$chr <- chr
  return(counts)
})))

# Make ids
window_snp_counts$window_id <- paste0(window_snp_counts$chr,":",as.integer(window_snp_counts$BP1),"-",as.integer(window_snp_counts$BP2))
window_snp_counts <- na.omit(window_snp_counts)

# Now we need to run all of the permutations
outlier_info$river_measure <- paste0(outlier_info$river,"-",outlier_info$measure)
outlier_sets <- unique(outlier_info$river_measure)

# Outlier set SNP perms
perms=10000

# First analyse SNP counts of overlapping outlier windows
rivers <- unique(outlier_info$river)
river_outliers <- unlist(lapply(rivers,function(river){
  
  # Just get the river outliers
  outlier_counts <- table(outlier_info[outlier_info$river == river,"window_id"])
  outliers <- names(outlier_counts)[outlier_counts >= 2]
  return(outliers)
}))
overlapping_outliers <- table(river_outliers)
overlapping_outliers <- names(overlapping_outliers)[overlapping_outliers >= 2]

# Compare these to a permuted null
outlierN <- length(overlapping_outliers)

# Draw from data
median_SNP_counts <- unlist(mclapply(1:perms,function(perm){
  
  perm_windows <- window_snp_counts[sample(1:nrow(window_snp_counts),outlierN),"BP"]
  return(median(perm_windows))
},mc.cores = detectCores()-1))

# Get a p-value
observed_windows <- window_snp_counts[window_snp_counts$window_id %in% overlapping_outliers,]
observed_median <- median(observed_windows$BP)

# Compare to the distribution
observed_p <- length(median_SNP_counts[median_SNP_counts >= observed_median])/perms
observed_p

# Make a plot
plot_dd <- data.frame(count=median_SNP_counts)
plot_dd$type <- "Permuted"
plot_dd <- rbind(plot_dd,data.frame(count=observed_windows$BP,type="Observed"))
OOA_snp_count_histo <- ggplot(plot_dd,aes(x=count))+
  geom_histogram(binwidth = 1)+
 # geom_density()+
  geom_vline(xintercept = observed_median)+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text =  element_text(size=18),
        title = element_text(size=20))+
  labs(x="Median SNP Count",y="Count")+
  facet_wrap(~type,ncol=1,scales = "free_y")+
  ggtitle("Overlapping Selection Scan Outlier Windows")

# Also fetch the baypass outliers
baypass_outliers <- readRDS("outputs/baypass_AUXmodel_windows_10000_processed_Scaf94_merged.rds")
baypass_outliers <- baypass_outliers[baypass_outliers$outlier == "999",]

outlierN <- nrow(baypass_outliers)

# Draw from data
median_SNP_counts <- unlist(mclapply(1:perms,function(perm){
  
  perm_windows <- window_snp_counts[sample(1:nrow(window_snp_counts),outlierN),"BP"]
  return(median(perm_windows))
},mc.cores = detectCores()-1))

# Get a p-value
observed_windows <- window_snp_counts[window_snp_counts$window_id %in% baypass_outliers$window_id,]
observed_median <- median(observed_windows$BP)

# Compare to the distribution
observed_p <- length(median_SNP_counts[median_SNP_counts >= observed_median])/perms
observed_p

# Make a plot
plot_dd <- data.frame(count=median_SNP_counts)
plot_dd$type <- "Permuted"
plot_dd <- rbind(plot_dd,data.frame(count=observed_windows$BP,type="Observed"))
baypass_snp_count_histo <- ggplot(plot_dd,aes(x=count))+
  geom_histogram(binwidth = 1)+
  #geom_density()+
  geom_vline(xintercept = observed_median)+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text =  element_text(size=18),
        title = element_text(size=20))+
  labs(x="Median SNP Count",y="Count")+
  facet_wrap(~type,ncol=1,scales = "free_y")+
  ggtitle("Baypass Outlier Windows")

# Combine these together for a supp
pdf("figs/FigureSX_snp_count_effects.pdf",width=16,height=8)
cowplot::plot_grid(OOA_snp_count_histo,baypass_snp_count_histo,
                   ncol=2,nrow=1,align = "h",axis = "tblr",labels = "AUTO",label_size = 20)
dev.off()
