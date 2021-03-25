#####################################################################################
# Taking outlier regions, perform pathway analysis
# Load packages
lib<-c("UpSetR","biomaRt","clusterProfiler","ggplot2","data.table","reshape2","parallel","dplyr","tidyr")
lapply(lib,library,character.only=TRUE)
source("~/Exeter/five_aside/five_aside_functions.R")

# Fetch our outlier regions...
wind_size<-10000
outliers<-readRDS(paste0("outputs/AFD,XTX,XPEHH_overlapping_outliers_",wind_size,".rds"))[,1]

#######################################################################################
# Set up Ensembl
ensembl <- useMart("ensembl")
guppy <- useDataset("preticulata_gene_ensembl",mart=ensembl)
human <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
zebrafish <- useDataset("drerio_gene_ensembl",mart=ensembl)

filters <- listFilters(guppy)

guppy_universe<-getBM(attributes = c("ensembl_gene_id","entrezgene_id","kegg_enzyme"),
                                          mart=guppy)
guppy_entrez_universe<-as.character(unique(guppy_universe$entrezgene_id))

# Which Biomart attributes do we want?
biomart_atts<-c("external_gene_name","ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id","kegg_enzyme")

# Also set up orthlogous universes. These are human/zebrafish genes that have a one-to-one guppy orthologue
human_universe<-getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_type"),
                      mart=guppy)
human_universe<-human_universe[human_universe$hsapiens_homolog_orthology_type == "ortholog_one2one",]
human_entrez_universe<-getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                             filters = "ensembl_gene_id",
                             values = human_universe$hsapiens_homolog_ensembl_gene,
                             mart=human)

danio_universe<-getBM(attributes = c("ensembl_gene_id","drerio_homolog_ensembl_gene","drerio_homolog_orthology_type"),
                      mart=guppy)

danio_universe<-danio_universe[danio_universe$drerio_homolog_orthology_type == "ortholog_one2one",]
danio_entrez_universe<-getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                             filters = "ensembl_gene_id",
                             values = danio_universe$drerio_homolog_ensembl_gene,
                             mart=zebrafish)

# Write zebrafish orthologue ensembl IDs to a file
write.table(unique(danio_universe$drerio_homolog_ensembl_gene),
            "outputs/danio_orthologue_universe.txt",quote = F,row.names = F,col.names = F)

########################################################################################
# Now run for individual outlier sets...
river_outliers<-readRDS("outputs/outlier_candidate_regions_for_all_rivers.rds")
rivers<-names(river_outliers)

# Fetch ensmebl info
river_enrichment<-lapply(rivers,function(x){
  tmp<-ensembl_digger(outliers=river_outliers[[x]],
                       genome1="/Users/jimwhiting/Exeter/Genomes/STAR.chromosomes.release.fasta",
                       genome2="/Users/jimwhiting/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa",
                       aln_block_length = 5000,
                       aln_qual=60,
                       n.cores=6)
  
  enrich_ncbi <- enrichKEGG(gene=na.omit(unique(tmp$entrezgene_id)),
                            keyType = "kegg",
                            organism = 'pret',
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "fdr",
                            universe = guppy_entrez_universe)
  
  # Enrichment on zebrafish...
  danio_outliers<-danio_universe[danio_universe$ensembl_gene_id %in% na.omit(unique(tmp$ensembl_gene_id)),"drerio_homolog_ensembl_gene"]
  # danio_outliers<-danio_entrez_universe[danio_entrez_universe$ensembl_gene_id %in% danio_outliers,"entrezgene_id"]
  # danio_enrichment <- enrichKEGG(gene=na.omit(danio_outliers),
  #                                keyType = "kegg",
  #                                organism = 'dre',
  #                                minGSSize = 1,
  #                                pvalueCutoff = 0.05,
  #                                pAdjustMethod = "fdr",
  #                                universe = as.character(na.omit(danio_entrez_universe$entrezgene_id)))
  # 
  # # Enrichment on humans...
  # human_outliers<-human_universe[human_universe$ensembl_gene_id %in% na.omit(unique(tmp$ensembl_gene_id)),"hsapiens_homolog_ensembl_gene"]
  # human_outliers<-human_entrez_universe[human_entrez_universe$ensembl_gene_id %in% human_outliers,"entrezgene_id"]
  # human_enrichment <- enrichKEGG(gene=na.omit(human_outliers),
  #                                keyType = "kegg",
  #                                organism = 'hsa',
  #                                minGSSize = 1,
  #                                pvalueCutoff = 0.05,
  #                                pAdjustMethod = "fdr",
  #                                universe = as.character(na.omit(human_entrez_universe$entrezgene_id)))
  
  return(list(tmp,enrich_ncbi,danio_outliers))
})

names(river_enrichment)<-rivers

# Save Danio orthologues to a file
lapply(rivers,function(x){
  write.table(river_enrichment[[x]][[3]],
              paste0("outputs/river_specific_outlier_danio_ensembl_",x,".txt"),
              quote = F,col.names = F,row.names = F)
})

# And do keggMENRICH
module_enrichment<-lapply(1:5,function(x){
  enrichMKEGG(gene=na.omit(unique(river_enrichment[[x]][[1]]$entrezgene_id)),
           keyType = "kegg",
           organism = 'pret',
           pvalueCutoff = 0.05,
           pAdjustMethod = "fdr")
})

########################################################################################
# How about scaffold 94 and chromosome 20?
caroni_region<-c("chr20:1-3164071","000094F_0:1-1800000")


caroni_ensembl<-ensembl_digger(outliers=caroni_region,
                               genome1="/Users/jimwhiting/Exeter/Genomes/STAR.chromosomes.release.fasta",
                               genome2="/Users/jimwhiting/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa",
                               aln_block_length = 1000,
                               aln_qual=60,
                               biomart_attributes=biomart_atts,
                               n.cores=6)

caroni_enrich<-enrichKEGG(gene=na.omit(unique(caroni_ensembl$entrezgene_id)),
                          keyType = 'kegg',
                          organism = 'pret',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "fdr",
                          universe = guppy_entrez_universe)

# Enrichment on zebrafish...
danio_outliers<-danio_universe[danio_universe$ensembl_gene_id %in% na.omit(unique(caroni_ensembl$ensembl_gene_id)),"drerio_homolog_ensembl_gene"]
danio_outliers<-danio_entrez_universe[danio_entrez_universe$ensembl_gene_id %in% danio_outliers,"entrezgene_id"]
danio_enrichment <- enrichKEGG(gene=na.omit(danio_outliers),
                               keyType = "kegg",
                               organism = 'dre',
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "fdr",
                               universe = as.character(na.omit(danio_entrez_universe$entrezgene_id)))

# Enrichment on humans...
human_outliers<-human_universe[human_universe$ensembl_gene_id %in% na.omit(unique(caroni_ensembl$ensembl_gene_id)),"hsapiens_homolog_ensembl_gene"]
human_outliers<-human_entrez_universe[human_entrez_universe$ensembl_gene_id %in% human_outliers,"entrezgene_id"]
human_enrichment <- enrichKEGG(gene=na.omit(human_outliers),
                               keyType = "kegg",
                               organism = 'hsa',
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "fdr",
                               universe = as.character(na.omit(human_entrez_universe$entrezgene_id)))

##########################################
# Panther Enrichment Analysis
panther_lists<-lapply(rivers,function(x){
  tmp<-read.table(paste0("outputs/pather_results_",x,"_outliers.tsv"),sep = "\t",header=T)
  colnames(tmp)<-c("pathway","universe_N","outlier_N","outlier_Exp","direction","fold_enrich","p","fdr")
  return(tmp)
})
names(panther_lists)<-rivers

# First, which pathways have more gene than expected in every river?
all_river_paths<-lapply(panther_lists,function(x){
  return(x[x$fold_enrich > 1,"pathway"])
})
counts<-table(unlist(all_river_paths))
observed_paths<-names(counts[counts > 3])
observed_paths<-observed_paths[observed_paths != "Unclassified (UNCLASSIFIED)"]
  
# Permutation tests
river_gene_counts<-sapply(panther_lists,function(x){return(sum(x$outlier_N))})

# Set up the universe to draw from
universe<-data.frame(panther_lists[["Tacarigua"]][,1:2])
universe$prob<-universe$universe_N/sum(universe$universe_N)
universe2<-unlist(sapply(1:nrow(universe),function(x){return(rep(universe$pathway[x],universe$universe_N[x]))}))

# Run the permutations
iterations<-10000

if(!(file.exists(paste0("outputs/panther_enrichment_permutations_",iterations,".rds")))){
paths_perm<-mclapply(1:iterations,function(iter){
  
  # For each river, draw an equivalent set of path representatives
  tmp_rivers<-lapply(rivers,function(x){
    
    # Draw lots
    random_draw<-sample(x=universe2,
                        size=river_gene_counts[x],
                        replace = F)
    
    # Count them
    random_counts<-data.frame(table(random_draw))
    
    # Expecteds
    universe_expected<-universe
    universe_expected$exp<-universe_expected$prob*river_gene_counts[x]
    
    # Get expected
    for(i in 1:nrow(random_counts)){
    random_counts$exp[i]<-universe_expected[universe_expected$pathway == random_counts$random_draw[i],"exp"]
    }
    
    # Only keep those with greater than expected
    to_keep<-random_counts[random_counts$Freq > random_counts$exp,"random_draw"]
    
    # Keep only uniques
    return(as.character(to_keep))
  })
  
  # Return those which appear in all 5 sets
  tmp_counts<-table(unlist(tmp_rivers))
  tmp_full<-names(tmp_counts[tmp_counts == 5])
  
},mc.cores=detectCores()-1)

# Save the permutations
saveRDS(paths_perm,paste0("outputs/panther_enrichment_permutations_",iterations,".rds"))

}

# Read in perms
paths_perm<-readRDS(paste0("outputs/panther_enrichment_permutations_",iterations,".rds"))

# How many permutations contain each of our observed pathways?
observed_perms<-data.frame(path=observed_paths)
for(i in 1:nrow(observed_perms)){
  observed_perms$perm_count[i]<-table(unlist(paths_perm))[observed_perms$path[i]]
}

observed_perms$p<-observed_perms$perm_count/iterations
observed_perms$p[is.na(observed_perms$p)]<-0
observed_perms$fdr<-p.adjust(observed_perms$p,method = "fdr")

# And just the significant ones...
significant_perms<-observed_perms[observed_perms$fdr < 0.05,]

# Plot the fold-enrichment for significant pathways...
to_plot<-data.frame(rbindlist(lapply(significant_perms$path,function(path){
  
  tmp<-rbindlist(lapply(panther_lists,function(x){
    tmp2<-x[x$pathway == path,]
  }))
  
  # Add rivers
  tmp$river<-rivers
return(tmp)
})))

# Split
to_plot<-separate(to_plot,"pathway",into = c("path1","path2"),sep = "pathway ")

# Tidy
to_plot$river_F<-factor(to_plot$river,levels=five_aside_colour_rivers$river)
to_plot$path1<-gsub("pathway","pathway\n",to_plot$path1)

# Plot the bars with asterisk
to_plot$fold_enrich<-as.numeric(to_plot$fold_enrich)
to_plot$path1_F <- factor(to_plot$path1,levels=rev(unique(to_plot$path1)))

label_dd <- data.frame(path1_F = to_plot$path1[1],
                       fold_enrich = 3.5,
                       label = "*",
                       river_F = "Aripo")

paths_plot<-ggplot(to_plot,aes(y=path1_F,x=fold_enrich,fill=river_F))+
  geom_bar(position="dodge",stat="identity")+
  geom_vline(xintercept = 1)+
  theme_bw()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(fill="River",x="Fold Enrichment")+
  scale_fill_manual(breaks=five_aside_colour_rivers$river_F,
                    values = five_aside_colour_rivers$colour)+
  geom_text(data=label_dd,aes(y=path1_F,x=fold_enrich),label=label_dd$label,size=16)

# Write to a file
saveRDS(paths_plot,"outputs/pathway_plot_figure.rds")


# Retrieve Cadherin-signalling genes --------------------------------------
cadherin_T <- list(c("ENSDARG00000055825","ENSDARG00000036175","ENSDARG00000024785","ENSDARG00000078088"),
                   c("Cadherin, EGF LAG seven-pass G-type receptor 3","Protocadherin 1b", "Catenin alpha-2","Si:ch211-186j3.6"))
cadherin_G <- list(c("ENSDARG00000007561","ENSDARG00000024785"),
                   c("Cadherin-23","Catenin alpha-2"))
cadherin_A <- list(c("ENSDARG00000078226","ENSDARG00000098652","ENSDARG00000036175","ENSDARG00000087709"),
                   c("Cadherin-12","Protocadherin-11 X-linked","Protocadherin-1","Protocadherin Fat 3"))
cadherin_O <- list(c("ENSDARG00000015002","ENSDARG00000060610","ENSDARG00000079850","ENSDARG00000007561","ENSDARG00000017803","ENSDARG00000103903","ENSDARG00000071107","ENSDARG00000018923"),
                   c("Cadherin-4","Protocadherin-7b","Daschsous cadherin-related 1b","Cadherin-23","Glycogen synthase kinase 3","Uncharacterised","Protein Wnt","FAT atypical cadherin 2"))
cadherin_M <- list(c("ENSDARG00000055825","ENSDARG00000015002","ENSDARG00000098824","ENSDARG00000007561","ENSDARG00000041117","ENSDARG00000102185","ENSDARG00000077996","ENSDARG00000078088","ENSDARG00000103903","ENSDARG00000031894"),
                  c("Cadherin, EGF LAG seven-pass G-type receptor 3","Cadherin-4","Cadherin-22","Cadherin-23","Protein Wnt-2","Protocadherin-8","Cadherin-24","Si:ch211-186j3.6","Cadherin-10","Lymphocyte enhancer binding factor 1"))

# Visualise
danio_list<-list(cadherin_T[[1]],cadherin_G[[1]],cadherin_A[[1]],cadherin_O[[1]],cadherin_M[[1]])
names(danio_list)<-rivers
upset(fromList(danio_list))

# Return to Guppy Cadherins with old genome co-ordinates...
guppy_cadherins <- mclapply(rivers,function(x){
  
  # Fetch from biomart
bm_out <- getBM(attributes = c("preticulata_homolog_ensembl_gene"),
          filters = "ensembl_gene_id",
          values = danio_list[[x]],
          mart=zebrafish)

  # Fetch common names
bm_out2 <- getBM(attributes = c("external_gene_name","description","ensembl_gene_id","chromosome_name","start_position","end_position"),
                filters = "ensembl_gene_id",
                values = bm_out[,1],
                mart=guppy)
bm_out2$river <- x
return(bm_out2)
},mc.cores=5)

# Now need to reverse transform each gene back to STAR co-ordinates...
cadherin_dd<-data.frame(rbindlist(guppy_cadherins))
for(i in 1:nrow(cadherin_dd)){
  chr <- cadherin_dd$chromosome_name[i]
  start <- cadherin_dd$start_position[i]
  end <- cadherin_dd$end_position[i]

  # Align
  system(paste0("samtools faidx ~/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa ",chr,":",start,"-",end," > outputs/temp.fa"))
  system(paste0("~/software/minimap2/minimap2 -t 6 ~/Exeter/Genomes/STAR.chromosomes.release.fasta outputs/temp.fa > outputs/temp.paf"))
  
  # Fetch alignment
  aln <- read.table("outputs/temp.paf",fill=T)
  
  # Keep good alignments...
  aln <- aln[aln$V12 == 60 &
               aln$V10 > 1000,]
  
  # Remove bum chr
  chr_to_keep <- names(table(aln$V6)[table(aln$V6) == max(table(aln$V6))])
  aln <- aln[aln$V6 == chr_to_keep,]
  
  # Any rearrangements?
  for(j in 1:nrow(aln)){
    if(aln$V9[j] < aln$V8[j]){
      save_this <- aln$V8[j]
      aln$V8[j] <- aln$V9[j]
      aln$V9[j] <- save_this
    }
  }
  
  # Now keep start and end...
  cadherin_dd$chromosome_name[i] <- chr_to_keep
  cadherin_dd$start_position[i] <- min(as.integer(aln$V8))
  cadherin_dd$end_position[i] <- max(as.integer(aln$V9))
  
}

# Save the results
write.table(cadherin_dd,
        "tables/cadherin_pathway_genes_by_river.txt",
            quote = F,row.names = F,sep = "\t")

