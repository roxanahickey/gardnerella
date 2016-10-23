## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information
## for downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## Last update: 2015-04-22

# This script summarizes OrthoMCL gene clusters and outputs FASTA IDs and gene annotations for each cluster and singleton coding sequence (CDS).
# It also performs hierarchical clustering of genomes based on the gene clusters and assigns genomes to groups or clades. It then extracts FASTA
# IDs and annotations for each clade. Finally, it produces summary figures of core and pan genome CDS counts. I have written this script to be
# sourced from another script that inputs the necessary data (e.g. run_orthomcl_summary_gv.R). That script should first be edited to specify the
# output directory, OrthoMCL cluster data, keys for FASTA IDs, gene features and annotations (in my case, from PATRIC CDS sequences), and custom
# color palettes.

library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(reshape)
source("code/ggplot_cust.R")

##################################################
## create directories for output data
dir.create(dir.out)
dir.create(dir.fig)

dir.create(paste(dir.out, "/clust-fasta-ids", sep=""))
dir.create(paste(dir.out, "/single-fasta-ids", sep=""))
dir.create(paste(dir.out, "/clade-info", sep=""))

dir.create(paste(dir.out, "/clust-anno", sep=""))
dir.create(paste(dir.out, "/single-anno", sep=""))
dir.create(paste(dir.out, "/clade-anno", sep=""))

for(i in c("ec", "go", "feat", "figfam", "path")){
  dir.create(paste(dir.out, "/clust-anno/", i, sep=""))
  dir.create(paste(dir.out, "/single-anno/", i, sep=""))
  dir.create(paste(dir.out, "/clade-anno/", i, sep=""))
}

##################################################
## get info from OrthoMCL clusters
ortho.clust <- list()
ortho.clust$orig <- ortho.clust.in

# list of clusters and fasta IDs per cluster
ortho.clust$ortho_id <- unlist(strsplit(ortho.clust$orig, split=": "))[seq(2, length(ortho.clust$orig)*2, by=2)]
names(ortho.clust$ortho_id) <- unlist(strsplit(ortho.clust$orig, split=": "))[seq(1, length(ortho.clust$orig)*2, by=2)]

## key for orthoMCL cluster IDs to CDS IDs
key.ortho <- data.frame(ortho_id=substr(names(unlist(strsplit(ortho.clust$ortho_id, " "))), start=1, stop=10),
                        cds_id=unlist(strsplit(ortho.clust$ortho_id, " "), use.names=F))

## merge ortho IDs with PATRIC features key
key.feat.cds <- merge(key.feat.cds, key.ortho, all.x=T)

# write function to extract annotation information by matching CDS IDs
get_anno <- function(x) {anno.tmp <- x[x$cds_id %in% unlist(cds.tmp),]
                         anno.tmp <- merge(anno.tmp, key.ortho, all.x=T)
                         anno.type <- unique(x$anno_type)
                         if(nrow(anno.tmp) > 0) {write.table(anno.tmp, 
                                                             paste(dir.out.anno, "/", anno.type, "/", i, "_", anno.type, ".txt", sep=""), 
                                                             sep="\t", row.names=F, quote=F)}}

# write fasta ID lists and annotation tables for each OrthoMCL cluster
for(i in unique(names(ortho.clust$ortho_id))){
  # extract CDS IDs for cluster
  cds.tmp <- strsplit(ortho.clust$ortho_id[i], " ")
  write.table(unlist(cds.tmp, use.names=F), paste(dir.out, "/clust-fasta-ids/", i, "_fasta_ids.txt", sep=""), row.names=F, col.names=F, quote=F)
  
  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/clust-anno", sep="")
  
  # extract PATRIC features for CDS IDs in cluster
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}

  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in cluster
  sapply(key.anno, get_anno)
}

# clean up
rm(i, dir.out.anno, cds.tmp, feat.tmp)

##################################################
## get info from OrthoMCL singletons
ortho.sing <- list()

# list of singleton fasta IDs
ortho.sing$cds_id <- ortho.sing.in

# strain names
ortho.sing$id_short <- unlist(strsplit(ortho.sing$cds_id, split="|", fixed=T))[seq(1, length(ortho.sing$cds_id)*2, by=2)]

# initialize empty data frame to keep track of counts of singleton clusters per strain (more info to be added later)
strain.uniq.df <- data.frame(matrix(0, nrow=length(unique(ortho.sing$id_short)), ncol=3))
rownames(strain.uniq.df) <- sort(unique(ortho.sing$id_short))
colnames(strain.uniq.df) <- c("uniq_clust_count", "total_clust_count", "id_abbr")

# write singleton fasta IDs and annotation tables per strain
for(i in unique(ortho.sing$id_short)){
  # extract CDS IDs for strain singletons
  cds.tmp <- ortho.sing$cds_id[grep(i, ortho.sing$cds_id)]
  write.table(cds.tmp, 
              paste(dir.out, "/single-fasta-ids/", i, "_uniq_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # count number of singletons per strain
  strain.uniq.df[i, "uniq_clust_count"] <- length(cds.tmp)

  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/single-anno", sep="")

  # extract PATRIC features for CDS IDs in cluster
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}

  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in cluster
  sapply(key.anno, get_anno)
}

# clean up
rm(i, cds.tmp, dir.out.anno, feat.tmp)

##################################################
## count clusters using presence/absence matrix
clust.mx <- clust.mx.in

# count occurrences in each cluster
clust.ct <- as.data.frame(table(colSums(clust.mx)))
colnames(clust.ct) <- c("no_genomes", "freq")

# add singletons to clust.ct
clust.ct$freq[1] <- clust.ct$freq[1] + length(ortho.sing$cds_id)

##################################################
## get clust IDs of core gene clusters (+ "relaxed" core)
core.clust <- list()
core.clust$all <- colnames(clust.mx)[colSums(clust.mx)==nrow(clust.mx)]
core.clust$less1 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 1]
core.clust$less2 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 2]

write.table(core.clust$all,
            file=paste(dir.out, "/clade-info/core_clust_all.txt", sep=""), 
            row.names=F, col.names=F, quote=F)

write.table(core.clust$less1,
            file=paste(dir.out, "/clade-info/core_clust_less1.txt", sep=""), 
            row.names=F, col.names=F, quote=F)

write.table(core.clust$less2,
            file=paste(dir.out, "/clade-info/core_clust_less2.txt", sep=""), 
            row.names=F, col.names=F, quote=F)

##################################################
## rename original clust.mx -- binary presence/absence
clust.mx.bin <- clust.mx

## convert clust.mx to show frequency of each gene cluster in each genome
clust.mx.freq <- matrix(nrow=nrow(clust.mx), ncol=ncol(clust.mx), dimnames=list(rownames(clust.mx), colnames(clust.mx)))

for(i in unique(colnames(clust.mx))){
  key.tmp <- subset(key.ortho, ortho_id==i)
  for(j in unique(rownames(clust.mx))){
    clust.mx.freq[j, i] <- length(grep(j, key.tmp$cds_id))
  }
}

clust.mx <- clust.mx.freq

rm(i, j, key.tmp, clust.mx.freq)

##################################################
## run hclust to identify groups
source("code/orthomcl_hclust_id_clades.R")

# add ortho clade id to genome.key
genome.key$ortho_clade <- cutg[genome.key$id_short]

##################################################
## clade core genes
# intialize empty data frame to fill as you calculate clusters
clust.smy <- data.frame(matrix(0, nrow=length(unique(cutg)), ncol=7))
rownames(clust.smy) <- paste("clade", unique(cutg), sep="")
colnames(clust.smy) <- c("all_core", "clade_core_not_uniq", "clade_core_uniq",
                         "clade_core_red", "clade_pan", "clade_sing", "clade_total")

# initialize empty lists to fill in
clade.id <- list()
clade.list <- list()

for(i in unique(cutg)){
  # list clade member ids
  clade.tmp <- names(cutg[cutg==i])
  j <- paste("clade", i, sep="")
  clade.id[j] <- list(clade.tmp)
  clade.list[i] <- paste("clade", i, sep="")
  
  # write list of strain IDs in cluster
  write.table(clade.tmp, 
              file=paste(dir.out, "/clade-info/clade", i, "_id.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

rm(i, j, clade.tmp)

clade.list <- unlist(clade.list, use.names=F)

clade.clust <- list()
clade.core.clust <- list()
clade.core.clust.red <- list()

# write new get_anno function to work with clades -- includes j term to specify cluster summarized
get_anno_clade <- function(x) {anno.tmp <- x[x$cds_id %in% unlist(cds.tmp, use.names=F),]
                               anno.tmp <- merge(anno.tmp, key.ortho, all.x=T)
                               anno.type <- unique(x$anno_type)
                               if(nrow(anno.tmp) > 0) {write.table(anno.tmp, 
                                                                   paste(dir.out.anno, "/", anno.type, "/", i, "_", j, "_", anno.type, ".txt", sep=""), 
                                                                   sep="\t", row.names=F, quote=F)}}

for(i in unique(clade.list)){

  id.tmp <- unlist(clade.id[i], use.names=F)
  
  ## ALL clusters in clade
  # reduce clust.mx.bin to clade members
  clade.clust.tmp <- clust.mx.bin[id.tmp,]
  clade.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)>=1]  
  
  # append list to clade.clust
  clade.clust[i] <- list(colnames(clade.clust.tmp))
  
  # write list of cluster IDs in clade
  write.table(colnames(clade.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/clade-anno", sep="")
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(clade.clust.tmp)], " ")
  j <- "clust"
  
  # extract PATRIC features for CDS IDs in clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in clusters
  sapply(key.anno, get_anno_clade)
  
  ## CORE clusters in clade
  # reduce clade.clust.tmp to core clusters
  clade.core.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)==nrow(clade.clust.tmp)]

  # append to clade.core.clust
  clade.core.clust[i] <- list(colnames(clade.core.clust.tmp))

  # write list of cluster IDs in core set in clade
  write.table(colnames(clade.core.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(clade.core.clust.tmp)], " ")
  j <- "clust_core"
  
  # extract PATRIC features for CDS IDs in core clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in core clusters
  sapply(key.anno, get_anno_clade)
    
  ## REDUCED CORE clusters in clade (i.e. not also present in core for all strains)
  # reduce clade core to those not in overall core
  clade.core.clust.red.tmp <- clade.core.clust.tmp[,!(colnames(clade.core.clust.tmp) %in% core.clust$all)]

  # append to clade.core.clust.red
  clade.core.clust.red[i] <- list(colnames(clade.core.clust.red.tmp))
  
  # write list of cluster IDs in reduced core set in clade
  write.table(colnames(clade.core.clust.red.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core_red.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(clade.core.clust.red.tmp)], " ")
  j <- "clust_core_red"
  
  # extract PATRIC features for CDS IDs in reduced core clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_red_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in reduced core clusters
  sapply(key.anno, get_anno_clade)
    
  # fill in summary table
  clust.smy[i, "all_core"] <- ncol(clade.core.clust.tmp) - ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_core_red"] <- ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_sing"] <- length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  clust.smy[i, "clade_total"] <- ncol(clade.clust.tmp) + length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  
}

rm(i, j, id.tmp, clade.clust.tmp, clade.core.clust.tmp, clade.core.clust.red.tmp, dir.out.anno, cds.tmp, feat.tmp)

clust.smy$clade_pan <- clust.smy$clade_total - (clust.smy$all_core + clust.smy$clade_core_red + clust.smy$clade_sing)

##################################################
## clade core unique genes
# initialize empty list
clade.core.clust.uniq <- list()

for(i in unique(clade.list)){
  # select focal clade, ignore others
  ignore <- clade.list[!(clade.list %in% i)]
  ignore.clust <- unique(unlist(clade.clust[ignore], use.names = F))
  query.clust <- unlist(clade.core.clust.red[i], use.names = F)
  
  # append unique core cluster IDs to clade.core.clust.uniq
  clade.core.clust.uniq[i] <- list(query.clust[!(query.clust %in% ignore.clust)])
  
  # count length of core unique cluster IDs, append to summary table
  clust.smy[i, "clade_core_uniq"] <- length(query.clust[!(query.clust %in% ignore.clust)])
  
  # write list of core unique cluster IDs in clade
  write.table(query.clust[!(query.clust %in% ignore.clust)], 
              file=paste(dir.out, "/clade-info/", i, "_clust_core_uniq.txt", sep=""),
              row.names=F, col.names=F, quote=F)
  
  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/clade-anno", sep="")
  
  # extract CDS IDs from core unique clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[query.clust[!(query.clust %in% ignore.clust)]], " ")
  j <- "clust_core_uniq"
  
  # extract PATRIC features for CDS IDs in clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_uniq_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in clusters
  sapply(key.anno, get_anno_clade)
}

rm(i, j, ignore, ignore.clust, query.clust, dir.out.anno, cds.tmp, feat.tmp)

# finish filling out summary table
clust.smy$clade_core_not_uniq <- clust.smy$clade_core_red - clust.smy$clade_core_uniq

# reduce summary table to make nonredundant
clust.smy.red <- clust.smy[,-c(4,7)]
clust.smy.red$clade_id <- rownames(clust.smy.red)
clust.smy.red.md <- melt(clust.smy.red, id="clade_id")

write.table(clust.smy.red, file=paste(dir.out, "/clade-info/clust_count_summary.txt", sep=""), sep="\t", quote=F)

##################################################
## singletons per genome -- complete data frame started earlier
strain.uniq.df$total_clust_count <- rowSums(clust.mx.bin)[match(rownames(strain.uniq.df), names(rowSums(clust.mx.bin)))] + strain.uniq.df$uniq_clust_count

strain.uniq.df$uniq_clust_prop <- strain.uniq.df$uniq_clust_count / strain.uniq.df$total_clust_count

strain.uniq.df$id_abbr <- genome.key$id_abbr[match(rownames(strain.uniq.df), genome.key$id_short)]

strain.uniq.df$ortho_clade <- genome.key$ortho_clade[match(rownames(strain.uniq.df), genome.key$id_short)]

##################################################
## plots!
# number of gene clusters present in exactly n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=freq, fill=no_genomes)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("No. genomes") + ylab("Gene cluster frequency") +
  geom_text(aes(label=freq, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, col.bw.seq))(length(clust.ct$no_genomes))) +
  theme(legend.position="none")
ggsave(paste(dir.fig, "/orthomcl_clusters_n_genome_exact.png", sep=""), width=6, height=4, units="in", scale=2)

# number of gene clusters present in at least n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=rev(cumsum(rev(clust.ct$freq))), fill=no_genomes)) + 
  geom_bar(stat="identity") + 
  xlab("No. genomes") + ylab("Gene cluster frequency") + theme_cust + 
  geom_text(aes(label=rev(cumsum(rev(clust.ct$freq))), angle=45, hjust=0, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, col.bw.seq))(length(clust.ct$no_genomes))) +
  theme(legend.position="none")
ggsave(paste(dir.fig, "/orthomcl_clusters_n_genome_atleast.png", sep=""), width=6, height=4, units="in", scale=2)

# cluster counts per clade
ggplot(clust.smy.red.md, aes(x=clade_id, y=value, fill=variable)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("Clade") + ylab("Gene cluster frequency") +
  labs(fill="") +
  scale_fill_manual(limits=c("all_core", "clade_core_not_uniq", "clade_core_uniq", "clade_pan", "clade_sing"),
                    labels=c("All core (all genomes, all clades)", "Clade core, non-unique (all genomes in clade)", 
                             "Clade core, unique (all genomes in clade)", "Clade accessory (2+ genomes in clade)", 
                             "Clade singletons (1 genome in clade)"), 
                    values=rev(brewer.pal(5, col.bw.seq)))
ggsave(paste(dir.fig, "/orthomcl_clusters_by_clade.png", sep=""), width=6, height=4, units="in", scale=2)

# number of singleton clusters per genome
ggplot(strain.uniq.df, aes(x=reorder(id_abbr, ortho_clade - uniq_clust_prop), y=uniq_clust_prop, fill=factor(ortho_clade))) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("Genome") + ylab("Proportion singleton gene clusters out of total per genome") + 
  theme(axis.text.x = element_text(angle=-30, hjust=0, vjust=1, size=6)) +
  labs(fill="Clade") +
  geom_text(aes(label=uniq_clust_count, vjust=1.5), size=2.5, color="white", show_guide=F) +
  scale_fill_manual(values=col.cust.hc)
ggsave(paste(dir.fig, "/orthomcl_clusters_uniq_per_strain.png", sep=""), width=6, height=4, units="in", scale=2)

##################################################
## clean up, save image
rm(ortho.clust.in, ortho.sing.in, clust.mx.in)
save.image(paste(dir.out, "/orthomcl_summary.RData", sep=""))
