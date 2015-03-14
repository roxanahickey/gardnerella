## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information
## for downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-13

setwd("~/Documents/research/gardnerella/")
library(ggplot2)
source("scripts/ggplot_cust.R")

## create directories for output data
dir.create("output/orthomcl-summary")
dir.create("output/RData")

## read in ID key containing short and full strain names
gv.id <- read.table("data/gv_id_key.txt", sep="\t", header=T)

##################################################
## read in OrthoMCL clusters
ortho.clust.orig <- readLines("data/orthomcl_gv_groups_70.txt")

# list of clusters and fasta IDs per cluster
ortho.clust <- unlist(strsplit(ortho.clust.orig, split=": "))[seq(2, length(ortho.clust.orig)*2, by=2)]
names(ortho.clust) <- unlist(strsplit(ortho.clust.orig, split=": "))[seq(1, length(ortho.clust.orig)*2, by=2)]

dir.create("output/orthomcl-summary/clust-fasta-ids")

# write lists of fasta IDs per cluster
for(i in unique(names(ortho.clust))){
  clust.tmp <- data.frame(strsplit(ortho.clust[i], " "))
  write.table(clust.tmp, paste("output/orthomcl-summary/clust-fasta-ids/", i, "_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

# clean up
rm(i, clust.tmp)

##################################################
## read in OrthoMCL singletons
ortho.sing.orig <- readLines("data/orthomcl_gv_singletons_70.txt")

# list of singleton fasta IDs
ortho.sing <- unlist(strsplit(ortho.sing.orig, split=": "))[seq(2, length(ortho.sing.orig)*2, by=2)]

# "cluster" names of singleton fasta IDs
ortho.sing.id <- unlist(strsplit(ortho.sing.orig, split=": "))[seq(1, length(ortho.sing.orig)*2, by=2)]

# GV strain names
ortho.sing.strain <- unlist(strsplit(ortho.sing, split="|", fixed=T))[seq(1, length(ortho.sing)*2, by=2)]

dir.create("output/orthomcl-summary/sing-fasta-ids")

# write lists of singleton fasta IDs per strain
for(i in unique(ortho.sing.strain)){
  sing.tmp <- ortho.sing[grep(i, ortho.sing)]
#   assign(paste(i, ".uniq", sep=""), sing.tmp)
  write.table(sing.tmp, paste("output/orthomcl-summary/sing-fasta-ids/", i, "_uniq_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  }

# clean up
rm(i, sing.tmp)

##################################################
## read in cluster presence-absence matrix (*only* values of 0 or 1 for absence/presence)
clust.mx <- read.table("data/orthomcl_gv_groups_70_mtx.txt", header=T)

# replace "single" colnames with singleton "cluster" IDs
colnames(clust.mx)[grep("single", colnames(clust.mx))] <- ortho.sing.id

# count occurrences in each cluster
clust.ct <- as.data.frame(table(colSums(clust.mx)))
colnames(clust.ct) <- c("no_genomes", "freq")

# plot number of gene clusters present in exactly n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=freq)) + 
  geom_bar(stat="identity") + 
  xlab("No. genomes") + ylab("No. gene clusters") + theme_cust +
  geom_text(aes(label=freq, vjust=-0.5, size=0.9)) +
  ggtitle("No. gene clusters in exactly n genomes") +
  theme(legend.position="none")
ggsave("figs/orthomcl_clusters_n_genome_exact.png", width=6, height=4, units="in", scale=2)

# plot number of gene clusters present in at least n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=rev(cumsum(rev(clust.ct$freq))))) + 
  geom_bar(stat="identity") + 
  xlab("No. genomes") + ylab("No. gene clusters") + theme_cust + 
  geom_text(aes(label=rev(cumsum(rev(clust.ct$freq))), vjust=-0.5, size=0.9)) +
  ggtitle("No. gene clusters in at least x genomes") +
  theme(legend.position="none")
ggsave("figs/orthomcl_clusters_n_genome_atleast.png", width=6, height=4, units="in", scale=2)

##################################################
# get clust IDs of core gene clusters (+ "relaxed" core)
clust.core.all <- colnames(clust.mx)[colSums(clust.mx)==nrow(clust.mx)]
clust.core.less1 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 1]
clust.core.less2 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 2]

##################################################
# run hclust to identify groups
source("scripts/orthomcl_hclust_id_clades.R")

##################################################
## clade-unique core
dir.create("output/orthomcl-summary/clade-info")

for(i in unique(cutg)){
  # list clade member ids
  clade.tmp <- names(cutg[cutg==i])
  assign(paste("clade", i, ".id", sep=""), clade.tmp)
  write.table(clade.tmp, 
              file=paste("output/orthomcl-summary/clade-info/clade", i, "_id.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # reduce clust to clade members
  clade.clust.tmp <- clust.mx[clade.tmp,]
  clade.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)>=1]
  assign(paste("clade", i, ".clust", sep=""), clade.clust.tmp)
  write.table(colnames(clade.clust.tmp), 
              file=paste("output/orthomcl-summary/clade-info/clade", i, "_clust.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract clade core clusters
  clade.clust.core.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)==nrow(clade.clust.tmp)]
  assign(paste("clade", i, ".clust.core", sep=""), clade.clust.core.tmp)
  write.table(colnames(clade.clust.core.tmp), 
              file=paste("output/orthomcl-summary/clade-info/clade", i, "_clust_core.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
    
  # reduce clade core to those not in overall core
  clade.clust.core.red.tmp <- clade.clust.core.tmp[,!(colnames(clade.clust.core.tmp) %in% clust.core.all)]
  assign(paste("clade", i, ".clust.core.red", sep=""), clade.clust.core.red.tmp)
  write.table(colnames(clade.clust.core.red.tmp), 
              file=paste("output/orthomcl-summary/clade-info/clade", i, "_clust_core_red.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  print(ncol(clade.clust.core.tmp) - ncol(clade.clust.core.red.tmp))
  
  rm(clade.tmp, clade.clust.tmp, clade.clust.core.tmp, clade.clust.core.red.tmp)
}

## conserved, unique to each clade
clade1.clust.core.uniq <- colnames(clade1.clust.core.red)[!(colnames(clade1.clust.core.red) %in% c(colnames(clade2.clust), 
                                                                                                   colnames(clade3.clust), 
                                                                                                   colnames(clade4.clust)))]
write.table(clade1.clust.core.uniq, file="output/orthomcl-summary/clade-info/clade1_clust_core_uniq.txt",
            row.names=F, col.names=F, quote=F)

clade2.clust.core.uniq <- colnames(clade2.clust.core.red)[!(colnames(clade2.clust.core.red) %in% c(colnames(clade1.clust), 
                                                                                                   colnames(clade3.clust), 
                                                                                                   colnames(clade4.clust)))]
write.table(clade2.clust.core.uniq, file="output/orthomcl-summary/clade-info/clade2_clust_core_uniq.txt",
            row.names=F, col.names=F, quote=F)

clade3.clust.core.uniq <- colnames(clade3.clust.core.red)[!(colnames(clade3.clust.core.red) %in% c(colnames(clade1.clust), 
                                                                                                   colnames(clade2.clust), 
                                                                                                   colnames(clade4.clust)))]
write.table(clade3.clust.core.uniq, file="output/orthomcl-summary/clade-info/clade3_clust_core_uniq.txt",
            row.names=F, col.names=F, quote=F)

clade4.clust.core.uniq <- colnames(clade4.clust.core.red)[!(colnames(clade4.clust.core.red) %in% c(colnames(clade1.clust), 
                                                                                                   colnames(clade2.clust), 
                                                                                                   colnames(clade3.clust)))]
write.table(clade4.clust.core.uniq, file="output/orthomcl-summary/clade-info/clade4_clust_core_uniq.txt",
            row.names=F, col.names=F, quote=F)

##################################################
## clean up, save image
rm(i)
save.image("output/RData/orthomcl_summary.RData")