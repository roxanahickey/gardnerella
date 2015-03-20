## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information
## for downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-13

library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(reshape)
source("scripts/ggplot_cust.R")

## create directories for output data
dir.create("output/orthomcl-summary")
dir.create("output/RData")

## read in ID key containing short and full strain names
gv.id <- read.table("data/key_gv_patric_all35_id.txt", sep="\t", header=T)

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

# initialize empty data frame to keep track of counts of singleton clusters per strain (more info to be added later)
strain.uniq.df <- data.frame(matrix(0, nrow=length(unique(ortho.sing.strain)), ncol=3))
rownames(strain.uniq.df) <- sort(unique(ortho.sing.strain))
colnames(strain.uniq.df) <- c("uniq_clust_count", "total_clust_count", "id_strain")

# write lists of singleton fasta IDs per strain
for(i in unique(ortho.sing.strain)){
  sing.tmp <- ortho.sing[grep(i, ortho.sing)]
  assign(paste(i, ".uniq", sep=""), sing.tmp)
  
  strain.uniq.df[i, "uniq_clust_count"] <- length(sing.tmp)
  
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

# because I'm obsessed with plot colors and wesanderson palettes are ah-mazing
wes.dj.orig <- wes_palette("Darjeeling", 5, type="discrete")[1:5]
# wes.dj.reorder <- wes_palette("Darjeeling", 5, type="discrete")[c(5,2,3,4,1)]

# plot number of gene clusters present in exactly n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=freq, fill=no_genomes)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("No. genomes") + ylab("No. gene clusters") +
  geom_text(aes(label=freq, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, "YlOrRd"))(35)) +
  theme(legend.position="none")
ggsave("figs/orthomcl_clusters_n_genome_exact.png", width=6, height=4, units="in", scale=2)

# plot number of gene clusters present in at least n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=rev(cumsum(rev(clust.ct$freq))), fill=no_genomes)) + 
  geom_bar(stat="identity") + 
  xlab("No. genomes") + ylab("No. gene clusters") + theme_cust + 
  geom_text(aes(label=rev(cumsum(rev(clust.ct$freq))), angle=45, hjust=0, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, "YlOrRd"))(35)) +
  theme(legend.position="none")
ggsave("figs/orthomcl_clusters_n_genome_atleast.png", width=6, height=4, units="in", scale=2)

##################################################
# get clust IDs of core gene clusters (+ "relaxed" core)
clust.core.all <- colnames(clust.mx)[colSums(clust.mx)==nrow(clust.mx)]
clust.core.less1 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 1]
clust.core.less2 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 2]

##################################################
# run hclust to identify groups
col.cust <- wes.dj.orig[-4]
source("scripts/orthomcl_hclust_id_clades.R")

# add ortho clade id to gv.id
gv.id$ortho_clade <- cutg[gv.id$id_short]

##################################################
## clade-unique core
dir.create("output/orthomcl-summary/clade-info")

# intialize empty data frame to fill as you calculate clusters
clust.smy <- data.frame(matrix(0, nrow=length(unique(cutg)), ncol=7))
rownames(clust.smy) <- unique(cutg)
colnames(clust.smy) <- c("species_core", "clade_core_not_uniq", "clade_core_uniq",
                         "clade_core_red", "clade_pan", "clade_sing", "clade_total")

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
    
  # fill in summary table
  clust.smy$species_core[i] <- ncol(clade.clust.core.tmp) - ncol(clade.clust.core.red.tmp)
  clust.smy$clade_core_red[i] <- ncol(clade.clust.core.red.tmp)
  clust.smy$clade_sing[i] <- length(clade.clust.tmp[colnames(clade.clust.tmp) %in% ortho.sing.id])
  clust.smy$clade_total[i] <- ncol(clade.clust.tmp)
  
  rm(clade.tmp, clade.clust.tmp, clade.clust.core.tmp, clade.clust.core.red.tmp)
}

clust.smy$clade_pan <- clust.smy$clade_total - (clust.smy$species_core + clust.smy$clade_core_red + clust.smy$clade_sing)

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

# finish filling out summary table
clust.smy$clade_core_uniq <- c(length(clade1.clust.core.uniq), length(clade2.clust.core.uniq),
                               length(clade3.clust.core.uniq), length(clade4.clust.core.uniq))
clust.smy$clade_core_not_uniq <- clust.smy$clade_core_red - clust.smy$clade_core_uniq

# reduce summary table to make nonredundant
clust.smy.red <- clust.smy[,-c(4,7)]
clust.smy.red$clade_id <- rownames(clust.smy.red)
clust.smy.red.md <- melt(clust.smy.red, id="clade_id")

# plot
ggplot(clust.smy.red.md, aes(x=clade_id, y=value, fill=variable)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("Cluster") + ylab("Number of gene clusters") +
  labs(fill="") +
  scale_fill_manual(limits=c("species_core", "clade_core_not_uniq", "clade_core_uniq", "clade_pan", "clade_sing"),
                    labels=c("Species core (35 genomes)", "Non-unique core (all genomes in clade)", "Unique core (all genomes in clade)",
                             "Pan (2+ genomes in clade)", "Singletons (1 genome)"), 
                    values=rev(brewer.pal(5, "YlOrRd")))
ggsave("figs/orthomcl_clusters_by_clade.png", width=6, height=4, units="in", scale=2)

##################################################
## singletons per genome -- complete data frame started earlier
strain.uniq.df$total_clust_count <- rowSums(clust.mx)[match(rownames(strain.uniq.df), names(rowSums(clust.mx)))] + strain.uniq.df$uniq_clust_count

strain.uniq.df$uniq_clust_prop <- strain.uniq.df$uniq_clust_count / strain.uniq.df$total_clust_count

strain.uniq.df$id_strain <- gv.id$id_strain[match(rownames(strain.uniq.df), gv.id$id_short)]

strain.uniq.df$ortho_clade <- gv.id$ortho_clade[match(rownames(strain.uniq.df), gv.id$id_short)]

# sort by clade, proportion of unique gene clusters
strain.uniq.df.sorted <- strain.uniq.df$id_strain[with(strain.uniq.df, order(ortho_clade, -uniq_clust_prop))]

# plot number of singleton clusters per genome
ggplot(strain.uniq.df, aes(x=reorder(id_strain, ortho_clade - uniq_clust_prop), y=uniq_clust_prop, fill=factor(ortho_clade))) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("Genome") + ylab("Proportion of total gene clusters per genome") + 
  theme(axis.text.x = element_text(angle=-30, hjust=0, vjust=1)) +
  labs(fill="Cluster") +
  geom_text(aes(label=uniq_clust_count, vjust=-0.5, size=1), show_guide=F) +
#   scale_fill_manual(values=wes_palette("Darjeeling", 4, type="discrete"))
  scale_fill_manual(values=wes.dj.orig[-4])
ggsave("figs/orthomcl_clusters_uniq_per_strain.png", width=6, height=4, units="in", scale=2)

##################################################
## clean up, save image
rm(i)
save.image("output/RData/orthomcl_summary.RData")