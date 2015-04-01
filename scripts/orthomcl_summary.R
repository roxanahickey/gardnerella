## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information
## for downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-13

library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(reshape)
source("scripts/ggplot_cust.R")

##################################################
## create directories for output data
dir.create(dir.out)
dir.create(dir.fig)

dir.create(paste(dir.out, "/clust-fasta-ids", sep=""))
dir.create(paste(dir.out, "/single-fasta-ids", sep=""))
dir.create(paste(dir.out, "/clade-info", sep=""))

##################################################
## get info from OrthoMCL clusters
ortho.clust <- list()
ortho.clust$orig <- ortho.clust.in

# list of clusters and fasta IDs per cluster
ortho.clust$ortho_id <- unlist(strsplit(ortho.clust$orig, split=": "))[seq(2, length(ortho.clust$orig)*2, by=2)]
names(ortho.clust$ortho_id) <- unlist(strsplit(ortho.clust$orig, split=": "))[seq(1, length(ortho.clust$orig)*2, by=2)]

# write lists of fasta IDs per cluster
for(i in unique(names(ortho.clust$ortho_id))){
  clust.tmp <- data.frame(strsplit(ortho.clust$ortho_id[i], " "))
  write.table(clust.tmp, paste(dir.out, "/clust-fasta-ids/", i, "_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

# clean up
rm(i, clust.tmp)

##################################################
## get info from OrthoMCL singletons
ortho.sing <- list()

# list of singleton fasta IDs
ortho.sing$cds_id <- ortho.sing.in

# GV strain names
ortho.sing$id_short <- unlist(strsplit(ortho.sing$cds_id, split="|", fixed=T))[seq(1, length(ortho.sing$cds_id)*2, by=2)]

# initialize empty data frame to keep track of counts of singleton clusters per strain (more info to be added later)
strain.uniq.df <- data.frame(matrix(0, nrow=length(unique(ortho.sing$id_short)), ncol=3))
rownames(strain.uniq.df) <- sort(unique(ortho.sing$id_short))
colnames(strain.uniq.df) <- c("uniq_clust_count", "total_clust_count", "id_abbr")

# write lists of singleton fasta IDs per strain
for(i in unique(ortho.sing$id_short)){
  sing.tmp <- ortho.sing$cds_id[grep(i, ortho.sing$cds_id)]
#   assign(paste(i, ".uniq", sep=""), sing.tmp)
  
  strain.uniq.df[i, "uniq_clust_count"] <- length(sing.tmp)
  
  write.table(sing.tmp, paste(dir.out, "/single-fasta-ids/", i, "_uniq_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  }

# clean up
rm(i, sing.tmp)

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

##################################################
## run hclust to identify groups
source("scripts/orthomcl_hclust_id_clades.R")

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
  
  write.table(clade.tmp, 
              file=paste(dir.out, "/clade-info/clade", i, "_id.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

rm(i, j, clade.tmp)

clade.list <- unlist(clade.list, use.names=F)

clade.clust <- list()
clade.core.clust <- list()
clade.core.clust.red <- list()

for(i in unique(clade.list)){

  id.tmp <- unlist(clade.id[i], use.names=F)
  
  # reduce clust.mx to clade members
  clade.clust.tmp <- clust.mx[id.tmp,]
  clade.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)>=1]  
  
  clade.clust[i] <- list(colnames(clade.clust.tmp))
  
  write.table(colnames(clade.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract clade core clusters
  clade.core.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)==nrow(clade.clust.tmp)]

  clade.core.clust[i] <- list(colnames(clade.core.clust.tmp))

  write.table(colnames(clade.core.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
    
  # reduce clade core to those not in overall core
  clade.core.clust.red.tmp <- clade.core.clust.tmp[,!(colnames(clade.core.clust.tmp) %in% core.clust$all)]

  clade.core.clust.red[i] <- list(colnames(clade.core.clust.red.tmp))
  
  write.table(colnames(clade.core.clust.red.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core_red.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
    
  # fill in summary table
  clust.smy[i, "all_core"] <- ncol(clade.core.clust.tmp) - ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_core_red"] <- ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_sing"] <- length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  clust.smy[i, "clade_total"] <- ncol(clade.clust.tmp) + length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  
}

rm(i, id.tmp, clade.clust.tmp, clade.core.clust.tmp, clade.core.clust.red.tmp)

clust.smy$clade_pan <- clust.smy$clade_total - (clust.smy$all_core + clust.smy$clade_core_red + clust.smy$clade_sing)

##################################################
## clade core unique genes
# initialize empty list
clade.core.clust.uniq <- list()

for(i in unique(clade.list)){
  ignore <- clade.list[!(clade.list %in% i)]
  ignore.clust <- unique(unlist(clade.clust[ignore], use.names = F))
  query.clust <- unlist(clade.core.clust.red[i], use.names = F)
  
  clade.core.clust.uniq[i] <- list(query.clust[!(query.clust %in% ignore.clust)])
  
  clust.smy[i, "clade_core_uniq"] <- length(query.clust[!(query.clust %in% ignore.clust)])
  
  write.table(query.clust[!(query.clust %in% ignore.clust)], 
              file=paste(dir.out, "/clade-info/", i, "_clust_core_uniq.txt", sep=""),
              row.names=F, col.names=F, quote=F)
  
  rm(i, ignore, ignore.clust, query.clust)
}

# finish filling out summary table
clust.smy$clade_core_not_uniq <- clust.smy$clade_core_red - clust.smy$clade_core_uniq

# reduce summary table to make nonredundant
clust.smy.red <- clust.smy[,-c(4,7)]
clust.smy.red$clade_id <- rownames(clust.smy.red)
clust.smy.red.md <- melt(clust.smy.red, id="clade_id")

write.table(clust.smy.red, file=paste(dir.out, "/clade-info/clust_count_summary.txt", sep=""), sep="\t", quote=F)

##################################################
## singletons per genome -- complete data frame started earlier
strain.uniq.df$total_clust_count <- rowSums(clust.mx)[match(rownames(strain.uniq.df), names(rowSums(clust.mx)))] + strain.uniq.df$uniq_clust_count

strain.uniq.df$uniq_clust_prop <- strain.uniq.df$uniq_clust_count / strain.uniq.df$total_clust_count

strain.uniq.df$id_abbr <- genome.key$id_abbr[match(rownames(strain.uniq.df), genome.key$id_short)]

strain.uniq.df$ortho_clade <- genome.key$ortho_clade[match(rownames(strain.uniq.df), genome.key$id_short)]

##################################################
## plots!
# number of gene clusters present in exactly n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=freq, fill=no_genomes)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("No. genomes") + ylab("Frequency of gene clusters") +
  geom_text(aes(label=freq, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, col.bw.seq))(length(clust.ct$no_genomes))) +
  theme(legend.position="none")
ggsave(paste(dir.fig, "/orthomcl_clusters_n_genome_exact.png", sep=""), width=6, height=4, units="in", scale=2)

# number of gene clusters present in at least n genomes
ggplot(clust.ct, aes(x=factor(no_genomes), y=rev(cumsum(rev(clust.ct$freq))), fill=no_genomes)) + 
  geom_bar(stat="identity") + 
  xlab("No. genomes") + ylab("Frequency of gene clusters") + theme_cust + 
  geom_text(aes(label=rev(cumsum(rev(clust.ct$freq))), angle=45, hjust=0, vjust=-0.5, size=0.8)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(5, col.bw.seq))(length(clust.ct$no_genomes))) +
  theme(legend.position="none")
ggsave(paste(dir.fig, "/orthomcl_clusters_n_genome_atleast.png", sep=""), width=6, height=4, units="in", scale=2)

# cluster counts per clade
ggplot(clust.smy.red.md, aes(x=clade_id, y=value, fill=variable)) + 
  geom_bar(stat="identity") + theme_cust +
  xlab("Clade") + ylab("Frequency of gene clusters") +
  labs(fill="") +
  scale_fill_manual(limits=c("all_core", "clade_core_not_uniq", "clade_core_uniq", "clade_pan", "clade_sing"),
                    labels=c("All core (all genomes, all clades)", "Clade core, non-unique (all genomes in clade)", 
                             "Clade core, unique (all genomes in clade)", "Clade pan (2+ genomes in clade)", 
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