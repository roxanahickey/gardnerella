## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information or downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## Last updated: 2016-11-30

# This script follows part 1 and the step in which genome groups are specified. This can be done using hierarchical clustering and/or phylogenetic assessment. In the case of 35 Gardnerella vaginalis genomes, I performed both hierarchical clustering of ortholog presence/absence and phylogenetic concordance analysis of the core genes, settling on a grouping of 3 major clades (6 subclades). This script follows that by extracting FASTA IDs and annotations for each clade. Finally, it produces summary figures of core and pan genome CDS counts.

## Note to self 5/26/15: did not account for clusters that were unique to a single genome when I parsed singletons vs. clusters. There are 7 clusters that are only found in one genome each (edited part1 to define this as clust.uniq.to.strain). May want to return to this later.

library(ggplot2)
library(RColorBrewer)
library(reshape)

##################################################
# add group and clade ids to genome.key
raxml.clade.id <- read.table("data/reference/raxml_consensus_clade.txt", 
                             sep="\t", header=T)

genome.key <- merge(genome.key, raxml.clade.id)

# list groups (clades)
cd <- unlist(raxml.clade.id$raxml_clade)
names(cd) <- unlist(raxml.clade.id$id_short)

sc <- unlist(raxml.clade.id$raxml_subclade)
names(sc) <- unlist(raxml.clade.id$id_short)

# pick clade or subclade
# gp <- cd
gp <- sc

##################################################
## group core genes
# intialize empty data frame to fill as you calculate clusters
clust.smy <- data.frame(matrix(0, nrow=length(unique(gp)), ncol=7))
rownames(clust.smy) <- paste("group", unique(gp), sep="")
colnames(clust.smy) <- c("all_core", "group_core_not_uniq", "group_core_uniq",
                         "group_core_red", "group_acc", "group_sing", "group_total")

# initialize empty lists to fill in
group.id <- list()
group.list <- list()

for(i in unique(gp)){
  # list group member ids
  group.tmp <- names(gp[gp==i])
  j <- paste("group", i, sep="")
  group.id[j] <- list(group.tmp)
  group.list[i] <- paste("group", i, sep="")
  
  # write list of strain IDs in cluster
  write.table(group.tmp, 
              file=paste(dir.out, "/group-info/group", i, "_id.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

rm(i, j, group.tmp)

group.list <- unlist(group.list, use.names=F)

group.clust <- list()
group.core.clust <- list()
group.core.clust.red <- list()

# write new get_anno function to work with groups -- includes j term to specify cluster summarized
get_anno_group <- function(x) {anno.tmp <- x[x$cds_id %in% unlist(cds.tmp, use.names=F),]
                               anno.tmp <- merge(anno.tmp, key.ortho, all.x=T)
                               anno.type <- unique(x$anno_type)
                               if(nrow(anno.tmp) > 0) {write.table(anno.tmp, 
                                                                   paste(dir.out.anno, "/", anno.type, "/", i, "_", j, "_", anno.type, ".txt", sep=""), 
                                                                   sep="\t", row.names=F, quote=F)}}

for(i in unique(group.list)){
  
  id.tmp <- unlist(group.id[i], use.names=F)
  
  ## ALL clusters in group
  # reduce clust.mx.bin to group members
  group.clust.tmp <- clust.mx.bin[id.tmp,]
  group.clust.tmp <- group.clust.tmp[,colSums(group.clust.tmp)>=1]  
  
  # append list to group.clust
  group.clust[i] <- list(colnames(group.clust.tmp))
  
  # write list of cluster IDs in group
  write.table(colnames(group.clust.tmp), 
              file=paste(dir.out, "/group-info/", i, "_clust.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/group-anno", sep="")
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(group.clust.tmp)], " ")
  j <- "clust"
  
  # extract PATRIC features for CDS IDs in clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract GO and KEGG pathways for CDS IDs in clusters
  sapply(key.anno, get_anno_group)
  
  ## CORE clusters in group
  # reduce group.clust.tmp to core clusters
  group.core.clust.tmp <- group.clust.tmp[,colSums(group.clust.tmp)==nrow(group.clust.tmp)]
  
  # append to group.core.clust
  group.core.clust[i] <- list(colnames(group.core.clust.tmp))
  
  # write list of cluster IDs in core set in group
  write.table(colnames(group.core.clust.tmp), 
              file=paste(dir.out, "/group-info/", i, "_clust_core.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(group.core.clust.tmp)], " ")
  j <- "clust_core"
  
  # extract PATRIC features for CDS IDs in core clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract GO and KEGG pathways for CDS IDs in core clusters
  sapply(key.anno, get_anno_group)
  
  ## REDUCED CORE clusters in group (i.e. not also present in core for all strains)
  # reduce group core to those not in overall core
  group.core.clust.red.tmp <- group.core.clust.tmp[,!(colnames(group.core.clust.tmp) %in% core.clust$all)]
  
  # append to group.core.clust.red
  group.core.clust.red[i] <- list(colnames(group.core.clust.red.tmp))
  
  # write list of cluster IDs in reduced core set in group
  write.table(colnames(group.core.clust.red.tmp), 
              file=paste(dir.out, "/group-info/", i, "_clust_core_red.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract CDS IDs from clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[colnames(group.core.clust.red.tmp)], " ")
  j <- "clust_core_red"
  
  # extract PATRIC features for CDS IDs in reduced core clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_red_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in reduced core clusters
  sapply(key.anno, get_anno_group)
  
  # fill in summary table
  clust.smy[i, "all_core"] <- ncol(group.core.clust.tmp) - ncol(group.core.clust.red.tmp)
  clust.smy[i, "group_core_red"] <- ncol(group.core.clust.red.tmp)
  clust.smy[i, "group_sing"] <- length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  clust.smy[i, "group_total"] <- ncol(group.clust.tmp) + length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  
}

rm(i, j, id.tmp, group.clust.tmp, group.core.clust.tmp, group.core.clust.red.tmp, dir.out.anno, cds.tmp, feat.tmp)

clust.smy$group_acc <- clust.smy$group_total - (clust.smy$all_core + clust.smy$group_core_red + clust.smy$group_sing)

##################################################
## group core unique genes
# initialize empty list
group.core.clust.uniq <- list()

for(i in unique(group.list)){
  # select focal group, ignore others
  ignore <- group.list[!(group.list %in% i)]
  ignore.clust <- unique(unlist(group.clust[ignore], use.names = F))
  query.clust <- unlist(group.core.clust.red[i], use.names = F)
  
  # append unique core cluster IDs to group.core.clust.uniq
  group.core.clust.uniq[i] <- list(query.clust[!(query.clust %in% ignore.clust)])
  
  # count length of core unique cluster IDs, append to summary table
  clust.smy[i, "group_core_uniq"] <- length(query.clust[!(query.clust %in% ignore.clust)])
  
  # write list of core unique cluster IDs in group
  write.table(query.clust[!(query.clust %in% ignore.clust)], 
              file=paste(dir.out, "/group-info/", i, "_clust_core_uniq.txt", sep=""),
              row.names=F, col.names=F, quote=F)
  
  # specify output directory for annotation tables
  dir.out.anno <- paste(dir.out, "/group-anno", sep="")
  
  # extract CDS IDs from core unique clusters and specify cluster summary type
  cds.tmp <- strsplit(ortho.clust$ortho_id[query.clust[!(query.clust %in% ignore.clust)]], " ")
  j <- "clust_core_uniq"
  
  # extract PATRIC features for CDS IDs in clusters
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(cds.tmp, use.names=F),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, 
                                      paste(dir.out.anno, "/feat/", i, "_clust_core_uniq_feat.txt", sep=""), 
                                      sep="\t", row.names=F, quote=F)}
  
  # extract EC, FigFam, GO and KEGG pathways for CDS IDs in clusters
  sapply(key.anno, get_anno_group)
}

rm(i, j, ignore, ignore.clust, query.clust, dir.out.anno, cds.tmp, feat.tmp)

# finish filling out summary table
clust.smy$group_core_not_uniq <- clust.smy$group_core_red - clust.smy$group_core_uniq

# reduce summary table to make nonredundant
clust.smy.red <- clust.smy[,-c(4,7)]
clust.smy.red$group_id <- rownames(clust.smy.red)
clust.smy.red.md <- melt(clust.smy.red, id="group_id")
write.table(clust.smy.red, file=paste(dir.out, "/group-info/clust_count_summary.txt", sep=""), sep="\t", quote=F)

##################################################
## singletons per genome -- complete data frame started earlier
strain.uniq.df$total_clust_count <- rowSums(clust.mx)[match(rownames(strain.uniq.df), names(rowSums(clust.mx)))] + strain.uniq.df$sing_clust_count

strain.uniq.df$sing_clust_prop <- strain.uniq.df$sing_clust_count / strain.uniq.df$total_clust_count

strain.uniq.df$id_abbr <- genome.key$id_abbr[match(rownames(strain.uniq.df), genome.key$id_short)]

strain.uniq.df$group_id <- gp[match(rownames(strain.uniq.df), names(gp))]

## obtain genome characteristics
genome.char <- read.table("data/reference/gv_genome_char_manual_calc.txt", header=T, sep="\t")
genome.key.char <- merge(genome.key, genome.char)

strain.uniq.df$id_short <- rownames(strain.uniq.df)

genome.key.char <- merge(genome.key.char, strain.uniq.df)

## accessory genome
group.acc.clust <- list()

for(i in names(group.clust)){
  group.acc.clust[i] <- list(unlist(group.clust[i], use.names = FALSE)[!(unlist(group.clust[i]) %in% unlist(group.core.clust[i]))])
}

## count number of core/accessory/unique orthologs per genome (count unique only, not number of occurrences)
strain.clust.count <-  data.frame(matrix(0, nrow=length(unique(genome.key$id_short)), ncol=6))
rownames(strain.clust.count) <- sort(unique(genome.key$id_short))
colnames(strain.clust.count) <- c("all_core", "group_core_not_uniq", "group_core_uniq", "accessory", "singleton", "group_id")

for(i in unique(genome.key$id_short)){
  ## reduce ortholog matrix to focal genome
  clust.tmp <- clust.mx[i,]
  
  ## keep orthologs present in focal genome
  pick <- clust.tmp[clust.tmp > 0]
  
  ## select group (clade) that contains focal genome
  group.tmp <- paste0("group", as.character(genome.key$raxml_subclade[genome.key$id_short==i]))
  
  ## number of accessory orthologs in focal genome
  acc.tmp <- pick[names(pick) %in% unlist(group.acc.clust[group.tmp], use.names = FALSE)]
  acc.len.tmp <- length(acc.tmp)
  
  ## number of unique core orthologs in focal genome
  core.uniq.tmp <- pick[unlist(group.core.clust.uniq[group.tmp], use.names = FALSE)]
  core.uniq.len.tmp <- length(core.uniq.tmp)
  
  ## number of non-unique core orthologs in focal genome
  core.red.tmp <- pick[unlist(group.core.clust.red[group.tmp], use.names = FALSE)]
  core.nonuniq.len.tmp <- length(core.red.tmp) - core.uniq.len.tmp
  
  ## number of singleton CDS in focal genome
  sing.len.tmp <- strain.uniq.df$sing_clust_count[rownames(strain.uniq.df)==i]
  
  ## number of core orthologs in species
  core.all.tmp <- length(core.clust$all)
  
  ## write to data frame
  strain.clust.count[i,"all_core"] <- core.all.tmp
  strain.clust.count[i,"group_core_not_uniq"] <- core.nonuniq.len.tmp
  strain.clust.count[i,"group_core_uniq"] <- core.uniq.len.tmp
  strain.clust.count[i,"accessory"] <- acc.len.tmp
  strain.clust.count[i,"singleton"] <- sing.len.tmp
  strain.clust.count[i,"group_id"] <- group.tmp
}

strain.clust.count$id_abbr <- genome.key$id_abbr[match(rownames(strain.clust.count), genome.key$id_short)]
strain.clust.count <- strain.clust.count[order(strain.clust.count$group_id),]
write.table(strain.clust.count, file=paste0(dir.out, "/group-info/strain_clust_count_summary.txt"), 
            sep="\t", quote=F)

strain.clust.count.md <- melt(strain.clust.count, id=c("id_abbr", "group_id"))

strain.clust.prop <- prop.table(as.matrix(strain.clust.count[,-c(6,7)]), 1)
write.table(strain.clust.prop, file=paste0(dir.out, "/group-info/strain_clust_prop_summary.txt"), 
            sep="\t", quote=F)

##################################################
## clean up, save image
save.image(paste(dir.out, "/orthomcl_summary_part2.RData", sep=""))
