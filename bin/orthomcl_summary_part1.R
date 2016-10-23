## Parse OrthoMCL cluster information to obtain group membership, fasta IDs, and other useful information for downstream core/pangenome analysis of Gardnerella vaginalis
## Roxana Hickey <roxana.hickey@gmail.com>
## Last updated: 2015-05-21

# This script summarizes OrthoMCL gene clusters and outputs FASTA IDs and gene annotations for each cluster and singleton coding sequence (CDS). I have written this script to be sourced from another script that inputs the necessary data (e.g. run_orthomcl_summary_gv.R). That script should first be edited to specify theoutput directory, OrthoMCL cluster data, keys for FASTA IDs, gene features and annotations (in my case, from PATRIC CDS sequences), and custom color palettes.

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
dir.create(paste(dir.out, "/group-info", sep=""))

dir.create(paste(dir.out, "/clust-anno", sep=""))
dir.create(paste(dir.out, "/single-anno", sep=""))
dir.create(paste(dir.out, "/group-anno", sep=""))

for(i in c("ec", "go", "feat", "figfam", "path")){
  dir.create(paste(dir.out, "/clust-anno/", i, sep=""))
  dir.create(paste(dir.out, "/single-anno/", i, sep=""))
  dir.create(paste(dir.out, "/group-anno/", i, sep=""))
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
colnames(strain.uniq.df) <- c("sing_clust_count", "total_clust_count", "id_abbr")

# write singleton fasta IDs and annotation tables per strain
for(i in unique(ortho.sing$id_short)){
  # extract CDS IDs for strain singletons
  cds.tmp <- ortho.sing$cds_id[grep(i, ortho.sing$cds_id)]
  write.table(cds.tmp, 
              paste(dir.out, "/single-fasta-ids/", i, "_uniq_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # count number of singletons per strain
  strain.uniq.df[i, "sing_clust_count"] <- length(cds.tmp)
  
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
core.clust$less5 <- colnames(clust.mx)[colSums(clust.mx)>=nrow(clust.mx) - 5]

write.table(core.clust$all,
            file=paste(dir.out, "/group-info/core_clust_all.txt", sep=""), 
            row.names=F, col.names=F, quote=F)

write.table(core.clust$less1,
            file=paste(dir.out, "/group-info/core_clust_less1.txt", sep=""), 
            row.names=F, col.names=F, quote=F)

write.table(core.clust$less5,
            file=paste(dir.out, "/group-info/core_clust_less5.txt", sep=""), 
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

## identify clusters present in only one genome (not same as singletons)
clust.uniq.to.strain <- colnames(clust.mx.bin[,colSums(clust.mx.bin)==1])

##################################################
## clean up, save image
rm(i, j, key.tmp, clust.mx.freq)
save.image(paste(dir.out, "/orthomcl_summary_part1.RData", sep=""))