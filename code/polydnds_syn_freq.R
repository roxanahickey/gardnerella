## polydnds synonymous site parsing
## 2015-04-29 RJH

#####################################
## OrthoMCL 50%, 35 Gardnerella + 20 Bifidobacterium
# setwd("~/Documents/research/gardnerella/gardnerella-git/results/polydnds-orthomcl50-gv35-bif20/")
# file.prefix <- "orthomcl50_gv35_bif20"
# poly.syn.gv.bif <- read.table(paste(file.prefix, "_core_syn_sites_bases.txt", sep=""), 
#                               sep="\t", header=F, row.names=1, stringsAsFactors = F)
# 
# poly.syn.gv <- poly.syn.gv.bif[grep("GV", rownames(poly.syn.gv.bif)),]
# poly.syn.bif <- poly.syn.gv.bif[grep("GV", rownames(poly.syn.gv.bif), invert=T),]

# pick one
# poly.syn <- poly.syn.gv
# file.prefix <- "orthomcl50_gv35"; col.fill <- "dodgerblue3"
# plot.title <- "35 Gardnerella genomes, OrthoMCL 50% (with Bifidobacterium spp.)"
# poly.syn <- poly.syn.bif

#####################################
## OrthoMCL 70%, 35 Gardnerella
setwd("~/Documents/research/gardnerella/gardnerella-git/results/polydnds-orthomcl70-gv35/")
file.prefix <- "orthomcl70_gv35"
poly.syn.gv <- read.table(paste(file.prefix, "_core_syn_sites_bases.txt", sep=""), 
                          sep="\t", header=F, row.names=1, stringsAsFactors = F)
poly.syn <- poly.syn.gv
plot.title <- "35 Gardnerella genomes, OrthoMCL 70%"; col.fill <- "gray20"

#####################################
## SEPARATE BY CLADE
load("~/Documents/research/gardnerella/gardnerella-git/results/orthomcl70-gv35/orthomcl_summary.RData")

# pick one
poly.syn <- poly.syn.gv[rownames(poly.syn.gv) %in% clade.id$clade1,]
file.prefix <- "orthomcl50_gv35_clade1"; col.fill <- col.cust.hc[1]
plot.title <- "11 Gardnerella genomes, OrthoMCL 70%, clade 1"

poly.syn <- poly.syn.gv[rownames(poly.syn.gv) %in% clade.id$clade2,]
file.prefix <- "orthomcl50_gv35_clade2"; col.fill <- col.cust.hc[2]
plot.title <- "8 Gardnerella genomes, OrthoMCL 70%, clade 2"

poly.syn <- poly.syn.gv[rownames(poly.syn.gv) %in% clade.id$clade3,]
file.prefix <- "orthomcl50_gv35_clade3"; col.fill <- col.cust.hc[3]
plot.title <- "14 Gardnerella genomes, OrthoMCL 70%, clade 3"

#####################################
## Ensure all sites have two alleles
allele.check <- apply(poly.syn, 2, FUN=function(x){length(unique(x))})
poly.syn.multi <- poly.syn[,allele.check>2]
poly.syn <- poly.syn[,allele.check==2]
sum(allele.check[allele.check==2]) == ncol(poly.syn)*2 ## check if TRUE
colnames(poly.syn) <- seq(1, ncol(poly.syn), by=1)

#####################################
## Identify minor/major alleles and replace with 0/1
df <- apply(poly.syn, 2, FUN = function(x) {
  allele <- unique(x)
  a <- length(x[x == allele[1]])
  b <- length(x[x == allele[2]])
  
  if(a < b){
    x[x == allele[1]] <- rep(1, length(x[x == allele[1]])) # replace minor with 1
    x[x == allele[2]] <- rep(0, length(x[x == allele[2]])) # replace major with 0
    x <- as.numeric(x)
  } else {
    x[x == allele[2]] <- rep(1, length(x[x == allele[2]])) # replace minor with 1
    x[x == allele[1]] <- rep(0, length(x[x == allele[1]])) # replace major with 0
    x <- as.numeric(x)
  }
})

rownames(df) <- rownames(poly.syn)
write.table(df, paste(file.prefix, "_core_syn_sites_binary.txt", sep=""), 
            sep="\t", row.names=T, quote=F, col.names=F)

## folded SFS
df.colsum <- colSums(df)
df.factor <- factor(df.colsum, levels=seq(1, nrow(df)/2, 1))
sfs.folded <- as.data.frame(table(df.factor))
colnames(sfs.folded) <- c("minor_count", "freq")

# plot(sfs.folded$freq ~ sfs.folded$minor_count)

library(ggplot2)
source("../../code/ggplot_cust.R")
ggplot(sfs.folded, aes(x=as.numeric(minor_count), y=freq)) + 
  geom_bar(stat="identity", fill = col.fill) +
  xlab("Minor allele frequency") +
  scale_x_continuous(breaks=seq_len(nrow(sfs.folded))) +
  ylab("Number of synonymous segregating sites") +
  ggtitle(print(plot.title)) +
  theme_cust_nominor
ggsave(paste(file.prefix, "_sfs_folded_syn.png", sep=""), 
       width = 8, height = 5, units = "in")

save.image(paste(file.prefix, "_polydnds_syn_freq.RData", sep=""))
