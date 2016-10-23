## Run orthomcl_summary.R for Gardnerella-only OrthoMCL (70% clustering)
## Roxana Hickey <roxana.hickey@gmail.com>
## Last update: 2015-04-23

## set up custom inputs and outputs
dir.out <- "results/orthomcl70-gv35"
dir.fig <- "results/figures/orthomcl70-gv35"

## read in ID key containing short and full strain names
genome.key <- read.table("data/reference/key_gv35_id.txt", sep="\t", header=T)

## read in orthoMCL groups
ortho.clust.in <- readLines("data/reference/orthomcl70_gv35_groups.txt")
ortho.sing.in <- readLines("data/reference/orthomcl70_gv35_single.txt")
clust.mx.in <- read.table("data/reference/orthomcl70_gv35_groups_binary.txt", header=T)

## key for fasta headers
key.fasta <- read.table("data/reference/key_gv35_fid_cds_fasta.txt", header=T, sep="\t")

## key for PATRIC gene features
key.feat <- read.table("data/reference/anno_gv35.features.tab", header=T, sep="\t", quote="")

# keep only CDS feature types
key.feat.cds <- subset(key.feat, feature_type=="CDS")
key.feat.cds <- merge(key.feat.cds, key.fasta[,c("cds_id", "locus_tag", "id_short")])

## key for EC, FigFam, GO and pathway annotations
key.anno <- list()
key.anno$ec <- read.table("data/reference/anno_gv35.ec", header=T, sep="\t", quote="")
key.anno$ec <- merge(key.fasta[,3:8], key.anno$ec)
key.anno$ec$anno_type <- rep("ec", nrow(key.anno$ec))

key.anno$figfam <- read.table("data/reference/anno_gv35.figfam", header=T, sep="\t", quote="")
key.anno$figfam <- merge(key.fasta[,3:8], key.anno$figfam)
key.anno$figfam$anno_type <- rep("figfam", nrow(key.anno$figfam))

key.anno$go <- read.table("data/reference/anno_gv35.go", header=T, sep="\t", quote="")
key.anno$go <- merge(key.fasta[,3:8], key.anno$go)
key.anno$go$anno_type <- rep("go", nrow(key.anno$go))

key.anno$path <- read.table("data/reference/anno_gv35.path", header=T, sep="\t", quote="")
key.anno$path <- merge(key.fasta[,3:8], key.anno$path)
key.anno$path$anno_type <- rep("path", nrow(key.anno$path))

## pick custom colors
library(wesanderson)
col.cust <- wes_palette("Darjeeling", 5, type="discrete")[1:5]
col.cust.hc <- col.cust[-4]
col.bw.seq <- "YlOrRd"

## run orthomcl_summary.R
source("code/orthomcl_summary.R")

## RData will be written to results/orthomcl-summary-gv/ortho_summary.RData