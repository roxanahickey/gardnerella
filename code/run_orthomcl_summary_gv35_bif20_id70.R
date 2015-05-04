## Run orthomcl_summary.R for Gardnerella (35 genomes) + Bifidobacterium (20 genomes) OrthoMCL (50% clustering)
## Roxana Hickey <roxana.hickey@gmail.com>
## Last updated: 2015-04-23

## set up custom inputs and outputs
dir.out <- "results/orthomcl70-gv35-bif20"
dir.fig <- "results/figures/orthomcl70-gv35-bif20"

## read in ID key containing short and full strain names
genome.key <- read.table("data/reference/key_gv35_bif20_id.txt", sep="\t", header=T)

## read in orthoMCL groups
ortho.clust.in <- readLines("data/reference/orthomcl70_gv35_bif20_groups.txt")
ortho.sing.in <- readLines("data/reference/orthomcl70_gv35_bif20_single.txt")
clust.mx.in <- read.table("data/reference/orthomcl70_gv35_bif20_groups_binary.txt", header=T)

## key for fasta headers
key.fasta.gv <- read.table("data/reference/key_gv35_fid_cds_fasta.txt", header=T, sep="\t")
key.fasta.bif <- read.table("data/reference/key_bif20_fid_cds_fasta.txt", header=T, sep="\t")
key.fasta <- rbind(key.fasta.gv, key.fasta.bif)

## key for PATRIC gene features
key.feat.gv <- read.table("data/reference/anno_gv35.features.tab", header=T, sep="\t", quote="")
key.feat.bif <- read.table("data/reference/anno_bif20.features.tab", header=T, sep="\t", quote="")
key.feat <- rbind(key.feat.gv, key.feat.bif)

# keep only CDS feature types
key.feat.cds <- subset(key.feat, feature_type=="CDS")
key.feat.cds <- merge(key.feat.cds, key.fasta[,c("cds_id", "locus_tag", "id_short")])

## key for EC, FigFam, GO and pathway annotations
key.anno <- list()
key.ec.gv <- read.table("data/reference/anno_gv35.ec", header=T, sep="\t", quote="")
key.ec.bif <- read.table("data/reference/anno_bif20.ec", header=T, sep="\t", quote="")
key.anno$ec <- rbind(key.ec.gv, key.ec.bif)
key.anno$ec <- merge(key.fasta[,3:8], key.anno$ec)
key.anno$ec$anno_type <- rep("ec", nrow(key.anno$ec))

key.figfam.gv <- read.table("data/reference/anno_gv35.figfam", header=T, sep="\t", quote="")
key.figfam.bif <- read.table("data/reference/anno_bif20.figfam", header=T, sep="\t", quote="")
key.anno$figfam <- rbind(key.figfam.gv, key.figfam.bif)
key.anno$figfam <- merge(key.fasta[,3:8], key.anno$figfam)
key.anno$figfam$anno_type <- rep("figfam", nrow(key.anno$figfam))

key.go.gv <- read.table("data/reference/anno_gv35.go", header=T, sep="\t", quote="")
key.go.bif <- read.table("data/reference/anno_bif20.go", header=T, sep="\t", quote="")
key.anno$go <- rbind(key.go.gv, key.go.bif)
key.anno$go <- merge(key.fasta[,3:8], key.anno$go)
key.anno$go$anno_type <- rep("go", nrow(key.anno$go))

key.path.gv <- read.table("data/reference/anno_gv35.path", header=T, sep="\t", quote="")
key.path.bif <- read.table("data/reference/anno_bif20.path", header=T, sep="\t", quote="")
key.anno$path <- rbind(key.path.gv, key.path.bif)
key.anno$path <- merge(key.fasta[,3:8], key.anno$path)
key.anno$path$anno_type <- rep("path", nrow(key.anno$path))

## pick custom colors
library(wesanderson)
col.cust <- c(wes_palette("Darjeeling", 5, type="discrete")[1:5],
              wes_palette("FantasticFox", 5, type="discrete")[2],
              wes_palette("Rushmore", 5, type="discrete")[4])
col.cust.hc <- col.cust
col.bw.seq <- "YlGnBu"

## run orthomcl_summary.R
source("code/orthomcl_summary.R")

## RData will be written to results/orthomcl-summary-gv-bif/ortho_summary.RData