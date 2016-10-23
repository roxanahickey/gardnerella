## Run orthomcl_to_fasta.R for Gardnerella + Bifidobacterium OrthoMCL
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-04-28

library(Biostrings)

## load output from orthomcl_summary.R
load("results/orthomcl50-gv35-bif20/orthomcl_summary.RData")

## clear existing output paths
rm(dir.fig, dir.out)

## set up custom inputs and outputs
dir.out <- "data/fasta-orthomcl/orthomcl50-clust-cds-fasta-gv35-bif20"

dna.gv <- readDNAStringSet("data/raw/cds_dna_short_id_gv35.fasta", format="fasta")
dna.bif <- readDNAStringSet("data/raw/cds_dna_short_id_bif20.fasta", format="fasta")

aa.gv <- readAAStringSet("data/raw/cds_aa_short_id_gv35.fasta", format="fasta")
aa.bif <- readAAStringSet("data/raw/cds_aa_short_id_bif20.fasta", format="fasta")

dna.fasta <- c(dna.gv, dna.bif)
aa.fasta <- c(aa.gv, aa.bif)

## run orthomcl_to_fasta.R
source("code/orthomcl_to_fasta.R")

## No RData will be written, only fasta files
## I moved the resulting fasta files to my local directory due to file size. They can be replicated using this script from the main directory.