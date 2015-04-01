## Run orthomcl_to_fasta.R for Gardnerella + Bifidobacterium OrthoMCL
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-25

library(Biostrings)

## load output from orthomcl_summary.R
load("output/orthomcl-summary-gv-bif/orthomcl_summary.RData")

## clear existing output paths
rm(dir.fig, dir.out)

## set up custom inputs and outputs
dir.out <- "output/orthomcl-fasta-gv-bif"

dna.gv <- readDNAStringSet("data/cds_dna_short_id_gv_all35.fasta", format="fasta")
dna.bif <- readDNAStringSet("data/cds_dna_short_id_bif_all22.fasta", format="fasta")

aa.gv <- readAAStringSet("data/cds_aa_short_id_gv_all35.fasta", format="fasta")
aa.bif <- readAAStringSet("data/cds_aa_short_id_bif_all22.fasta", format="fasta")

dna.fasta <- c(dna.gv, dna.bif)
aa.fasta <- c(aa.gv, aa.bif)

## run orthomcl_to_fasta.R
source("scripts/orthomcl_to_fasta.R")

## No RData will be written, only fasta files
## I moved the resulting fasta files to my local directory due to file size. They can be replicated using this script from the main directory.