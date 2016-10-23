## Run orthomcl_to_fasta.R for Gardnerella OrthoMCL
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-25

library(Biostrings)

## load output from orthomcl_summary.R
load("results/orthomcl-summary-gv/orthomcl_summary.RData")

## clear existing output paths
rm(dir.fig, dir.out)

## set up custom inputs and outputs
dir.out <- "data/fasta-orthomcl/clust-cds-fasta-gv"

dna.fasta <- readDNAStringSet("data/raw/cds_dna_short_id_gv_all35.fasta", format="fasta")

aa.fasta <- readAAStringSet("data/raw/cds_aa_short_id_gv_all35.fasta", format="fasta")

## run orthomcl_to_fasta.R
source("code/orthomcl_to_fasta.R")

## No RData will be written, only fasta files
## I moved the resulting fasta files to my local directory due to file size. They can be replicated using this script from the main directory.