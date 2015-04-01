## Run orthomcl_summary.R for Gardnerella + Bifidobacterium OrthoMCL
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-20

## set up custom inputs and outputs
dir.out <- "output/orthomcl-summary-gv-bif"
dir.fig <- "figs/orthomcl-summary-gv-bif"

## read in ID key containing short and full strain names
genome.key <- read.table("data/key_gv_bif_id.txt", sep="\t", header=T)

## read in input files
ortho.clust.in <- readLines("data/orthomcl_gv_bif_groups_70.txt")
ortho.sing.in <- readLines("data/orthomcl_gv_bif_single_70.txt")
clust.mx.in <- read.table("data/orthomcl_gv_bif_groups_70_mtx.txt", header=T)

## pick custom colors
col.cust <- wes_palette("Darjeeling", 5, type="discrete")[1:5]
col.cust.hc <- col.cust
col.bw.seq <- "YlGnBu"

## run orthomcl_summary.R
source("scripts/orthomcl_summary.R")

## RData will be written to output/orthomcl-summary-gv-bif/ortho_summary.RData