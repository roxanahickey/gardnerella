## Get fasta sequences from OrthoMCL clusters
## Roxana Hickey
## 2015-03-12

## This script extracts fasta DNA and AA sequences for each cluster or singleton from the OrthoMCL clusters.
## It is a time-consuming script and should only be re-run if OrthoMCL clusters have changed (see below):
## OrthoMCL was performed on 35 genomes, PATRIC AA CDS sequences, 70% identity on March 2 2015. 
## This resulted in 2399 clusters and 1495 singletons.

library(Biostrings)
setwd("~/Documents/research/gardnerella/")
load("output/RData/orthomcl_summary.RData")

## read in DNA and AA fasta files (all 35 genomes combined)
dna.fasta <- readDNAStringSet("data/all_gv_orthomcl_compliant_dna.fasta", format="fasta")
aa.fasta <- readAAStringSet("data/all_gv_orthomcl_compliant_aa.fasta", format="fasta")

# create directories for output
dir.create("output/orthomcl-fasta")

dir.create("output/orthomcl-fasta/clust-fasta-dna")
dir.create("output/orthomcl-fasta/clust-fasta-aa")

dir.create("output/orthomcl-fasta/sing-fasta-dna")
dir.create("output/orthomcl-fasta/sing-fasta-aa")

dir.create("output/orthomcl-fasta/uniq-fasta-dna")
dir.create("output/orthomcl-fasta/uniq-fasta-aa")

# write cluster fasta files
for(i in 1:length(ortho.clust)){
  clust.tmp <- names(ortho.clust)[i]
  seq.tmp <- unlist(strsplit(ortho.clust[i], split=" "), use.names = FALSE)
  
  # DNA fasta
  dna.tmp <- dna.fasta[names(dna.fasta) %in% seq.tmp]
  writeXStringSet(dna.fasta, paste("output/orthomcl-fasta/clust-fasta-dna/", clust.tmp, "_dna.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  # AA fasta
  aa.tmp <- aa.fasta[names(aa.fasta) %in% seq.tmp]
  writeXStringSet(aa.fasta, paste("output/orthomcl-fasta/clust-fasta-aa/", clust.tmp, "_aa.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta") 
  
  # clean up
  rm(clust.tmp, seq.tmp, dna.tmp, aa.tmp)
}

# write singleton fasta files
for(i in 1:length(ortho.sing)){
  sing.tmp <- ortho.sing.id[i]
  seq.tmp <- ortho.sing[i]
  
  # DNA fasta
  dna.tmp <- dna.fasta[seq.tmp]
  writeXStringSet(dna.tmp, paste("output/orthomcl-fasta/sing-fasta-dna/", sing.tmp, "_dna.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")

  # AA fasta  
  aa.tmp <- aa.fasta[seq.tmp]
  writeXStringSet(aa.tmp, paste("output/orthomcl-fasta/sing-fasta-aa/", sing.tmp, "_aa.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  # clean up
    rm(sing.tmp, seq.tmp, dna.tmp, aa.tmp)
}

# write unique fasta files (per GV strain)
for(i in unique(ortho.sing.strain)){
  uniq.tmp <- ortho.sing[grep(i, ortho.sing)]

  # DNA fasta
  dna.tmp <- dna.fasta[names(dna.fasta) %in% uniq.tmp]
  writeXStringSet(dna.tmp, paste("output/orthomcl-fasta/uniq-fasta-dna/", i, "_uniq_dna.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  # AA fasta
  aa.tmp <- aa.fasta[names(aa.fasta) %in% uniq.tmp]
  writeXStringSet(aa.tmp, paste("output/orthomcl-fasta/uniq-fasta-aa/", i, "_uniq_aa.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  # clean up
  rm(uniq.tmp, dna.tmp, aa.tmp)
  }