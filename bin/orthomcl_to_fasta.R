## Get fasta sequences from OrthoMCL clusters
## Roxana Hickey
## 2015-03-12

## This script extracts fasta DNA and AA sequences for each cluster or singleton from the OrthoMCL clusters.
## It is a time-consuming script and should only be re-run if OrthoMCL clusters have changed (see below):
## OrthoMCL was performed on 35 genomes, PATRIC AA CDS sequences, 70% identity on March 2 2015. 
## This resulted in 2399 clusters and 1495 singletons.

# create directories for output
dir.create(dir.out)

dir.create(paste(dir.out, "/clust-dna", sep=""))
dir.create(paste(dir.out, "/clust-aa", sep=""))

dir.create(paste(dir.out, "/single-dna", sep=""))
dir.create(paste(dir.out, "/single-aa", sep=""))

# write cluster fasta files
for(i in 1:length(ortho.clust$ortho_id)){
  clust.tmp <- names(ortho.clust$ortho_id)[i]
  seq.tmp <- unlist(strsplit(ortho.clust$ortho_id[i], split=" "), use.names = FALSE)
  
  # DNA fasta
  dna.tmp <- dna.fasta[names(dna.fasta) %in% seq.tmp]
  writeXStringSet(dna.tmp, paste(dir.out, "/clust-dna/", clust.tmp, "_dna.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  # AA fasta
  aa.tmp <- aa.fasta[names(aa.fasta) %in% seq.tmp]
  writeXStringSet(aa.tmp, paste(dir.out, "/clust-aa/", clust.tmp, "_aa.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta") 
  
  # clean up
  rm(clust.tmp, seq.tmp, dna.tmp, aa.tmp)
}

# write singleton fasta files (per GV strain)
for(i in unique(ortho.sing$id_short)){
  seq.tmp <- ortho.sing$cds_id[ortho.sing$id_short == i]
  
  # DNA fasta
  dna.tmp <- dna.fasta[seq.tmp]
  writeXStringSet(dna.tmp, paste(dir.out, "/single-dna/", i, "_uniq_dna.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")

  # AA fasta  
  aa.tmp <- aa.fasta[seq.tmp]
  writeXStringSet(aa.tmp, paste(dir.out, "/single-aa/", i, "_uniq_aa.fasta", sep=""), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  # clean up
    rm(seq.tmp, dna.tmp, aa.tmp)
}