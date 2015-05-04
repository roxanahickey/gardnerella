## Limbo for code that I may or may not need
## Last update: 2015-04-02

##################################################
## orthomcl_summary.R

# write lists of fasta IDs per cluster
for(i in unique(names(ortho.clust$ortho_id))){
  clust.tmp <- data.frame(strsplit(ortho.clust$ortho_id[i], " "))
  write.table(clust.tmp, paste(dir.out, "/clust-fasta-ids/", i, "_fasta_ids.txt", sep=""), row.names=F, col.names=F, quote=F)
  
  feat.tmp <- key.feat.cds[key.feat.cds$cds_id %in% unlist(clust.tmp),]
  if(nrow(feat.tmp) > 0) {write.table(feat.tmp, paste(dir.out, "/clust-anno/feat/", i, "_feat.txt", sep=""), row.names=F, quote=F)}
  
  # in newer version, replaced this with sapply and function get_anno
  ec.tmp <- key.anno$ec[key.anno$ec$cds_id %in% unlist(clust.tmp),]
  if(nrow(ec.tmp) > 0) {write.table(ec.tmp, paste(dir.out, "/clust-anno/ec/", i, "_ec.txt", sep=""), row.names=F, quote=F)}
  
  figfam.tmp <- key.anno$figfam[key.anno$figfam$cds_id %in% unlist(clust.tmp),]
  if(nrow(figfam.tmp) > 0) {write.table(figfam.tmp, paste(dir.out, "/clust-anno/figfam/", i, "_figfam.txt", sep=""), row.names=F, quote=F)}
  
  go.tmp <- key.anno$go[key.anno$go$cds_id %in% unlist(clust.tmp),]
  if(nrow(go.tmp) > 0) {write.table(go.tmp, paste(dir.out, "/clust-anno/go/", i, "_go.txt", sep=""), row.names=F, quote=F)}
  
  path.tmp <- key.anno$path[key.anno$path$cds_id %in% unlist(clust.tmp),]
  if(nrow(path.tmp) > 0) {write.table(path.tmp, paste(dir.out, "/clust-anno/path/", i, "_path.txt", sep=""), row.names=F, quote=F)}
}

# write lists of singleton fasta IDs per strain
for(i in unique(ortho.sing$id_short)){
  sing.tmp <- ortho.sing$cds_id[grep(i, ortho.sing$cds_id)]
  #   assign(paste(i, ".uniq", sep=""), sing.tmp)
  
  strain.uniq.df[i, "uniq_clust_count"] <- length(sing.tmp)
  
  write.table(sing.tmp, paste(dir.out, "/single-fasta-ids/", i, "_uniq_fasta_ids.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
}

for(i in unique(clade.list)){
  
  id.tmp <- unlist(clade.id[i], use.names=F)
  
  # reduce clust.mx to clade members
  clade.clust.tmp <- clust.mx[id.tmp,]
  clade.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)>=1]  
  
  clade.clust[i] <- list(colnames(clade.clust.tmp))
  
  write.table(colnames(clade.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # extract clade core clusters
  clade.core.clust.tmp <- clade.clust.tmp[,colSums(clade.clust.tmp)==nrow(clade.clust.tmp)]
  
  clade.core.clust[i] <- list(colnames(clade.core.clust.tmp))
  
  write.table(colnames(clade.core.clust.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # reduce clade core to those not in overall core
  clade.core.clust.red.tmp <- clade.core.clust.tmp[,!(colnames(clade.core.clust.tmp) %in% core.clust$all)]
  
  clade.core.clust.red[i] <- list(colnames(clade.core.clust.red.tmp))
  
  write.table(colnames(clade.core.clust.red.tmp), 
              file=paste(dir.out, "/clade-info/", i, "_clust_core_red.txt", sep=""), 
              row.names=F, col.names=F, quote=F)
  
  # fill in summary table
  clust.smy[i, "all_core"] <- ncol(clade.core.clust.tmp) - ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_core_red"] <- ncol(clade.core.clust.red.tmp)
  clust.smy[i, "clade_sing"] <- length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  clust.smy[i, "clade_total"] <- ncol(clade.clust.tmp) + length(ortho.sing$id_short[ortho.sing$id_short %in% id.tmp])
  
}

