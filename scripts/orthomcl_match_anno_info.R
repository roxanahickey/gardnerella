## Extracting protein info!!
group.df <- data.frame(cds_id=unlist(strsplit(ortho.clust, " ")),
                       clust_id=substr(names(unlist(strsplit(ortho.clust, " "))), 1, 10),
                       row.names=NULL)
sing.df <- data.frame(cds_id=ortho.sing,
                      clust_id=ortho.sing.id,
                      row.names=NULL)
clust.df <- rbind(group.df, sing.df)

fasta.key <- read.table("fasta_headers_key.txt", sep="\t", header=T)

fasta.superkey <- merge(clust.df, fasta.key)

## GV id key
gv.id.key <- read.table("~/Documents/research/gardnerella/data/GV_id_name_clade_key.txt", header=T, sep="\t")
gv.id.key$ortho_clade <- cutg[gv.id.key$id_short]

## all core genes (plus all-1)
# dir.create("core-all-clades")

clust.core.all.df <- fasta.superkey[fasta.superkey$clust_id %in% colnames(clust.core.all),]
write.table(clust.core.all.df, file="core-all-clades/core_all35_info.txt", sep="\t", row.names=F)
write.table(clust.core.all.df$short_id, file="core-all-clades/core_all35_fid_list.txt", row.names=F, col.names=F, quote=F)

clust.core.less1.df <- fasta.superkey[fasta.superkey$clust_id %in% colnames(clust.core.less1),]
write.table(clust.core.less1.df, file="core-all-clades/core_34_info.txt", sep="\t", row.names=F)
write.table(clust.core.less1.df$short_id, file="core-all-clades/core_34_fid_list.txt", row.names=F, col.names=F, quote=F)

pick <- colnames(clust.core.less1)[!(colnames(clust.core.less1) %in% colnames(clust.core.all))]
clust.core.less1.uniq.df <- clust.core.less1.df[clust.core.less1.df$clust_id %in% pick,]
write.table(clust.core.less1.uniq.df, file="core-all-clades/core_34_uniq_info.txt", sep="\t", row.names=F)
write.table(clust.core.less1.uniq.df$short_id, file="core-all-clades/core_34_uniq_fid_list.txt", row.names=F, col.names=F, quote=F)

## clade core unique dataframes
clade1.clust.core.uniq.df <- fasta.superkey[fasta.superkey$clust_id %in% clade1.clust.core.uniq,]
write.table(clade1.clust.core.uniq.df, file="clade-core-uniq/clade1_clus_core_uniq_info.txt", sep="\t", row.names=F)

clade2.clust.core.uniq.df <- fasta.superkey[fasta.superkey$clust_id %in% clade2.clust.core.uniq,]
write.table(clade2.clust.core.uniq.df, file="clade-core-uniq/clade2_clus_core_uniq_info.txt", sep="\t", row.names=F)

clade3.clust.core.uniq.df <- fasta.superkey[fasta.superkey$clust_id %in% clade3.clust.core.uniq,]
write.table(clade3.clust.core.uniq.df, file="clade-core-uniq/clade3_clus_core_uniq_info.txt", sep="\t", row.names=F)

clade4.clust.core.uniq.df <- fasta.superkey[fasta.superkey$clust_id %in% clade4.clust.core.uniq,]
write.table(clade4.clust.core.uniq.df, file="clade-core-uniq/clade4_clus_core_uniq_info.txt", sep="\t", row.names=F)