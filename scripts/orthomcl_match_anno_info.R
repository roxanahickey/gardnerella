## Match OrthoMCL clusters to original fasta IDs to get gene annotations
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-14

load("output/RData/orthomcl_summary.RData")

## combine clusters and singletons into one dataframe
group.df <- data.frame(cds_id=unlist(strsplit(ortho.clust, " ")),
                       ortho_id=substr(names(unlist(strsplit(ortho.clust, " "))), 1, 10),
                       row.names=NULL)
sing.df <- data.frame(cds_id=ortho.sing,
                      ortho_id=ortho.sing.id,
                      row.names=NULL)
clust.df <- rbind(group.df, sing.df)

## read in CDS key, combine with clust.df to make master key
fasta.key <- read.table("data/key_gv_patric_all35_cds_fasta.txt", sep="\t", header=T)
master.key <- merge(clust.df, fasta.key)

## read in EC, FigFam, GO and Path keys
ec.key <- read.table("data/gv_patric_all35.ec", sep="\t", header=T, quote="")
figfam.key <- read.table("data/gv_patric_all35.figfam", sep="\t", header=T, quote="")
go.key <- read.table("data/gv_patric_all35.go", sep="\t", header=T, quote="")
path.key <- read.table("data/gv_patric_all35.path", sep="\t", header=T, quote="")

# merge annotation keys
anno.key <- merge(figfam.key, ec.key, all.x=T)
anno.key <- merge(anno.key, go.key, all.x=T)
anno.key <- merge(anno.key, path.key, all.x=T)

# add short_id variable in same format as master.key (used in fasta headers)
anno.key$short_id <- paste("fid", anno.key$na_feature_id, "locus", anno.key$locus_tag, sep="|")

# match anno.key with OrthoMCL cluster IDs from master.key (~700 rows lost in process)
anno.master <- merge(anno.key, master.key[,c("short_id","ortho_id")])

# add ortho_id to individual keys
ec.key$ortho_id <- anno.master$ortho_id[match(ec.key$locus_tag, anno.master$locus_tag)]
figfam.key$ortho_id <- anno.master$ortho_id[match(figfam.key$locus_tag, anno.master$locus_tag)]
go.key$ortho_id <- anno.master$ortho_id[match(go.key$locus_tag, anno.master$locus_tag)]
path.key$ortho_id <- anno.master$ortho_id[match(path.key$locus_tag, anno.master$locus_tag)]

##################################################
## write annotation information in tables
dir.create("output/orthomcl-match-anno")

clust.core.all.df <- master.key[master.key$ortho_id %in% clust.core.all,]
write.table(clust.core.all.df, file="output/orthomcl-match-anno/core_all35_info.txt", sep="\t", row.names=F)
write.table(clust.core.all.df$short_id, file="output/orthomcl-match-anno/core_all35_fid_list.txt", 
            row.names=F, col.names=F, quote=F)

clust.core.less1.df <- master.key[master.key$ortho_id %in% clust.core.less1,]
write.table(clust.core.less1.df, file="output/orthomcl-match-anno/core_34_info.txt", sep="\t", row.names=F)
write.table(clust.core.less1.df$short_id, file="output/orthomcl-match-anno/core_34_fid_list.txt", 
            row.names=F, col.names=F, quote=F)

pick <- clust.core.less1[!(clust.core.less1 %in% clust.core.all)]
clust.core.less1.uniq.df <- clust.core.less1.df[clust.core.less1.df$ortho_id %in% pick,]
write.table(clust.core.less1.uniq.df, file="output/orthomcl-match-anno/core_34_uniq_info.txt", sep="\t", row.names=F)
write.table(clust.core.less1.uniq.df$short_id, file="output/orthomcl-match-anno/core_34_uniq_fid_list.txt", 
            row.names=F, col.names=F, quote=F)

## clade core unique dataframes (protein annotations from fasta headers)
clade1.clust.core.uniq.df <- master.key[master.key$ortho_id %in% clade1.clust.core.uniq,]
write.table(clade1.clust.core.uniq.df, file="output/orthomcl-match-anno/clade1_core_uniq_info.txt", 
            sep="\t", row.names=F)
write.table(clade1.clust.core.uniq.df$short_id, file="output/orthomcl-match-anno/clade1_core_uniq_fid.txt",
            row.names=F, col.names=F, quote=F)

clade2.clust.core.uniq.df <- master.key[master.key$ortho_id %in% clade2.clust.core.uniq,]
write.table(clade2.clust.core.uniq.df, file="output/orthomcl-match-anno/clade2_core_uniq_info.txt", 
            sep="\t", row.names=F)
write.table(clade2.clust.core.uniq.df$short_id, file="output/orthomcl-match-anno/clade2_core_uniq_fid.txt",
            row.names=F, col.names=F, quote=F)

clade3.clust.core.uniq.df <- master.key[master.key$ortho_id %in% clade3.clust.core.uniq,]
write.table(clade3.clust.core.uniq.df, file="output/orthomcl-match-anno/clade3_core_uniq_info.txt", 
            sep="\t", row.names=F)
write.table(clade3.clust.core.uniq.df$short_id, file="output/orthomcl-match-anno/clade3_core_uniq_fid.txt",
            row.names=F, col.names=F, quote=F)

clade4.clust.core.uniq.df <- master.key[master.key$ortho_id %in% clade4.clust.core.uniq,]
write.table(clade4.clust.core.uniq.df, file="output/orthomcl-match-anno/clade4_core_uniq_info.txt", 
            sep="\t", row.names=F)
write.table(clade4.clust.core.uniq.df$short_id, file="output/orthomcl-match-anno/clade4_core_uniq_fid.txt",
            row.names=F, col.names=F, quote=F)

##################################################
## write tables for GO/Path/EC/FigFam annotations
# core 35
clust.core.all.anno <- anno.master[anno.master$ortho_id %in% clust.core.all,]
write.table(clust.core.all.anno, file="output/orthomcl-match-anno/core_all35_anno.txt", sep="\t", row.names=F)

# core 34
clust.core.less1.anno <- anno.master[anno.master$ortho_id %in% clust.core.less1,]
write.table(clust.core.less1.anno, file="output/orthomcl-match-anno/core_34_anno.txt", sep="\t", row.names=F)

# no anno info for clade 1 unique (clust_1265)

# clade 2
clade2.clust.core.uniq.anno <- anno.master[anno.master$ortho_id %in% clade2.clust.core.uniq,]
write.table(clade2.clust.core.uniq.anno, file="output/orthomcl-match-anno/clade2_core_uniq_anno.txt", 
            sep="\t", row.names=F)

# clade 3
clade3.clust.core.uniq.anno <- anno.master[anno.master$ortho_id %in% clade3.clust.core.uniq,]
write.table(clade3.clust.core.uniq.anno, file="output/orthomcl-match-anno/clade3_core_uniq_anno.txt", 
            sep="\t", row.names=F)

# clade 4
clade4.clust.core.uniq.anno <- anno.master[anno.master$ortho_id %in% clade4.clust.core.uniq,]
write.table(clade4.clust.core.uniq.anno, file="output/orthomcl-match-anno/clade4_core_uniq_anno.txt", 
            sep="\t", row.names=F)

##################################################
## clean up, save image
rm(pick)
save.image("output/RData/orthomcl_match_anno.RData")
