## GO term enrichment analysis
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-18

load("output/RData/orthomcl_match_anno.RData")
library(reshape)

dir.create("output/orthomcl-gene-enrich")


## set up a master cont tables (start with all zeros)
ortho.id.all <- c(names(ortho.clust), ortho.sing.id)

ortho.cont <- data.frame(matrix(0, nrow=length(ortho.id.all), ncol=8))
rownames(ortho.cont) <- ortho.id.all
colnames(ortho.cont) <- c("clade1.pres", "clade1.abs", "clade2.pres", "clade2.abs",
                          "clade3.pres", "clade3.abs", "clade4.pres", "clade4.abs")

# fill in counts for each clade
ortho.cont[colnames(clade1.clust), "clade1.pres"] <- colSums(clade1.clust)
ortho.cont[colnames(clade1.clust), "clade1.abs"] <- nrow(clade1.clust) - colSums(clade1.clust)

ortho.cont[colnames(clade2.clust), "clade2.pres"] <- colSums(clade2.clust)
ortho.cont[colnames(clade2.clust), "clade2.abs"] <- nrow(clade2.clust) - colSums(clade2.clust)

ortho.cont[colnames(clade3.clust), "clade3.pres"] <- colSums(clade3.clust)
ortho.cont[colnames(clade3.clust), "clade3.abs"] <- nrow(clade3.clust) - colSums(clade3.clust)

ortho.cont[colnames(clade4.clust), "clade4.pres"] <- colSums(clade4.clust)
ortho.cont[colnames(clade4.clust), "clade4.abs"] <- nrow(clade4.clust) - colSums(clade4.clust)

write.table(ortho.cont, "output/orthomcl-gene-enrich/ortho_clust_contingency.txt", sep="\t", row.names=F, quote=F)

go.cont <- match(rownames(ortho.cont), anno.master)


##### back burner
# make separate cont table for each clade vs. all others
ortho.cont.1 <- data.frame(clade1.pres = ortho.cont$clade1.pres,
                           clade1.abs = ortho.cont$clade1.abs,
                           other.pres = rowSums(ortho.cont[,c("clade2.pres", "clade3.pres", "clade4.pres")]),
                           other.abs = rowSums(ortho.cont[,c("clade2.abs", "clade3.abs", "clade4.abs")]))

ortho.cont.2 <- data.frame(clade2.pres = ortho.cont$clade2.pres,
                           clade2.abs = ortho.cont$clade2.abs,
                           other.pres = rowSums(ortho.cont[,c("clade1.pres", "clade3.pres", "clade4.pres")]),
                           other.abs = rowSums(ortho.cont[,c("clade1.abs", "clade3.abs", "clade4.abs")]))

ortho.cont.3 <- data.frame(clade3.pres = ortho.cont$clade3.pres,
                           clade3.abs = ortho.cont$clade3.abs,
                           other.pres = rowSums(ortho.cont[,c("clade2.pres", "clade1.pres", "clade4.pres")]),
                           other.abs = rowSums(ortho.cont[,c("clade2.abs", "clade1.abs", "clade4.abs")]))

ortho.cont.4 <- data.frame(clade4.pres = ortho.cont$clade4.pres,
                           clade4.abs = ortho.cont$clade4.abs,
                           other.pres = rowSums(ortho.cont[,c("clade2.pres", "clade3.pres", "clade1.pres")]),
                           other.abs = rowSums(ortho.cont[,c("clade2.abs", "clade3.abs", "clade1.abs")]))

## calculate odds ratios
or.1 <- (ortho.cont.1$clade1.pres * ortho.cont.1$other.abs) / (ortho.cont.1$clade1.abs * ortho.cont.1$other.pres)
names(or.1) <- rownames(ortho.cont)

or.2 <- (ortho.cont.2$clade2.pres * ortho.cont.2$other.abs) / (ortho.cont.2$clade2.abs * ortho.cont.2$other.pres)
names(or.2) <- rownames(ortho.cont)

or.3 <- (ortho.cont.3$clade3.pres * ortho.cont.3$other.abs) / (ortho.cont.3$clade3.abs * ortho.cont.3$other.pres)
names(or.3) <- rownames(ortho.cont)

or.4 <- (ortho.cont.4$clade4.pres * ortho.cont.4$other.abs) / (ortho.cont.4$clade4.abs * ortho.cont.4$other.pres)
names(or.4) <- rownames(ortho.cont)

or.1.over <- na.exclude(or.1[or.1>1])
or.1.under <- na.exclude(or.1[or.1<1])

or.2.over <- na.exclude(or.2[or.2>1])
or.2.under <- na.exclude(or.2[or.2<1])

or.3.over <- na.exclude(or.3[or.3>1])
or.3.under <- na.exclude(or.3[or.3<1])

or.4.over <- na.exclude(or.4[or.4>1])
or.4.under <- na.exclude(or.4[or.4<1])


##################################################
# # get lists of clade cluster IDs
# ortho.clade.df <- master.key[,c("cds_id", "ortho_id")]
# ortho.clade.df$id_short <- substr(ortho.clade.df$cds_id, 1, 4)
# ortho.clade.df <- merge(ortho.clade.df, gv.id[,c("id_short", "ortho_clade")])
# ortho.clade.df$count <- rep(1, nrow(ortho.clade.df))
# 
# # recast into count table
# ortho.clade.count <- cast(ortho.clade.df, ortho_id ~ ortho_clade, sum)
# colnames(ortho.clade.count) <- c("ortho_id", "clade1", "clade2", "clade3", "clade4")
# 
# ortho.clade.count.nosing <- subset(ortho.clade.count, ortho_id %in% names(ortho.clust))
# 
# # scale counts by number of genomes in each clade
# ortho.clade.count.scaled <- ortho.clade.count
# ortho.clade.count.scaled$clade1 <- ortho.clade.count.scaled$clade1 / length(clade1.id)
# ortho.clade.count.scaled$clade2 <- ortho.clade.count.scaled$clade2 / length(clade2.id)
# ortho.clade.count.scaled$clade3 <- ortho.clade.count.scaled$clade3 / length(clade3.id)
# ortho.clade.count.scaled$clade4 <- ortho.clade.count.scaled$clade4 / length(clade4.id)
# 
# write.table(ortho.clade.count.scaled, 
#             "output/orthomcl-gene-enrich/orthomcl_clade_cluster_count_scaled.txt", 
#             sep="\t", row.names=F, quote=F)
# 
# ortho.clade.count.nosing.scaled <- subset(ortho.clade.count.scaled, ortho_id %in% names(ortho.clust))
