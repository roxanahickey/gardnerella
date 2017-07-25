get.mode <- function(x){
  y <- which.max(tabulate(as.factor(x)))
  x[y]
}

ortho.id <- sort(unique(key.feat.cds$ortho_id))

key.ortho.feat <- data.frame(ortho_id = character(),
                             top_hit = character(),
                             stringsAsFactors = F)

for(i in 1:length(ortho.id)){
  key.ortho.feat[i,"ortho_id"] <- ortho.id[i]
  tmp <- key.feat.cds$product[key.feat.cds$ortho_id == ortho.id[i]]
  key.ortho.feat$top_hit[i] <- get.mode(as.character(tmp))
}

write.csv(key.ortho.feat, "data/reference/ortho_id_match_annotation.csv")

ab.cd.ef.res <- read.csv("_sandbox/functional-anno-gv/20150517-gsea-HS/table_enrich_raxml_clades_ab_cd_ef_2016_qvalues.csv")
ab.cd.ef.ortho <- subset(ab.cd.ef.res, Database == "ortho_id")

colnames(ab.cd.ef.ortho) <- c("Group", "Database", "ortho_id", "a", "b", "c", "d", 
                              "Odds.ratio", "p.value", "q.value")

ab.cd.ef.ortho.mrg <- merge(ab.cd.ef.ortho, key.ortho.feat)

write.csv(ab.cd.ef.ortho.mrg, "_sandbox/functional-anno-gv/20150517-gsea-HS/table_enrich_raxml_clades_ab_cd_ef_orthoID_anno_2016_qvalues.csv",
          row.names = FALSE)
