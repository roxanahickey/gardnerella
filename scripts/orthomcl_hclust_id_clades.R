## Perform hierarchical clustering of genomes based on gene cluster presence/absence
## Roxana Hickey <roxana.hickey@gmail.com>
## 2015-03-14

## This script is sourced from within orthomcl_summarize.R. If running alone, first set appropriate working directory and load Rdata.

library(cluster)
library(gclus)
source("scripts/hcoplot.R")

# hclust based on gene clusters
hc.a <- hclust(dist(clust.mx), method="average")
plot(hc.a, main="Average (UPGMA)")

# average silhouette width
asw <- numeric(nrow(clust.mx))

for (k in 2:(nrow(clust.mx)-1)) {
  sil <- silhouette(cutree(hc.a, k=k), dist(clust.mx))
  asw[k] <- summary(sil)$avg.width
}

rm(k)

k.best <- which.max(asw)

# plot silhouette widths indicating optimal cluster size
png(paste(dir.fig, "/orthomcl_hclust_sil_avg.png", sep=""), width=6, height=4, units="in", res=300, pointsize=8)
plot(1:nrow(clust.mx), asw, type="h", 
     main="Silhouette-optimal number of clusters, UPGMA",
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2, col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
dev.off()

# cut the tree and assign samples to each of the groups
cutg <- cutree(hc.a, k=k.best)
sil <- silhouette(cutg, dist(clust.mx))
rownames(sil) <- row.names(clust.mx)

# plot silhouette partition
png(paste(dir.fig, "/orthomcl_hclust_sil_part.png", sep=""), width=6, height=4, units="in", res=300, pointsize=8)
plot(sil, main="Silhouette plot - UPGMA", cex.names=0.5, col=col.cust.hc, nmax=100)
dev.off()

# plot dendrogram with group labels
png(paste(dir.fig, "/orthomcl_hclust_dendro.png", sep=""), width=8, height=6, units="in", res=300, pointsize=8)
hcoplot(hc.a, dist(clust.mx), k=k.best)
dev.off()

# add original strain names
hc.a.2 <- hc.a
hc.a.2$labels <- genome.key$id_abbr[order(genome.key$id_short)]

# plot dendrogram with group labels + original strain names
png(paste(dir.fig, "/orthomcl_hclust_dendro_genome_ids.png", sep=""), width=8, height=6, units="in", res=300, pointsize=8)
hcoplot(hc.a.2, dist(clust.mx), k=k.best)
dev.off()

# clean up, save image
save.image(paste(dir.out, "/orthomcl_hclust_id_clades.RData", sep=""))
rm(asw, k.best, sil)