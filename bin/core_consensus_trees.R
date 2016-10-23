## Visualize gene trees from MrBayes analysis of 664 core genes
## Roxana Hickey <roxana.hickey@gmail.com>
## Last updated: 2015-05-21

setwd("~/Documents/research/gardnerella/gardnerella-git/sandbox/")

path <- "../data/mrbayes-trees/consensus-trees/"

library(ape)
library(dplyr)

file.list <- dir(path)

check_clades <- function(t, clades){
  tmp <- lapply(clades, function(x) is.monophyletic(t, x))
  do.call(cbind.data.frame, tmp)
}

load("~/Documents/research/gardnerella/gardnerella-git/results/orthomcl70-gv35/orthomcl_summary_part2.RData")

group.id.2 <- group.id # change to group.id
group.id.2$group2 <- c(group.id.2$group2a, group.id.2$group2b)

res <- lapply(file.list, function(x) {
  t <- read.nexus(paste0(path, x))
  check_clades(t, group.id.2)
})

clust_names <- as.character(sapply(file.list, function(x) substr(x,1,10)))

out <- do.call(rbind, res)
rownames(out) <- clust_names

apply(out, 2, function(x) length(which(x))/nrow(out)) ## boring apply method
out %>% summarise_each(funs(mean)) ## fancy dplyr method

## extra playing around with plots
plot(read.nexus(paste0(path, file.list[100])))

phy <- read.nexus(paste0(path,file.list[[500]]))

cl <- genome.key.char$group_id
names(cl) <- genome.key.char$id_short

cv <- cl[phy$tip.label]
plot.phylo(phy,tip.color=col.cust[factor(cv)])

## finding close trees
phy <- lapply(file.list, function(x) read.nexus(paste0(path,x)))
class(phy) <- "multiPhylo"

t <- phy[[556]]
t <- multi2di(t)
rf <- sapply(phy, function(x) RF.dist(tree1 = t, multi2di(x)))

ids <- vector()
for (i in 1:length(names(clade.id.2))){
  ids <- c(ids, rep(names(clade.id.2[[i]]), length(clade.id.2[[i]])))
}

## BUCKy concordance tree
phy <- read.nexus("../results/bucky-orthomcl70-gv35-core-singlecopy-664/primary_conc.tre")
plot.phylo(phy,tip.color=col.cust[factor(cv)])
