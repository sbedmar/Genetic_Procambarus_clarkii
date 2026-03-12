library(tidyverse)

dist.m <- read.table("data/fastreer/pclarkii.qc_ac_bial.lowmiss_maf.dist", skip = 1) |> 
  column_to_rownames("V1")
colnames(dist.m) <- rownames(dist.m)
dist.d <- as.dist(as.matrix(dist.m))

# ape
library(ape)
nj.tree <- ape::nj(dist.d)
# remove EX03 and LA107 from the tree
nj.tree <- ape::drop.tip(nj.tree, c("EX03", "LA107", "LA305"))
# find what edge contains all samples for rooting
root.tips <- grep("^LA3|^LA4|^LA1", nj.tree$tip.label, value = TRUE)
root.mrca <- ape::getMRCA(nj.tree, extremadura.tips)
# reroot at root.mrca 
nj.tree <- ape::root(nj.tree, node = root.mrca, resolve.root = TRUE)

# plot circular tree with constant branch lengths
pdf("plots/fastreer/nj.tree.pdf", height = 20, width = 20)
plot.phylo(nj.tree, type = "fan", cex = 0.5, no.margin = TRUE, use.edge.length = TRUE)
dev.off()
# plot tree with constant branch lengths
pdf("plots/fastreer/nj.tree.pdf", height = 30, width = 10)
plot.phylo(nj.tree, no.margin = TRUE, use.edge.length = FALSE)
dev.off()

# export to newick
ape::write.tree(nj.tree, "data/fastreer/pclarkii.nj.tree.newick")

# fancy tree
