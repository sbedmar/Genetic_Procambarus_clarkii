library(tidyverse)

dist.m <- read.table("data/fastreer/pclarkii.qc_ac_bial.lowmiss_maf.dist", skip = 1) |> 
  column_to_rownames("V1")
colnames(dist.m) <- rownames(dist.m)
dist.d <- as.dist(as.matrix(dist.m))

# ape
library(ape)
nj.tree <- ape::nj(dist.d)

pdf("plots/fastreer/nj.tree.pdf", height = 20, width = 20)
plot(nj.tree, "unrooted")
dev.off()

# export to newick
ape::write.tree(nj.tree, "data/fastreer/pclarkii.nj.tree.newick")

# fancy tree
