library(ape)
library(pheatmap)
library(viridis)

fst_data <- read.table("data/fst/all_pairs.weighted.fst", header=FALSE, col.names=c("pair", "fst"))
fst_matrix <- matrix(0, nrow=length(unique(unlist(strsplit(fst_data$pair, "_")))), ncol=length(unique(unlist(strsplit(fst_data$pair, "_")))))
rownames(fst_matrix) <- colnames(fst_matrix) <- unique(unlist(strsplit(fst_data$pair, "_")))
for (i in 1:nrow(fst_data)) {
  pops <- unlist(strsplit(fst_data$pair[i], "_"))
  fst_matrix[pops[1], pops[2]] <- fst_data$fst[i]
  fst_matrix[pops[2], pops[1]] <- fst_data$fst[i]
}

fst_dist <- as.dist(fst_matrix)
nj_tree <- nj(fst_dist)
plot(nj_tree, "unrooted")
# root at middle point or specific outgroup
nj_tree_rooted <- root(nj_tree, outgroup = c("EX"))
pdf("plots/fst/nj_tree_rooted.pdf", height=10, width=6)
plot(nj_tree_rooted, main="Rooted NJ Tree based on Fst", use.edge.length=FALSE)
dev.off()

# draw heatmap of fst matrix with distance rows ordered by fst tree
# order matrix according to tree
is_tip <- nj_tree_rooted$edge[,2] <= length(nj_tree_rooted$tip.label)
ordered_tips <- nj_tree_rooted$edge[is_tip, 2]
ordered_pops <- nj_tree_rooted$tip.label[ordered_tips]
fst_matrix_ordered <- fst_matrix[ordered_pops, ordered_pops]
pheatmap(fst_matrix_ordered, 
         color = viridis::rocket(20),
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         display_numbers=FALSE,
         main="Fst Heatmap Ordered by NJ Tree",
         filename="plots/fst/fst_heatmap_ordered.pdf",
         width=8, height=8)


