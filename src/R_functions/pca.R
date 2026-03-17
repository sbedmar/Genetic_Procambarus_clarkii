library(tidyverse)
library(ggrepel)
library(cowplot)

read_vec <- function(name, dataset){
  # read eigenvectors
  vec <- paste0("data/pca/", name, ".", dataset, ".eigenvec")
  eigenvecs <- read.table(vec, header = TRUE)[, -1]
  eigenvecs$pop <- paste0(sapply(strsplit(eigenvecs$IID, ""), `[`, 1),
                          sapply(strsplit(eigenvecs$IID, ""), `[`, 2))
  return(eigenvecs)
}

get_percents <- function(name, dataset){
  # read eigenvalues
  val <- paste0("data/pca/", name, ".", dataset, ".eigenval")
  eigenvals <- read.table(val, col.names = c("eigenval"))
  percents <- data.frame(
    percentage = eigenvals$eigenval/sum(eigenvals$eigenval)*100
    )
  return(percents)
}

plot_percents <- function(percents){
  # scree-plot of eigenvalues
  scree <- ggplot(data = percents) +
    geom_line(aes(x = as.numeric(row.names(percents)), y = percentage)) +
    geom_point(aes(x = as.numeric(row.names(percents)), y = percentage),
               fill = "grey", shape = 21, size = 2.5) +
    xlab(paste0("Number of Principal Component")) +
    ylab(paste0("Percentage of variance Explained")) +
    scale_x_continuous(n.breaks = 10) +
    theme_bw()
  return(scree)
}

plot_pca <- function(eigenvecs, percents, x = 1, y = 2){
  labels <- eigenvecs |>
    select(!"IID") |>
    group_by(pop) |>
    summarize(across(starts_with("PC"), mean))
  labels <- data.frame(labels)
  pca_plot <- ggplot() +
    geom_point(data = eigenvecs,
               aes(x = eigenvecs[, colnames(eigenvecs)[x+1]],
                   y = eigenvecs[, colnames(eigenvecs)[y+1]],
                   fill = pop),
               shape = 21, size = 2.5) +
    geom_label_repel(
      aes(x = labels[, colnames(eigenvecs)[x+1]],
          y = labels[, colnames(eigenvecs)[y+1]],
          label = labels$pop),
      size = 2.5, max.overlaps = 500,
      force = 1, segment.size = 0.1,
      box.padding = 0.05,
      label.padding = 0.05,
      label.size = 0.05,       # border thickness
      fill = "white",         # background color
      color = "black",       # text color
      alpha = 0.85
    ) +
    xlab(paste0(colnames(eigenvecs)[x+1], " - ", round(percents[x, 1], 2), "%")) +
    ylab(paste0(colnames(eigenvecs)[y+1], " - ", round(percents[y, 1], 2), "%")) +
    theme_bw() + theme(legend.position="none")
  return(pca_plot)
}


### prepare datasets for plot
plot_pca_admix <- function(name, dataset, x = 1, y = 2, k = 8){
  eigenvecs <- read_vec(name, dataset)
  percents <- get_percents(name, dataset)
  labels <- eigenvecs |>
    select(!"IID") |>
    group_by(pop) |>
    summarize(across(starts_with("PC"), mean))
  labels <- data.frame(labels)
  
  ### plot!
  baseplot <- ggplot() +
    geom_point(data = eigenvecs,
               aes(x = eigenvecs[, colnames(eigenvecs)[x+1]],
                   y = eigenvecs[, colnames(eigenvecs)[y+1]]), size = 0.001)
  
  for (sample in eigenvecs$IID){
    filename <- paste0("plots/admixture/", k, "/", sample, ".png")
    x1 <- eigenvecs[eigenvecs$IID == sample, ][,x+1]
    y1 <- eigenvecs[eigenvecs$IID == sample, ][,y+1]
    baseplot <- baseplot +
      draw_image(filename, x = x1, y = y1, vjust = 0.5, hjust = 0.5, scale = 0.017)
  }
  
  baseplot <- baseplot +
    geom_label_repel(
      aes(x = labels[, colnames(eigenvecs)[x+1]],
          y = labels[, colnames(eigenvecs)[y+1]],
          label = labels$pop),
      size = 2.5, max.overlaps = 500,
      force = 1, segment.size = 0.1,
      box.padding = 0.05,
      label.padding = 0.05,
      label.size = 0.05,       # border thickness
      fill = "white",         # background color
      color = "black",       # text color
      alpha = 0.85
    ) +
    xlab(paste0(colnames(eigenvecs)[x+1], " - ", round(percents[x, 1], 2), "%")) +
    ylab(paste0(colnames(eigenvecs)[y+1], " - ", round(percents[y, 1], 2), "%")) +
    theme_bw() + theme(legend.position="none")
  
  ggsave(
    plot = baseplot,
    filename = paste0("plots/pca/", name, ".", dataset, ".pc", x, "pc", y, ".k", k, ".png")
  )
  
}

### args
# name <- "pclarkii.qc_ac_bial.lowmiss_maf"
# dataset <- "ALL"
# k <- 8
# 
# plot_pca_admix(name, dataset, x = 1, y = 2, k = k)
# plot_pca_admix(name, dataset, x = 3, y = 4, k = k)
# plot_pca_admix(name, dataset, x = 5, y = 6, k = k)
# plot_pca_admix(name, dataset, x = 7, y = 8, k = k)
