library(tidyverse)
library(ggrepel)

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
      size = 2, max.overlaps = 500,
      box.padding = 0.1,
      label.size = 0.1,       # border thickness
      fill = "white",         # background color
      color = "black"         # text color
    ) +
    xlab(paste0(colnames(eigenvecs)[x+1], " - ", round(percents[x, 1], 2), "%")) +
    ylab(paste0(colnames(eigenvecs)[y+1], " - ", round(percents[y, 1], 2), "%")) +
    theme_bw() + theme(legend.position="none")
  return(pca_plot)
}
