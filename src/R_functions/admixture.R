library(tidyverse)
library(RColorBrewer)

# get samples from fam file
get_samples_from_fam <- function(famfile){
  fam <- read.table(famfile)
  return(as.character(fam[, 1]))
}

# read a .meanQ file into a table
read_qvals <- function(meanQ_file, samples){
  qvals <- read.table(meanQ_file)
  colnames(qvals) <- paste0("Q", seq(1:ncol(qvals)))
  qvals$IID <- samples
  qvals <- qvals %>% gather(Q, value, -IID)
  return(qvals)
}

# plot admixture from df generated with read_qvals(meanQ_file, samples)
plot_admix <- function(qvals, k){
  admix_plot <- ggplot(qvals, aes(x = IID, y = value, fill = factor(Q))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colpal_f(k)) +
    xlab("Sample") + ylab("Ancestry") +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 2))
  return(admix_plot)
}

colpal_f <- function(k){
  
  colpal <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  
  if (k == 3){
    cp <- colpal[c(3, 6, 4)]
  } else if (k == 4){
    cp <- colpal[c(2, 6, 4, 3)]
  } else if (k == 5){
    cp <- colpal[c(4, 6, 3, 8, 2)]
  } else if (k == 6){
    cp <- colpal[c(3, 4, 8, 5, 6, 2)]
  } else if (k == 7){
    cp <- colpal[c(2, 5, 4, 8, 1, 6, 3)]
  } else if (k == 8){
    cp <- colpal[c(5, 8, 4, 1, 7, 3, 6, 2)]
  }
  return(cp)
}

dataname <- "pclarkii.qc_ac_bial.lowmiss_maf"
# read in samples from a fam file (from plink --make-bed)
famfile <- paste0("data/vcfs/", dataname,".fam")
samples <- get_samples_from_fam(famfile)

for (k in c(3:8)){
  print(k)
  meanqfile <- paste0("data/admixture/", dataname,".admixture.", k, ".Q")
  #
  ggsave(
    filename = paste0("plots/admixture/", dataname,".admixture.", k, ".png"),
    plot = plot_admix(read_qvals(meanQ_file = meanqfile, samples = samples), k),
    height = 5, width = 15
  )
  
  ##########
  
  admix <- read_qvals(meanQ_file = meanqfile, samples = samples)
  admix$pop <- paste0(sapply(strsplit(admix$IID, ""), `[`, 1),
                      sapply(strsplit(admix$IID, ""), `[`, 2))
  
  for (population in unique(admix$pop)){
    
    P <- admix |>
      filter(pop == population) |>
      group_by(Q) |>
      summarize(Qsum = sum(value)) |>
      ungroup() |>
      mutate(Qprop = Qsum/sum(Qsum)) 
    
    PP <- data.frame(
      group = P$Q,
      value = P$Qprop
    )

    PPP <- ggplot(PP, aes(x = 1, y = value, fill = group)) +
      geom_col(width = 0.1, color = "grey20", linewidth = 2) +
      coord_polar("y", start=0) +
      scale_fill_manual(values = colpal_f(k)) +
      #scale_color_manual(values = "grey20") +
      theme_void() + theme(legend.position="none")

    ggsave(
      filename = paste0("plots/admixture/", k, "/", population, ".png"),
      plot = PPP,
      height = 3, width = 3
    )
    
  }

}

####
####

library(sf)
library(tidyverse)
library(cowplot)
library(igraph)

sample_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE)
sample_table_NOLA <- sample_table |> 
  filter(Site != "Louisiana_1" & Site != "Louisiana_2" & Site != "Louisiana_3" & Site != "Louisiana_4")

plot_base_map <- function(){

  cuenc <- sf::read_sf("data/maps/hi_cuenca_s_ES050.shp")
  c <- st_as_sf(cuenc)
  
  guad <- sf::read_sf(paste0("data/maps/", "hi_tramocurso_s_ES050.shp"))
  g <- st_as_sf(guad)
  
  estanc <- sf::read_sf(paste0("data/maps/", "hi_aguaestanc_s_ES050.shp"))
  e <- st_as_sf(estanc)
  
  # Plot
  p <- ggplot() +
    geom_sf(data = c, fill = "grey90", color = "grey20",  linewidth = 0.2) +
    geom_sf(data = g, fill = "#70a4ba", color = "#396ecb", linewidth = 0.32) +
    geom_sf(data = e, fill = "#70a4ba", color = "#396ecb", linewidth = 0.32) +
    geom_point(data = sample_table_NOLA, aes(x=LONG, y=LAT), shape = 1, stroke = 0.4) +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey20")
    ) +
    xlab("") + ylab("")
  
  return(p)
}
# plot_size = 5
ggsave(plot = plot_base_map(), filename = "plots/admixture/basemap.png",
       width = plot_size*2, height = plot_size*1.2)

separate_close_points <- function(coords, min_dist, max_iter = 100, tol = 1e-9) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  
  if (ncol(coords) != 2) stop("coords must have 2 columns")
  
  orig <- coords
  n <- nrow(coords)
  
  for (iter in seq_len(max_iter)) {
    disp <- matrix(0, n, 2)
    any_overlap <- FALSE
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        dx <- coords[j, 1] - coords[i, 1]
        dy <- coords[j, 2] - coords[i, 2]
        d  <- sqrt(dx * dx + dy * dy)
        if (d < min_dist) {
          any_overlap <- TRUE
          if (d < tol) {
            # identical points: choose a fixed tiny direction
            ux <- 1
            uy <- 0
            d <- 0
          } else {
            ux <- dx / d
            uy <- dy / d
          }
          overlap <- min_dist - d
          
          # split correction equally between the two points
          shift_x <- 0.5 * overlap * ux
          shift_y <- 0.5 * overlap * uy
          disp[i, 1] <- disp[i, 1] - shift_x
          disp[i, 2] <- disp[i, 2] - shift_y
          disp[j, 1] <- disp[j, 1] + shift_x
          disp[j, 2] <- disp[j, 2] + shift_y
        }
      }
    }
    if (!any_overlap) break
    coords <- coords + disp
  }
  coords
}

coords <- as.matrix(data.frame(x = sample_table_NOLA$LONG, y = sample_table_NOLA$LAT))
new_coords <- separate_close_points(coords, min_dist = 0.2)
ct <- data.frame(cbind(sample_table_NOLA$CODE, new_coords))
colnames(ct) <- c("CODE", "LONG", "LAT")

plot_size <- 5.5

for (k in c(3:8)){
  print(k)
  map <- plot_base_map()
  for (pop in unique(sample_table_NOLA$CODE)){
    filename <- paste0("plots/admixture/", k, "/", pop, ".png")
    x1 <- sample_table_NOLA[sample_table_NOLA$CODE == pop, ]$LONG
    y1 <- sample_table_NOLA[sample_table_NOLA$CODE == pop, ]$LAT
    x2 <- as.numeric(ct[ct$CODE == pop, ]$LONG)
    y2 <- as.numeric(ct[ct$CODE == pop, ]$LAT)
    dxy <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
    map <- map +
      geom_segment(data = dxy, aes(x = x1, xend = x2, y = y1, yend = y2), color = "grey20",  linewidth = 0.3) +
      draw_image(filename, x = x2, y = y2, vjust = 0.5, hjust = 0.5, scale = 0.2) 
  }
  map <- map +
    geom_label(data = ct, aes(x = as.numeric(LONG), y = as.numeric(LAT), label = CODE),
               size = 2, label.padding = unit(0.1, "lines"), linewidth = 0.01)
  ggsave(paste0("plots/admixture/map_k", k, ".png"), map, width = plot_size*2, height = plot_size*1.2)  
   # ggsave(paste0("plots/admixture/map_k", k, ".pdf"), map, height = 8, width = 12)  
}

