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
plot_admix <- function(qvals){
  admix_plot <- ggplot(qvals, aes(x = IID, y = value, fill = factor(Q))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set1") +
    xlab("Sample") + ylab("Ancestry") +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 2))
  return(admix_plot)
}

dataname <- "pclarkii.qc_ac_bial.lowmiss_maf"
# read in samples from a fam file (from plink --make-bed)
famfile <- paste0("data/vcfs/", dataname,".fam")
samples <- get_samples_from_fam(famfile)

for (k in c(3:8)){
  meanqfile <- paste0("data/admixture/", dataname,".admixture.", k, ".Q")
  #
  ggsave(
    filename = paste0("plots/admixture/", dataname,".admixture.", k, ".png"),
    plot = plot_admix(read_qvals(meanQ_file = meanqfile, samples = samples)),
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
    
    PPP <- ggplot(PP, aes(x = "", y = value, fill = group)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) +
      scale_fill_brewer(palette = "Set1") +
      theme_void() + theme(legend.position="none")
    
    ggsave(
      filename = paste0("plots/admixture/", k, "/", population, ".png"),
      plot = PPP,
      height = 3, width = 3
    )
    
  }

}

#######
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

# 1) get Spain admin-1 (states/regions) and select Andalucía
spain_states <- ne_states(country = "spain", returnclass = "sf")
andalucia <- spain_states  |> 
  filter(name == "Sevilla" | name == "Huelva" | name == "Badajoz" | name == "Córdoba" | name == "Jaén" | name == "Ciudad Real")
andalucia <- st_make_valid(andalucia)
spain_states$name
sample_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE)
sample_table_NOLA <- sample_table |> 
  filter(Site != "Louisiana_1" & Site != "Louisiana_2" & Site != "Louisiana_3" & Site != "Louisiana_4")

# 5) Simple ggplot
map <- ggplot() +
  geom_sf(data = andalucia, fill = "grey95", color = "grey20", size = 0.2) +
  #geom_sf(data = rivers_andalucia, size = 0.4, alpha = 0.9) +
  coord_sf(xlim = st_bbox(andalucia)[c("xmin","xmax")],
           ylim = st_bbox(andalucia)[c("ymin","ymax")],
           expand = FALSE) +
  geom_text(data = sample_table_NOLA,
            aes(x=LONG, y=LAT, label=CODE)) + 
  theme(panel.background = element_rect(fill = "white"), #, color = "grey20"
          panel.grid = element_line(color = "grey90", linewidth = 0.3)) + 
  labs(y = NULL, x = NULL)
map 
# ggsave(paste0("plots/admixture/map.png"), map, height = 12, width = 12)  
for (k in c(3:8)){
  print(k)
  map2 <- map
  for (pop in unique(sample_table_NOLA$CODE)){
    filename <- paste0("plots/admixture/", k, "/", pop, ".png")
    popx <- sample_table_NOLA[sample_table_NOLA$CODE == pop, ]$LONG
    popy <- sample_table_NOLA[sample_table_NOLA$CODE == pop, ]$LAT
    # img <- magick::image_read(filename)
    map2 <- map2 + 
      draw_image(filename, x = popx, y = popy, vjust = 0.5, hjust = 0.5, scale = 0.3) 
  }
  # map2
  ggsave(paste0("plots/admixture/map_k", k, ".png"), map2, height = 8, width = 12)  
}

