library(tidyverse)
library(adegenet)
library(RColorBrewer)

name <- "pclarkii.qc_ac_bial.lowmiss_maf"
dataset <- "NOLA"

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

read_df_from_raw <- function(name, dataset){
  gtdf <- data.frame(as.matrix(read.PLINK(paste0("data/dapc/", name, ".", dataset, ".raw"))))
  gtdf_imputed <- as.data.frame(
    lapply(gtdf, function(x) {
      x[is.na(x)] <- get_mode(x)
      x
    })
  )
  rownames(gtdf_imputed) <- rownames(gtdf)
  return(gtdf_imputed)
}

gtdata <- read_df_from_raw(name, dataset)

grp <- find.clusters(gtdata)
for (n in c(1:20)){
  grp <- find.clusters(gtdata, n.pca = 100,
                       choose.n.clust = F, criterion = "diffNgroup")
  print(length(levels(grp$grp)))
}
grp <- find.clusters(gtdata, n.pca = 100,
                     choose.n.clust = F, criterion = "diffNgroup")
length(levels(grp$grp))


candidate_grps <- grp$grp

xval <- xvalDapc(gtdata, candidate_grps, n.pca.max = 100,
                 result = "groupMean",
                 n.pca = NULL, n.rep = 10, xval.plot = TRUE)

DAPC <- dapc(gtdata, candidate_grps,
             n.pca = xval$DAPC$n.pca,
             n.da = xval$DAPC$n.da)

dapc_matrix <- data.frame(DAPC$ind.coord)

myCol <- c(brewer.pal(n = 8, name = "Dark2"))

dapc_post <- data.frame(DAPC$posterior)
dapc_post <- dapc_post |> rownames_to_column("IID")
admix <- dapc_post |> gather(X, value, -IID)
admix_plot <- ggplot(admix, aes(x = IID, y = value, fill = factor(X))) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 2))
admix_plot

dapc_post$pop <- paste0(sapply(strsplit(dapc_post$IID, ""), `[`, 1),
                        sapply(strsplit(dapc_post$IID, ""), `[`, 2))

for (population in unique(dapc_post$pop)){
  P <- dapc_post |>
    filter(pop == population) |>
    summarize(across(starts_with("X"), sum))
  P <- P / sum(P)
  Ppie <- data.frame(
    group = colnames(P),
    value = as.numeric(P[1,])
  )
  ppp <- ggplot(Ppie, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+
    theme_void() + theme(legend.position="none")
  ggsave(paste0("plots/dapc/", population, ".png"), ppp)
}


########

for (k in c(2:10)){
  print(k)
  grp <- find.clusters(gtdata, n.pca = 100, n.clust = k)
  candidate_grps <- grp$grp
  
  DAPC <- dapc(gtdata, candidate_grps,
               n.pca = 40,
               n.da = 5)
  dapc_post <- data.frame(DAPC$posterior)
  dapc_post <- dapc_post |> rownames_to_column("IID")
  dapc_post$pop <- paste0(sapply(strsplit(dapc_post$IID, ""), `[`, 1),
                          sapply(strsplit(dapc_post$IID, ""), `[`, 2))
  
  for (population in unique(dapc_post$pop)){
    P <- dapc_post |>
      filter(pop == population) |>
      summarize(across(starts_with("X"), sum))
    P <- P / sum(P)
    Ppie <- data.frame(
      group = colnames(P),
      value = as.numeric(P[1,])
    )
    ppp <- ggplot(Ppie, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0)+
      theme_void() + theme(legend.position="none")
    ggsave(paste0("plots/dapc/", k, ".", population, ".png"), ppp)
  }
}
