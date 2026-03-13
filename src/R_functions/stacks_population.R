library(tidyverse)

read_sumstats <- function(sumstats_file){
  sumstats <- read.table(
    sumstats_file,
    header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 38, nrows = 36
  )
  colnames(sumstats) <- c(
    "Pop_ID", "Private", "Sites", "Variant_Sites", 
    "Polymorphic_Sites", "Percent_Polymorphic_Loci", 
    "Num_Indv", "Num_Indv_Var", "Num_Indv_StdErr", 
    "P", "P_Var", "P_StdErr", 
    "Obs_Het", "Obs_Het_Var", "Obs_Het_StdErr", 
    "Obs_Hom", "Obs_Hom_Var", "Obs_Hom_StdErr", 
    "Exp_Het", "Exp_Het_Var", "Exp_Het_StdErr", 
    "Exp_Hom", "Exp_Hom_Var", "Exp_Hom_StdErr", 
    "Pi", "Pi_Var", "Pi_StdErr", 
    "Fis", "Fis_Var", "Fis_StdErr"
  )
  return(sumstats)
}

read_fst_data <- function(fst_table){
  
  fst_data <- read.table(
    fst_table,
    header=FALSE, col.names=c("pair", "fst")
  )
  
  fst_matrix <- matrix(NA, nrow=length(unique(unlist(strsplit(fst_data$pair, "_")))), ncol=length(unique(unlist(strsplit(fst_data$pair, "_")))))
  
  rownames(fst_matrix) <- colnames(fst_matrix) <- unique(unlist(strsplit(fst_data$pair, "_")))
  
  for (i in 1:nrow(fst_data)) {
    pops <- unlist(strsplit(fst_data$pair[i], "_"))
    fst_matrix[pops[1], pops[2]] <- fst_data$fst[i]
    fst_matrix[pops[2], pops[1]] <- fst_data$fst[i]
  }
  
  return(fst_matrix)
}

get_mean_fst <- function(pop_id){
  return(mean(na.omit(fst_matrix[pop_id,])))
}

get_ex_q <- function(k){
  if (k == 8){
    ex_q <- "Q6"
  }
  return(ex_q)
}

get_la_q <- function(k){
  if (k == 8){
    la_q <- "Q2"
  }
  return(la_q)
}

get_ex_la_props <- function(dataname, k){
  famfile <- paste0("data/vcfs/", dataname, ".fam")
  meanqfile <- paste0("data/admixture/", dataname,".admixture.", k, ".Q")
  fam <- read.table(famfile)
  samples <- as.character(fam[, 1])
  
  qvals <- read.table(meanqfile)
  colnames(qvals) <- paste0("Q", seq(1:ncol(qvals)))
  qvals$IID <- samples
  
  qvals_id <- data.frame(qvals[c("IID", get_ex_q(k), get_la_q(k))])
  colnames(qvals_id) <- c("IID", "EX_proportion", "LA_proportion")
  
  qvals_pop <- qvals_id |> 
    mutate(Pop_ID = substr(IID, 1, 2)) |> 
    group_by(Pop_ID) |>
    summarise(
      mean_EX = mean(EX_proportion, na.rm = TRUE),
      mean_LA = mean(LA_proportion, na.rm = TRUE)
    )
  return(qvals_pop)
}

### args

fst_table <- "data/fst/all_pairs.weighted.fst"
sumstats_file <- "data/stacks_populations/pclarkii.qc_ac_bial.lowmiss_maf.p.sumstats_summary.tsv"
dataname <- "pclarkii.qc_ac_bial.lowmiss_maf"
k <- 8

### analysis

fst_matrix <- read_fst_data(fst_table)
sumstats <- read_sumstats(sumstats_file)

sumstats$AvgFst <- sapply(sumstats$Pop_ID, get_mean_fst)

sumstats <- left_join(
  sumstats,
  get_ex_la_props(dataname, k),
  by = "Pop_ID"
  )

pihetfstqprop <- sumstats[, c("Pop_ID", "Pi", "Obs_Het", "AvgFst", "mean_EX", "mean_LA", "Fis")]

ggplot(data = pihetfstqprop |> filter(Pop_ID != "LA" & Pop_ID != "EX"),
       aes(x = mean_EX, y = Pi)) +
  #geom_smooth(method = "gam",
  #            method.args = list(family = mgcv::betar(link = "logit"))) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_log10() +
  xlab("Extremadura ancestry") + ylab("Nucleotide diversity π") +
  theme_minimal()

ggplot(data = pihetfstqprop |> filter(Pop_ID != "LA" & Pop_ID != "EX"),
       aes(x = mean_EX, y = AvgFst)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_log10() +
  xlab("Extremadura ancestry") + ylab("Average Fst") +
  theme_minimal()

ggplot(data = pihetfstqprop |> filter(Pop_ID != "LA" & Pop_ID != "EX"),
       aes(x = mean_LA, y = AvgFst)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_log10() +
  xlab("Lousiana ancestry") + ylab("Average Fst") +
  theme_minimal()

ggplot(data = pihetfstqprop |> filter(Pop_ID != "LA" & Pop_ID != "EX"),
       aes(x = mean_LA, y = Pi)) +
  geom_smooth(method = "lm") +
  geom_label(aes(label = Pop_ID), size = 2, label.padding = unit(0.01, "lines"), linewidth = 0.01) +
  scale_x_log10() +
  xlab("Lousiana ancestry") + ylab("Nucleotide diversity π") +
  theme_minimal()

ggplot(data = pihetfstqprop |> filter(Pop_ID != "LA" & Pop_ID != "EX"),
       aes(x = Fis, y = Obs_Het)) +
  geom_smooth(method = "lm") +
  geom_label(aes(label = Pop_ID), size = 2, label.padding = unit(0.01, "lines"), linewidth = 0.01) +
  #scale_x_log10() +
  #xlab("Lousiana ancestry") + ylab("Nucleotide diversity π") +
  theme_minimal()

write.table(fst_matrix, "data/fst_matrix.txt", quote = F)
write.table(pihetfstqprop, "data/pi_het_fst_qprop.txt", quote = F, row.names = F)

a <- pihetfstqprop[which(pihetfstqprop$mean_EX < 0.01 & pihetfstqprop$mean_LA < 0.001 & pihetfstqprop$Pop_ID != "EX" & pihetfstqprop$Pop_ID != "LA"), ]
b <- pihetfstqprop[which(pihetfstqprop$mean_EX > 0.01 & pihetfstqprop$Pop_ID != "EX" & pihetfstqprop$Pop_ID != "LA"), ]
c <- pihetfstqprop[which(pihetfstqprop$mean_EX < 0.01 & pihetfstqprop$mean_LA > 0.001 & pihetfstqprop$Pop_ID != "EX" & pihetfstqprop$Pop_ID != "LA"), ]


ggsave(plot = ggplot() +
  geom_boxplot(data = b, aes(y = Obs_Het, x = 1), alpha = 0.5, color = "darkgreen", outliers = F) +
  geom_boxplot(data = a, aes(y = Obs_Het, x = 2), alpha = 0.5, color = "blue", outliers = F) +
  geom_boxplot(data = c, aes(y = Obs_Het, x = 3), alpha = 0.5, color = "red", outliers = F) +
  geom_jitter(data = b, aes(y = Obs_Het, x = 1), alpha = 0.5, color = "darkgreen") +
  geom_jitter(data = a, aes(y = Obs_Het, x = 2), alpha = 0.5, color = "blue") +
  geom_jitter(data = c, aes(y = Obs_Het, x = 3), alpha = 0.5, color = "red") +
  theme_minimal() + xlab("EX mixed (green) - Isolated (blue) - Connected (red)") + ylab("Heterozygosity"),
  file = "plots/stacks_population/heterozygosity.png"
)

ggsave(plot = ggplot() +
         geom_boxplot(data = b, aes(y = Pi, x = 1), alpha = 0.5, color = "darkgreen", outliers = F) +
         geom_boxplot(data = a, aes(y = Pi, x = 2), alpha = 0.5, color = "blue", outliers = F) +
         geom_boxplot(data = c, aes(y = Pi, x = 3), alpha = 0.5, color = "red", outliers = F) +
         geom_jitter(data = b, aes(y = Pi, x = 1), alpha = 0.5, color = "darkgreen") +
         geom_jitter(data = a, aes(y = Pi, x = 2), alpha = 0.5, color = "blue") +
         geom_jitter(data = c, aes(y = Pi, x = 3), alpha = 0.5, color = "red") +
         theme_minimal() + xlab("EX mixed (green) - Isolated (blue) - Connected (red)") + ylab("Nucleotide Diversity π"),
       file = "plots/stacks_population/pi.png"
)

ggsave(plot = ggplot() +
         geom_boxplot(data = b, aes(y = AvgFst, x = 1), alpha = 0.5, color = "darkgreen", outliers = F) +
         geom_boxplot(data = a, aes(y = AvgFst, x = 2), alpha = 0.5, color = "blue", outliers = F) +
         geom_boxplot(data = c, aes(y = AvgFst, x = 3), alpha = 0.5, color = "red", outliers = F) +
         geom_jitter(data = b, aes(y = AvgFst, x = 1), alpha = 0.5, color = "darkgreen") +
         geom_jitter(data = a, aes(y = AvgFst, x = 2), alpha = 0.5, color = "blue") +
         geom_jitter(data = c, aes(y = AvgFst, x = 3), alpha = 0.5, color = "red") +
         theme_minimal() + xlab("EX mixed (green) - Isolated (blue) - Connected (red)") + ylab("Average Fst"),
       file = "plots/stacks_population/fst.png"
)

# EX ancestry has strong effect on both Pi and Het
# Populations with EX ancestry:
# "AR", "CM", "DA", "FA", "FM", "GA", "HM", "JA"
# Populations with no EX ancestry can be divided using LA ancestry as a proxy for isolation:
# isolated (a; ~ 0 LA ancestry):
# "AM", "BA", "BM", "CA", "DM", "HA", "IB", "JM"
# and connected (c; some LA ancestry)
# "AA", "AB", "BB", "CB", "CQ", "CS", "DB", "EA", "EB", "EC", "EM", "FB", "GB", "GM", "HB", "IA", "IM", "JB"
# Isolated populations have lower Pi, higher Fst but not lower Het, meaning that:
# rare variants were lost, they are not very connected to other populations, but continue random mating internally?
