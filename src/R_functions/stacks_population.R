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

### args

fst_table <- "data/fst/all_pairs.weighted.fst"
sumstats_file <- "data/stacks_populations/pclarkii.qc_ac_bial.lowmiss_maf.p.sumstats_summary.tsv"

### analysis

fst_matrix <- read_fst_data(fst_table)
sumstats <- read_sumstats(sumstats_file)

sumstats$AvgFst <- sapply(sumstats$Pop_ID, get_mean_fst)

ggplot() +
  geom_label(data = sumstats, aes(x = AvgFst, y = Pi, label = Pop_ID))

pihetfst <- sumstats[, c("Pop_ID", "Pi", "Obs_Het", "AvgFst")]

write.table(fst_matrix, "data/fst_matrix.txt", quote = F)
write.table(pihetfst, "data/pi_het_fst.txt", quote = F, row.names = F)

