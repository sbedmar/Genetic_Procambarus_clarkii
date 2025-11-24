library(tidyverse)

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
    xlab("Sample") + ylab("Ancestry") +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 2))
  return(admix_plot)
}


