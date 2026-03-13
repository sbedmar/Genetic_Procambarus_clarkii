Sys.setenv(PATH=paste("/Users/enrico/Documents/software/AdmixTools/src", Sys.getenv("PATH"), sep=":"))
library(admixr)
library(tidyverse)

prefix <- "data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat"
# list.files(path = dirname(prefix), pattern = basename(prefix), full.names = TRUE)

snps <- eigenstrat(prefix)

result <- qpAdm(
  target = c("AR", "CM", "DA", "FA", "FM", "GA", "HM", "JA"),
  sources = c("EX", "AA", "JB"),
  outgroups = c("LA1", "LA2", "LA3", "LA4"),
  data = snps,
  params = list(inbreed = "YES") # forced by new ADMIXTOOLS qpfstats
)

result <- qpAdm(
  target = c("HB"),
  sources = c("EX", "HA", "JB"),
  outgroups = c("LA1", "LA2", "LA3", "LA4"),
  data = snps,
  params = list(inbreed = "YES") # forced by new ADMIXTOOLS qpfstats
)

result$proportions
