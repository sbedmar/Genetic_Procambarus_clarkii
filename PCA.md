# PCA

choose bcf (see [snp_filtering.md](./snp_filtering.md)):
```
# all samples filtered bcf
bcf=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf
name=$(basename -s .bcf ${bcf})
```

prune snps:
```
plink \
 --bcf ${bcf} \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.5 \
 --out data/pca/${name}.ALL

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA" data/pca/pclarkii.qc_ac_bial.lowmiss_maf.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.5 \
 --out data/pca/${name}.NOLA

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA|EX" data/pca/pclarkii.qc_ac_bial.lowmiss_maf.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.5 \
 --out data/pca/${name}.NOLAEX
```

run PCA on the three datasets
```
plink \
 --bcf ${bcf} \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.ALL.prune.out \
 --pca 'header' \
 --out data/pca/${name}.ALL

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA" data/pca/${name}.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.NOLA.prune.out \
 --pca 'header' \
 --out data/pca/${name}.NOLA

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA|EX" data/pca/${name}.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.NOLAEX.prune.out \
 --pca 'header' \
 --out data/pca/${name}.NOLAEX
```
plot in R:
```{r}
source("src/R_functions/pca.R", print.eval=TRUE)

name <- "pclarkii.qc_ac_bial.lowmiss_maf"

for (dataset in c("ALL", "NOLA", "NOLAEX")){
  # read in data
  eigenvecs <- read_vec(name, dataset)
  percents <- get_percents(name, dataset)
  # scree-plot of eigenvalues
  scree <- plot_percents(percents)
  ggsave(paste0("plots/pca/", name, ".", dataset, ".scree.png"),
         plot_percents(percents), height = 4, width = 6)
  # plot pca
  ggsave(paste0("plots/pca/", name, ".", dataset, ".pc1pc2.png"),
         plot_pca(eigenvecs, percents, x = 1, y = 2),
         height = 3, width = 4.5)
  ggsave(paste0("plots/pca/", name, ".", dataset, ".pc3pc4.png"),
         plot_pca(eigenvecs, percents, x = 3, y = 4),
         height = 3, width = 4.5)
  ggsave(paste0("plots/pca/", name, ".", dataset, ".pc5pc6.png"),
         plot_pca(eigenvecs, percents, x = 5, y = 6),
         height = 3, width = 4.5)
  ggsave(paste0("plots/pca/", name, ".", dataset, ".pc7pc8.png"),
         plot_pca(eigenvecs, percents, x = 7, y = 8),
         height = 3, width = 4.5)
}
```