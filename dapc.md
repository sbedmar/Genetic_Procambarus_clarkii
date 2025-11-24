# dapc

choose bcf (see [snp_filtering.md](./snp_filtering.md)):
```
# all samples filtered bcf
bcf=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf
name=$(basename -s .bcf ${bcf})
```

get RAW from pruned snps from [pca](./PCA.md)
```
plink \
 --bcf ${bcf} \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.ALL.prune.out \
 --recode A \
 --out data/dapc/${name}.ALL

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA" data/pca/${name}.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.NOLA.prune.out \
 --recode A \
 --out data/dapc/${name}.NOLA

plink \
 --bcf ${bcf} \
 --remove <(grep -E "LA|EX" data/pca/${name}.ALL.nosex) \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --exclude data/pca/${name}.NOLAEX.prune.out \
 --recode A \
 --out data/dapc/${name}.NOLAEX
```
run dapc in R
```{r}

```
