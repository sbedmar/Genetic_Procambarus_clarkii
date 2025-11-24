# fastStructure

[fastStructure](https://rajanil.github.io/fastStructure/) software from [Raj et al. 2014](https://doi.org/10.1534/genetics.114.164350)

First choose bcf:
```
# all samples filtered bcf
bcf=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf

name=$(basename -s .bcf ${bcf})
```

make bed from bcf
```
plink --make-bed \
 --bcf ${bcf} \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --out data/vcfs/${name}
```

run fastStructure 10 times at different ks:
```
for i in {1..10}; do
    for k in {2..8}; do
        fastStructure -K ${k} \
            --input=data/vcfs/${name} \
            --output=data/fastStructure/${name}.i${i} \
            --cv=5
    done
    cat data/fastStructure/${name}.i${i}.*.log | grep -h CV | tr -d ',' | cut -d' ' -f4,5 \
        > data/fastStructure/${name}.i${i}.cv_table.txt
done
```

analyze in R:
```{r}
source("src/fastStructure/functions.R", print.eval=TRUE)
dataname <- "pclarkii.qc_ac_bial.lowmiss_maf"

# read in samples from a fam file (from plink --make-bed)
famfile <- paste0("data/vcfs/", dataname,".fam")
samples <- get_samples_from_fam(famfile)

# qvalues -> admix dataframe for plotting
for(i in 1:10){
  for (k in 2:8){
    meanqfile <- paste0("data/fastStructure/", dataname,".i", i, ".", k, ".meanQ")
    admix <- read_qvals(meanqfile, samples)
    # plot!
    admix_plot <- plot_admix(admix)
    ggsave(
      filename = paste0("plots/fastStructure/k", k,"/", dataname,".i", i, ".", k, ".pdf"),
      plot = admix_plot,
      width = 25, height = 6, units = "cm"
    )
  }
}

# what k?
cv_table <- data.frame()
for (k in 2:8){
  for(i in 1:10){
    cv <- read.table(paste0(
      "data/fastStructure/pclarkii.qc_ac_bial.lowmiss_maf.i",
      i, ".cv_table.txt"
      ))
    cv_table <- rbind(cv_table, data.frame(
      "K" = k,
      "CVerr" = cv[k - 1, 2],
      "i" = i
    ))
  }
}
cv_table$i <- factor(cv_table$i)
# plot cv values
for (n in c(1:10)){
  cvplot <- ggplot() +
  geom_point(data = cv_table |> filter(i==n),
    aes(x = K, y = CVerr, fill = i),
    shape = 21, size = 2.5) +
  xlab(paste0("Number of Populations")) +
  ylab(paste0("Cross Validation Error")) +
  scale_fill_discrete() +
  scale_x_continuous(n.breaks = 8) +
  theme_bw()
  ggsave(
    filename = paste0("plots/fastStructure/", dataname,".i", n, ".cvplot.pdf"),
    plot = cvplot,
    width = 6, height = 6, units = "cm"
    )

}
cvplot <- ggplot() +
  geom_point(data = cv_table,
    aes(x = K, y = CVerr, fill = i),
    shape = 21, size = 2.5) +
  xlab(paste0("Number of Populations")) +
  ylab(paste0("Cross Validation Error")) +
  scale_fill_discrete() +
  scale_x_continuous(n.breaks = 8) +
  theme_bw()
```
