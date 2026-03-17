# admixture

[admixture](https://rajanil.github.io/fastStructure/) software from [Alexander et al. 2009](10.1101/gr.094052.109)

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

plink --make-bed \
 --bcf ${bcf} \
 --double-id --allow-extra-chr 0 \
 --out data/vcfs/${name}.admixture

```

run admixture 10 times at different ks:
```
for i in {1..9}; do
    for k in {3..10}; do
        admixture -j10 --cv data/vcfs/${name}.admixture.bed ${k} \
            > data/admixture/${name}.admixture.log.${k}.${i}.out
        mv ${name}.admixture.${k}.P data/admixture/${name}.admixture.${k}.${i}.P
        mv ${name}.admixture.${k}.Q data/admixture/${name}.admixture.${k}.${i}.Q
    done
done

paste <(grep -h CV data/admixture/${name}.admixture.log.*.out | cut -d'=' -f2 | cut -d')' -f1) \
    <(grep -h CV data/admixture/${name}.admixture.log.*.out | cut -d' ' -f4) \
    > data/admixture/cv_table.txt
```
analyze in R:
```{r}

```
