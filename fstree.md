# Fst with plink + NJ tree

First choose bcf:
```
# all samples filtered bcf
bcf=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf
name=$(basename -s .bcf ${bcf})
```

Build cluster file with populations for [plink --within](https://www.cog-genomics.org/plink/1.9/input#within) 
```
paste <(cat data/pop_*) <(cat data/pop_*) <(cat data/pop_* | cut -c 1,2) | \
    grep -f <(cut -d' ' -f1 data/vcfs/${name}.fam) \
    > data/fid_iid_pop.txt
```

Get [fst from plink](https://www.cog-genomics.org/plink/1.9/basic_stats)
```
plink \
 --bcf ${bcf} \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --within data/fid_iid_pop.txt --fst \
 --out data/vcfs/${name}
```

# Fst with vcftools then tree

```
pops=($(ls data/pop_* | cut -d'_' -f2 | cut -d'.' -f1))

for ((i=0; i<${#pops[@]}; i++)); do
  for ((j=i+1; j<${#pops[@]}; j++)); do
    echo "${pops[i]}_${pops[j]}"
    vcftools --vcf data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.vcf \
        --weir-fst-pop data/pop_${pops[i]}.txt \
        --weir-fst-pop data/pop_${pops[j]}.txt \
        --out data/fst/${pops[i]}_${pops[j]}
  done
done

grep "weighted" data/fst/*.log | \
    cut -d' ' -f1,7 | cut -d'/' -f3 | sed 's/.log:Weir//' \
    > data/fst/all_pairs.weighted.fst

grep "mean" data/fst/*.log | \
    cut -d' ' -f1,7 | cut -d'/' -f3 | sed 's/.log:Weir//' \
    > data/fst/all_pairs.mean.fst
```
analyze results in R:
```{r}
```