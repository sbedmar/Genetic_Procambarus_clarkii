# fastree R

fastreeR version: 2.1.0 (backend: 0.9.1)

Citation:
Anestis Gkanogiannis (2016).
A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes.
BMC Bioinformatics 17, 311 (2016)
https://doi.org/10.1186/s12859-016-1186-3
https://github.com/gkanogiannis/fastreeR

## install

```
conda create -n fastreer python=3.8
conda activate fastreer
pip install fastreer
```

## run

```
fastreeR VCF2DIST \
    -i data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.vcf \
    -o data/fastreer/output_file data/fastreer/pclarkii.qc_ac_bial.lowmiss_maf.dist
```

## analyze results in R

```{r}



```