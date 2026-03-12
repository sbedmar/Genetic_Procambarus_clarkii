## admixtools

install from conda:
```
conda activate clarkii
conda install -c bioconda eigensoft
```

### prep input

filter out non "nc_" chromosomes (smaller ones):
```
bcftools view -h data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf.gz | grep -v "contig=<ID=NW_" > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf
bcftools view -H data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf.gz | awk '$1 ~ /^NC_/' >> data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf
```
make plink file
```
plink \
  --vcf data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf \
  --double-id \
  --allow-extra-chr \
  --chr-set 94 \
  --make-bed \
  --out data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr
```
make a chromosome to number map and change bim chromosome names to numbers:
```
# map
paste -d ' ' \
    <(grep -v "#" data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf | cut -f1 | uniq) \
    <(for n in {1..94}; do echo $n; done) \
    > data/admixtools/chrom_map.txt
# change names:
awk 'NR==FNR {map[$1]=$2; next} {$1=map[$1]; print}' data/admixtools/chrom_map.txt \
    data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.bim > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bim
cp data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.bed data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bed
cp data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.fam data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.fam
```

from and Sonia's code to generate input:
```
nano data/admixtools/par.PLINK2EIGEN

# paste this in file:
genotypename:    data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bed
snpname:         data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bim
indivname:       data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.fam
outputformat:    EIGENSTRAT
genotypeoutname: data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.geno
snpoutname:      data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.snp
indoutname:      data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.ind
familynames:     NO
####### 

conda activate clarkii
convertf -p data/admixtools/par.PLINK2EIGEN
```

for ind file:
```
paste -d ' ' \
    <(bcftools query -l data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf) \
    <(yes "U" | head -n 328) \
    <(bcftools query -l data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf | cut -c 1,2) \
    > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.ind
```


[admixr](https://cran.r-project.org/web/packages/admixr/vignettes/vignette-01-tutorial.html)