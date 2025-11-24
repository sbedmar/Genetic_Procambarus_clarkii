##  Filter variants from the calling 

```
# change genotypes with <5 reads to missing data 
# (just realized that GQ is 3xDP so filter is coded weird and shoudld be changed)
# probably just to GQ < X
bcftools +setGT data/vcfs/pclarkii_allsamples.vcf -- \
    -t q -n . -i 'FORMAT/DP<5 | FORMAT/GQ<20' \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.vcf

# recalculate allele frequency and allele number
bcftools +fill-tags data/vcfs/pclarkii_allsamples.min5dp20gq.vcf -- \
    -t AF,AN,AC \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.vcf

# remove allele count=0 and low quality variants
bcftools view -e \
    'INFO/MQRankSum < -12.5 | INFO/ReadPosRankSum < -8.0 | INFO/QD < 6.0 | INFO/FS > 60.0 | INFO/MQ < 40.0 | INFO/AC < 1 | INFO/ExcessHet > 10' \
    data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.vcf \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.vcf

# remove INDELs and non-biallelic SNPs
bcftools view -v snps -m2 -M2 \
    data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.vcf \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.vcf

# filter missing data by population
for pop in $(grep -m1 "^#CHROM" data/vcfs/pclarkii_allsamples.min5dp20gq.vcf | cut -f10- | tr '\t' '\n' | cut -c 1-2 | sort -u); do
    echo $pop
    grep -m1 "^#CHROM" data/vcfs/pclarkii_allsamples.min5dp20gq.vcf | cut -f10- | tr '\t' '\n' | grep "^${pop}" > data/pop_${pop}.txt
    bcftools view -S data/pop_${pop}.txt data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.vcf.gz | \
        bcftools view -i 'COUNT(GT!="mis")>=3' -Oz -o data/vcfs/${pop}_filtered.vcf.gz
    tabix data/vcfs/${pop}_filtered.vcf.gz
done
N=$(ls data/vcfs/*_filtered.vcf.gz | wc -l)
bcftools isec -n=$N -w1 -Oz \
    -o data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min3perpop.vcf.gz \
    data/vcfs/*_filtered.vcf.gz
gunzip data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min3perpop.vcf.gz
grep -v "#" data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min3perpop.vcf |
    cut -f 1,2 | awk '{print $1, $2-1, $2}' | tr ' ' '\t' > data/vcfs/min3perpop_snps.bed
bedtools intersect -header \
    -a data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.vcf.gz \
    -b data/vcfs/min3perpop_snps.bed \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min3perpop.vcf

# make raw
plink --vcf data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min3perpop.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --recode A --out data/vcfs/pclarkii_allsamples.min3perpop

##### MIN 4 BELOW

# filter missing data by population
for pop in $(grep -m1 "^#CHROM" data/vcfs/pclarkii_allsamples.min5dp20gq.vcf | cut -f10- | tr '\t' '\n' | cut -c 1-2 | sort -u); do
    echo $pop
    grep -m1 "^#CHROM" data/vcfs/pclarkii_allsamples.min5dp20gq.vcf | cut -f10- | tr '\t' '\n' | grep "^${pop}" > data/pop_${pop}.txt
    bcftools view -S data/pop_${pop}.txt data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.vcf.gz | \
        bcftools view -i 'COUNT(GT!="mis")>=4' -Oz -o data/vcfs/${pop}_filtered4.vcf.gz
    tabix data/vcfs/${pop}_filtered4.vcf.gz
done
N=$(ls data/vcfs/*_filtered.vcf.gz | wc -l)
bcftools isec -n=$N -w1 -Oz \
    -o data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min4perpop.vcf.gz \
    data/vcfs/*_filtered4.vcf.gz
gunzip data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min4perpop.vcf.gz
grep -v "#" data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min4perpop.vcf |
    cut -f 1,2 | awk '{print $1, $2-1, $2}' | tr ' ' '\t' > data/vcfs/min4perpop_snps.bed
bedtools intersect -header \
    -a data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.vcf.gz \
    -b data/vcfs/min3perpop_snps.bed \
    > data/vcfs/pclarkii_allsamples.min5dp20gq.afanac.qual.bial_snps.min4perpop.vcf
```

## BCFTOOLS
```
# filter to include high quality and biallelic variants
bcftools +setGT -Ou data/vcfs/pclarkii_allsamples.vcf -- \
    -t q -n . -i 'FORMAT/GQ < 15' |
bcftools +fill-tags -Ou -- \
    -t AF,AN,AC |
bcftools view -Ou -e \
    'INFO/MQRankSum < -12.5 | INFO/ReadPosRankSum < -8.0 | INFO/QD < 6.0 | INFO/FS > 60.0 | INFO/MQ < 40.0 | INFO/AC < 1 | INFO/ExcessHet > 10' |
bcftools annotate -Ou \
    -x ^INFO/AC,^INFO/AN,^INFO/AF |
bcftools view -Ou \
   -v snps -m2 -M2 \
> data/vcfs/pclarkii.qc_ac_bial.bcf

# calculate individual missingness
vcftools --bcf data/vcfs/pclarkii.qc_ac_bial.bcf --missing-indv
mv out.imiss data/

# remove high missing individuals
# remove high missing and low frequency snps
bcftools view -Ou \
    -S <(awk '$5 < 0.95' data/out.imiss | cut -f1) data/vcfs/pclarkii.qc_ac_bial.bcf |
bcftools view -Ou \
    -i 'F_MISSING < 0.2 & MAF > 0.02' \
> data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf

# make bed
plink --make-bed \
 --bcf data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --out data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf

# Total genotyping rate is 0.853602.
# 40335 variants and 328 people pass filters and QC.

plink --recode A \
 --bcf data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --out data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf

plink --bcf data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --pca 'header' \
    --out pclarkii.qc_ac_bial.lowmiss_maf.plink.pca

```

## stacks: populations module

`populations -V vcf -O dir [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)`

```
# input vcf:
invcf=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.vcf
# popmap file:
paste \
  <(bcftools query -l data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf) \
  <(bcftools query -l data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf | cut -c 1,2) \
> data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.popmap
popmap=data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.popmap
# outputdirectory
odir=data/stacks_populations

populations \
  -V $invcf \
  -O $odir \
  -M $popmap \
  --fstats
```