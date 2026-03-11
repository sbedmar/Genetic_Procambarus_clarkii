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