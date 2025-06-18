## Variant calling from aligned reads

### Haplotype caller

To generate a VCF file from the BAMs, we performed variant calling using [GATK v4.2.6.1](https://gatk.broadinstitute.org/hc/en-us).

The first step is to run [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) using the -ERC GVCF option, to generate genome VCFs of each sample (using the [run_gatk_sample.sh](./src/calling/run_gatk_sample.sh) script). 

```
samples=($(cat Re-sequencing_Reads_tar/Analysis_4229123/sample_ids.txt))

for sample in ${samples[@]:342:38}; do
    echo "calling - ${sample} -"
    bash src/calling/run_gatk_sample.sh ${sample}
done
```

### Create GenomicsDB

To create a GenomicsDB before joint genotyping (see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)) we use the [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport) option of GATK.

I generate a sample map file mapping sample names to the corresponding gvcf and the run GenomicsDBImport:

```
samples=($(cat Re-sequencing_Reads_tar/Analysis_4229123/sample_ids.txt))

for sample in ${samples[@]}; do
    name=$(echo ${sample} | cut -d'_' -f1)
    echo -e "${name}\tdata/gvcfs/${sample}.g.vcf" >> data/gvcfs.sample_map
done

gatk GenomicsDBImport \
   --genomicsdb-workspace-path data/genomicsdb \
   --sample-name-map data/gvcfs.sample_map \
   --intervals data/reference_genome/GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.all_scaffolds.bed \
   --merge-input-intervals \
   --tmp-dir tmp \
   --reader-threads 10
```

### Genotype GVCFs / GenomicsDB

Final step is to call variants using the tool [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs):

```
gatk GenotypeGVCFs \
    -R data/reference_genome/GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna \
    -V gendb://data/genomicsdb \
    -O data/vcfs/pclarkii_allsamples.vcf
```
