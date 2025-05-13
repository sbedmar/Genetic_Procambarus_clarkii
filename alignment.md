## Alignment to reference genomes

I align the sequencing reads to the [*Procambarus clarkii* reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040958095.1).

The reference genome is located in the `data/reference_genome` folder.

A script will be run for each sample that will use bwa v0.7.17, samtools v1.9, picard v2.25.5 and gatk v3.7-0 that will do the following:

- Align each of the sample's R1-R2 fastq pairs to the reference genome using BWA-MEM and Samtools view
- samtools sort to sort the reads in the bam file
- picardtools AddOrReplaceReadGroups to add read groups to the reads in the bam file
- samtools merge to merge the bam files from the multiple R1-R2 pairs of the sample if necessary
- picardtools MarkDuplicates to mark duplicate reads in the bam
- realign indels with GATK which comprises:
  - gatk RealignerTargetCreator to identify realign targets
  - gatk IndelRealigner to realign the indels
- samtools index to index the indel realigned bam for downstream analyses

To prepare the reference genome of the alignment pipeline I run the following bwa and samtools commands:
```
cd data/reference_genome
bwa index GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna
samtools faidx GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna
samtools dict GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna -o GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.dict
```

Run alignments:
```
mkdir data/bams

samples=($(cat Re-sequencing_Reads_tar/Analysis_4229123/sample_ids.txt))
for sample in ${samples[@]}; do
    echo "${sample}:"
    bwa mem \
        data/reference_genome/GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna \
        Re-sequencing_Reads_tar/Analysis_4229123/fastp/${sample}_R1_001.fastp.fastq.gz \
        Re-sequencing_Reads_tar/Analysis_4229123/fastp/${sample}_R2_001.fastp.fastq.gz \
        -t 10 |
        samtools view -hbS -@ 10 - -o data/bams/${sample}.bam
done
```

Sort bam and add read groups:
```
samples=($(cat Re-sequencing_Reads_tar/Analysis_4229123/sample_ids.txt))
for sample in ${samples[@]; do
    
    echo " - sorting ${sample}.bam"
    samtools sort \
        -@ 10 \
        data/bams/${sample}.bam \
        -o data/bams/${sample}.sorted.bam
    
    id1=$(echo ${sample} | cut -d'_' -f1)
    id2=$(echo ${sample} | cut -d'_' -f2)

    echo " - adding read groups to ${sample}.sorted.bam"
    java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups \
        I=data/bams/${sample}.sorted.bam \
        O=data/bams/${sample}.sorted.rg.bam \
        RGID=${sample} \
        RGLB=${id2} \
        RGPL=ILLUMINA \
        RGPU=${sample} \
        RGSM=${id1} \
        VALIDATION_STRINGENCY=SILENT

done
```

Index:
```
samples=($(cat Re-sequencing_Reads_tar/Analysis_4229123/sample_ids.txt))
for sample in ${samples[@]}; do

    # Indexing for GATK
    echo " - Indexing ${sample} -"
    samtools index data/bams/${sample}.sorted.rg.bam
    
done
```