#!/bin/bash

# This script runs GATK HaplotypeCaller on a single sample.
# It requires the following argument:
# 1. Sample name (e.g., "sample1")

sample=${1}

REF=data/reference_genome/GCF_040958095.1_FALCON_Pclarkii_2.0_genomic.fna
IBAM=data/bams/${sample}.sorted.rg.bam
OGVCF=data/gvcfs/${sample}.g.vcf


gatk HaplotypeCaller \
    -R ${REF} \
    -I ${IBAM} \
    -O ${OGVCF} \
    -ERC GVCF \
    --native-pair-hmm-threads 1
