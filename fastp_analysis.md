# fastp


```
# generate file with list of IDs
ls Reads_unpacked/Reads/*.fastq.gz | rev | cut -d'_' -f3- | rev | sort -u | rev | cut -d'/' -f1 | rev > Reads_unpacked/Reads/sample_ids.txt

# loop through samples
for n in {1..381}; do
    sample=$(sed -n "${n}p" Reads_unpacked/Reads/sample_ids.txt)
    fastq_dir=Reads_unpacked/Reads/
    r1_fastq=${sample}_R1_001.fastq.gz
    r2_fastq=${sample}_R2_001.fastq.gz

    src/fastp/run_fastp.sh \
        ${fastq_dir} \
        ${r1_fastq} \
        ${r2_fastq}
done    
```



