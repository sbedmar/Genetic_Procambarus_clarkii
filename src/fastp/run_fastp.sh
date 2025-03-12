#!/bin/bash

source activate clarkii

# Leer el directorio de FASTQ y los archivos R1 y R2
fastq_dir=${1}
r1_fastq=${2}
r2_fastq=${3}

# Convertir la extensión de los archivos de salida a .fastp.fastq.gz
if [[ ${r1_fastq} == *.fastq.gz ]]; then
    r1_fastp=${r1_fastq/.fastq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fastq.gz/.fastp.fastq.gz}
elif [[ ${r1_fastq} == *.fq.gz ]]; then
    r1_fastp=${r1_fastq/.fq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fq.gz/.fastp.fastq.gz}
else
    echo "Error: Los archivos FASTQ deben tener extensión .fastq.gz o .fq.gz"
    exit 1
fi

# Crear el directorio para los archivos procesados si no existe
mkdir -p ${fastq_dir}/fastp
chmod g+w ${fastq_dir}/fastp

# Ejecutar fastp con los parámetros optimizados para tus datos:
fastp \
    -i ${fastq_dir}/${r1_fastq} -I ${fastq_dir}/${r2_fastq} \
    -o ${fastq_dir}/fastp/${r1_fastp} -O ${fastq_dir}/fastp/${r2_fastp} \
    -h ${fastq_dir}/fastp/${r1_fastp/.fastq.gz/_fastp.html} -j ${fastq_dir}/fastp/${r1_fastp/.fastq.gz/_fastp.json} \
    --unpaired1 ${fastq_dir}/fastp/${r1_fastp/.fastq.gz/_unpaired.fastq.gz} \
    --unpaired2 ${fastq_dir}/fastp/${r2_fastp/.fastq.gz/_unpaired.fastq.gz} \
    --failed_out ${fastq_dir}/fastp/${r1_fastp/.fastq.gz/_failed.fastq.gz} \
    --dont_overwrite \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --trim_poly_x \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 50 \
    --thread 1
