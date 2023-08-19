#!/bin/bash

# Loop through each FASTQ file in the current directory
for fastq_file in *.fastq.gz; do
    trim_galore "$fastq_file"

    if [ $? -eq 0 ]; then
        echo "Trimmed $fastq_file"
    else
        echo "Error trimming $fastq_file"
    fi
done

# Loop through each FASTQ file in the current directory
for trimmed_fastq_file in *trimmed.fq.gz; do
    fastqc "$fastq_file"

    if [ $? -eq 0 ]; then
        echo "Analysed $trimmed_fastq_file"
    else
        echo "Error analysing $trimmed_fastq_file"
    fi
done