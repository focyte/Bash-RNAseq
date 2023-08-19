#!/bin/bash

# Loop through each fq.gz file in the current directory
for fastq_file in *.fq.gz; do
    fastqc "$fastq_file"

    if [ $? -eq 0 ]; then
        echo "Analysed $fastq_file"
    else
        echo "Error analysing $fastq_file"
    fi
done