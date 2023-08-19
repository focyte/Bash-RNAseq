#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Loop through each fastq.gz file in the input directory
for fastq_file in "$INPUT_DIR"/*.fastq.gz; do
    fastqc "$fastq_file" -o "$OUTPUT_DIR"

    if [ $? -eq 0 ]; then
        echo "Analysed $fastq_file"
    else
        echo "Error analysing $fastq_file"
    fi
done