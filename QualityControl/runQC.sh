#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"

bash fastqc.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash trim_fastq.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash fastqcTrimmed.sh "$INPUT_DIR" "$OUTPUT_DIR"