#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <index_path> <splice_sites_file>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"
INDEX_PATH="$3"
SPLICE_SITES="$4"

# Run individual scripts
bash fastqc.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash trim_fastq.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash fastqcTrimmed.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash map.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
bash samToBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash indexBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
bash featureCount.sh "$INPUT_DIR" "$OUTPUT_DIR"