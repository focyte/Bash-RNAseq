#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <index_path> <splice_sites_file>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
INDEX_PATH="$3"
SPLICE_SITES="$4"

for file in "$INPUT_DIR"/*trimmed.fq.gz; do
    hisat2 -x "$INDEX_PATH" --known-splicesite-infile "$SPLICE_SITES" -p 8 -U "$file" -S "${OUTPUT_DIR}/${file##*/trimmed.fq.gz}.sam"
done
