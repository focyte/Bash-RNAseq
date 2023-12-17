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

for file in "$INPUT_DIR"/*.fastq.gz; do
    echo "Processing file: $file"
    hisat2 -x "${INDEX_PATH}/GCA_000001405.15_GRCh38_full_analysis_set" --known-splicesite-infile "${SPLICE_SITES}/human_splice_sites.txt" -p 8 -U "$file" -S "${OUTPUT_DIR}/${file##*/}_trimmed.sam"
done

