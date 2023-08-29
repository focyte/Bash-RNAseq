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

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

for forward_read in "$INPUT_DIR"/*_1.fastq.gz; do
    reverse_read="${forward_read/_1/_2}"
    output_name="$(basename "$forward_read" _1.fastq.gz)"

    hisat2 -x "${INDEX_PATH}/GCA_000001405.15_GRCh38_full_analysis_set" --known-splicesite-infile "${SPLICE_SITES}/human_splice_sites.txt" -p 8 -1 "$forward_read" -2 "$reverse_read" -S "${OUTPUT_DIR}/${output_name}.sam"
done
