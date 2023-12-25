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

for file1 in "$INPUT_DIR"/*_R1.fastq.gz; do
    file2="${file1/_R1/_R2}"

    echo "Processing paired-end files: $file1 and $file2"

    sample_name=$(basename "$file1" | cut -d '_' -f 1)

    hisat2 -x "${INDEX_PATH}/GCA_000001405.15_GRCh38_full_analysis_set" \
           --known-splicesite-infile "${SPLICE_SITES}/human_splice_sites.txt" \
           -p 8 -1 "$file1" -2 "$file2" -S "${OUTPUT_DIR}/${sample_name}_trimmed.sam"
done

