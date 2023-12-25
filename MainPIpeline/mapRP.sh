#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <index_path> <splice_sites_file>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"
INDEX_PATH="$3"
SPLICE_SITES="$4"

# Generate paired-end file names (file2) by replacing "_R1" with "_R2" in the filenames of the first set of 
# input files (file1) for correct pairing in the subsequent loop iterations
for file1 in "$INPUT_DIR"/*_R1.fastq.gz; do
    file2="${file1/_R1/_R2}"

    echo "Processing raw paired-end files: $file1 and $file2"

    # Extract the sample name from the base filename of file1
    sample_name=$(basename "$file1" | cut -d '_' -f 1)

    # Run hisat2 with the specified options.
    hisat2 -x "${INDEX_PATH}/GCA_000001405.15_GRCh38_full_analysis_set" \
           --known-splicesite-infile "${SPLICE_SITES}/human_splice_sites.txt" \
           -p 8 -1 "$file1" -2 "$file2" -S "${OUTPUT_DIR}/${sample_name}_trimmed.sam"
done
