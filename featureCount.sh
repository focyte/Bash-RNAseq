#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <gtf_file>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
GTF_FILE="$3"

# Loop through each sorted BAM file in the input directory
for sorted_bam_file in "$INPUT_DIR"/*_sorted.bam; do
    # Generate the output featureCounts file name by replacing "_sorted.bam" with "_featurecounts.txt"
    featurecounts_file="${sorted_bam_file%_sorted.bam}_featurecounts.txt"

    # Run featureCounts with the specified options
    featureCounts -T 4 -s 2 -a "$GTF_FILE"/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf -o "$featurecounts_file" "$sorted_bam_file"


    # Check if featureCounts was successful
    if [ $? -eq 0 ]; then
        echo "Feature counts calculated for $sorted_bam_file and saved in $OUTPUT_DIR/$featurecounts_file"
    else
        echo "Error calculating feature counts for $sorted_bam_file"
    fi
done