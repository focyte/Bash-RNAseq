#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Loop through each BAM file in the input directory
for bam_file in "$INPUT_DIR"/*.bam; do
    # Extract the filename from the full path
    bam_filename=$(basename "$bam_file")

    # Generate the output sorted BAM filename by adding "_sorted" to the base filename
    sorted_bam_file="${bam_filename%.bam}_sorted.bam"

    # Run samtools 'sort' to sort the BAM file
    samtools sort "$bam_file" -o "$OUTPUT_DIR/$sorted_bam_file"

    # Check if the sorting was successful
    if [ $? -eq 0 ]; then
        echo "Sorted $bam_file to $OUTPUT_DIR/$sorted_bam_file"

        # Run samtools 'index' to create an index for the sorted BAM file
        samtools index "$OUTPUT_DIR/$sorted_bam_file"

        # Check if indexing was successful
        if [ $? -eq 0 ]; then
            echo "Indexed $OUTPUT_DIR/$sorted_bam_file"
        else
            echo "Error indexing $OUTPUT_DIR/$sorted_bam_file"
        fi
    else
        echo "Error sorting $bam_file"
    fi
done
