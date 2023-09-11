#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Loop through each SAM file in the input directory
for sam_file in "$INPUT_DIR"/*.fq.gz_trimmed.sam; do
    # Generate the output BAM filename by replacing the .sam extension with .bam
    sam_filename=$(basename "$sam_file")
    bam_filename="${sam_filename%.sam}.bam"

    # Run samtools view to convert SAM to BAM
    samtools view -bS "$sam_file" > "$OUTPUT_DIR/$bam_filename"

    # Optionally, you can check if the command was successful and provide some feedback
    if [ $? -eq 0 ]; then
        echo "Converted $sam_file to $OUTPUT_DIR/$bam_filename"
    else
        echo "Error converting $sam_file to $OUTPUT_DIR/$bam_filename"
    fi
done