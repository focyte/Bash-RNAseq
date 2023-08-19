#!/bin/bash

# Loop through each BAM file in the current directory
for bam_file in *.bam; do
    # Generate the output sorted BAM filename by adding "_sorted" to the base name
    sorted_bam_file="${bam_file%.bam}_sorted.bam"

    # Run samtools sort to sort the BAM file
    samtools sort "$bam_file" -o "$sorted_bam_file"

    # Check if the sorting was successful
    if [ $? -eq 0 ]; then
        echo "Sorted $bam_file to $sorted_bam_file"

        # Run samtools index to create an index for the sorted BAM file
        samtools index "$sorted_bam_file"

        # Check if indexing was successful
        if [ $? -eq 0 ]; then
            echo "Indexed $sorted_bam_file"
        else
            echo "Error indexing $sorted_bam_file"
        fi
    else
        echo "Error sorting $bam_file"
    fi
done
