#!/bin/bash

# Loop through each SAM file in the current directory
for sam_file in *.sam; do
    # Generate the output BAM filename by replacing the .sam extension with .bam
    bam_file="${sam_file%.sam}.bam"

    # Run samtools view to convert SAM to BAM
    samtools view -bS "$sam_file" > "$bam_file"

    # Optionally, you can check if the command was successful and provide some feedback
    if [ $? -eq 0 ]; then
        echo "Converted $sam_file to $bam_file"
    else
        echo "Error converting $sam_file to $bam_file"
    fi
done