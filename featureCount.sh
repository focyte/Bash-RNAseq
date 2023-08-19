#!/bin/bash

# Provide the GTF file name
gtf_file="./GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"

# Loop through each sorted BAM file in the current directory
for sorted_bam_file in *_sorted.bam; do
    # Generate the output featureCounts file name by replacing "_sorted.bam" with "_featurecounts.txt"
    featurecounts_file="${sorted_bam_file%_sorted.bam}_featurecounts.txt"

    # Run featureCounts with the specified options
    featureCounts -T 4 -s 2 -a "$gtf_file" -o "$featurecounts_file" "$sorted_bam_file"

    # Check if featureCounts was successful
    if [ $? -eq 0 ]; then
        echo "Feature counts calculated for $sorted_bam_file and saved in $featurecounts_file"
    else
        echo "Error calculating feature counts for $sorted_bam_file"
    fi
done
