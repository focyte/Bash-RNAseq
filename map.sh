#!/bin/bash
index_path="./Index/GCA_000001405.15_GRCh38_full_analysis_set"
splice_sites="./Index/human_splice_sites.txt"

for file in ./*trimmed.fq.gz; do
    hisat2 -x "$index_path" --known-splicesite-infile "$splice_sites" -p 8 -U "$file" -S "${file%trimmed.fq.gz}.sam"
done