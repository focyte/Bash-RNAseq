#!/bin/bash
for file in ./*.fastq; do
    hisat2 [options] -x ./Index/ -p 8 -U "$file" -S "${file%.fastq}.sam" --known-splicesite-infile ./Index/human_splice_sites.txt
done
