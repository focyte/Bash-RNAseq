for file in /home/john/Documents/project2/*.fastq; do
    hisat2 [options] -x /home/john/Documents/project2/Index -p 8 -U "$file" -S "${file%.fastq}.sam" --known-splicesite-infile /path/to/splicing_annotation.gtf
done
