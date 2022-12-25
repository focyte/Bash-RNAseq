#!/bin/bash

# use file yeast_genome.fa to creat an index file for hisat2
hisat2-build yeast_genome.fa yeast_index

# create a list of splice sites (yeast_splice_sites.txt) to help with mapping using the yeast_gene_models.gtf input file
hisat2_extract_splice_sites.py yeast_gene_models.gtf > yeast_splice_sites.txt

# check the file has been created 
head yeast_splice_sites.txt

# run hisat
hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 2 -U SRR453566_yeast_rnaseq_trimmed.fq.gz -S SRR453566_yeast_rnaseq_trimmed.sam 

# convert the very large .sam file into a smaller .bam file
samtools view -bS SRR453566_yeast_rnaseq_trimmed.sam > SRR453566_yeast_rnaseq_trimmed.bam

# sort the bam file and make an index
samtools sort SRR453566_yeast_rnaseq_trimmed.bam -o SRR453566_yeast_rnaseq_trimmed_sorted.bam 
samtools index SRR453566_yeast_rnaseq_trimmed_sorted.bam 
