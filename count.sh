#!/bin/bash

sudo apt install subread

featureCounts -T 4 -s 2 \
  -a yeast_gene_models.gtf \
  -o featurecounts.txt \
  SRR453566_yeast_rnaseq_trimmed_sorted.bam
