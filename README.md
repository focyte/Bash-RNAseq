# RNA Sequencing Analysis Pipelines

This GitHub repository contains two pipelines for RNA Sequencing analysis: one for initial anlysis of RNA sequencing read data (**Quality Control**) and the other for alignment and mapping of reads to reference genome and counting of features (genes) (**Main Pipeline**). Each pipeline consists of a series of Bash scripts that automate key steps in RNA sequencing data analysis, along with additional Python and R scripts for downstream analysis.

---

## Table of Contents
1. [Requirements](#requirements)
2. [Pipeline](#pipeline)
    - [Quality Control Usage](#quality-control-usage)
    - [Main Pipeline Usage](#main-pipeline-usage)
    - [Downstream Analysis](#downstream-analysis)
    - [Results](#results)
3. [How to Run](#how-to-run)
4. [Results](#results)

---


## Requirements

### Software
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Hisat2](https://daehwankimlab.github.io/hisat2/)
- [Samtools](http://www.htslib.org/)
- [FeatureCounts](http://subread.sourceforge.net/)

### Files
- Sequencing read data in the fastq.gz format
- Index files for the reference genome of interest, in this case Human Genome hg38
- Ideally perform your own indexing using software such as [STAR](https://github.com/alexdobin/STAR) aligner
- A .gtf file of annotated features for your indexed genome 
- Splice Site file for your indexed genome to improve alignment accuracy across exon-exon boundaries

## Pipeline

## Quality Control Usage

```bash
./runQC.sh <input_dir> <output_dir>
```
### Individual Steps

1. **FastQC Analysis**

   - Script: `fastqc.sh`
   - Usage: 
    ```bash 
    fastqc.sh "$INPUT_DIR" "$OUTPUT_DIR"
    ```

2. **Trimming with Trimgalore**

   - Script: `trim_fastq.sh`
   - Usage: 
   ```bash 
   trim_fastq.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

3. **FastQC Analysis on Trimmed Data**

   - Script: `fastqcTrimmed.sh`
   - Usage: 
   ```bash 
   fastqcTrimmed.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

## Main Pipeline Usage

```bash
./runPipeline.sh <input_dir> <output_dir> <index_path> <splice_sites_file> <gtf_file> <read_type> <data_type>
```

1. **Mapping to Human Genome using Hisat2**

   - Script: `mapPP.sh, mapPU.sh, mapRP.sh, mapRE.sh`
   - When specifying the read_type and data_type in runPipeline.sh, IF statements determine which mapping script to use
   - Read_types = Unpaired OR Paired
   - Data_types = Raw OR Processed (Raw will use files processed by runQC.sh in the Quality Control step)
   - Usage: 
   ```bash 
   map.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
   ```

2. **Conversion of SAM to BAM**

   - Script: `samToBam.sh`
   - Usage: 
   ```bash 
   samToBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

3. **Indexing BAM Files**

   - Script: `indexBam.sh`
   - Usage: 
   ```bash 
   indexBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

4. **Counting Reads for Each Gene Feature using FeatureCounts**

   - Script: `featureCount.sh`
   - Usage: 
   ```bash 
   featureCount.sh "$INPUT_DIR" "$OUTPUT_DIR" "$GTF_FILE"
   ```

## Downstream Analysis

### Python Script for Merging FeatureCounts Results

   - Script: `merge_featureCounts.py`
   - Usage: 
   ```python 
    merge_featureCounts.py file_paths output_path
   ```

### R Script for DSeq2 Analysis

   - Script: `DSeq2_analysis.R`
   - Usage: `Execute in an R environment`

## Results

![Figure](https://github.com/focyte/Bash-RNAseq/blob/main/Figure.png) 
