# RNA Sequencing Analysis Pipelines

This GitHub repository contains two distinct pipelines for RNA Sequencing analysis: one for raw RNA sequencing read data (**RawReads**) and the other for pre-processed data (**PreProcessed**). Each pipeline consists of a series of Bash scripts that automate key steps in RNA sequencing data analysis, along with additional Python and R scripts for downstream analysis.

## RawReads Pipeline

### Requirements
#### Software
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Hisat2](https://daehwankimlab.github.io/hisat2/)
- [Samtools](http://www.htslib.org/)
- [FeatureCounts](http://subread.sourceforge.net/)

#### Files
- Sequencing read data in the fastq.gz format
- Index files for the reference genome of interest, in this case Human Genome hg38.
- Ideally perform your own indexing using STAR aligner or similar tool.
- A .gtf file of annotated features of your indexed genome 
- human_splice_sites file for your indexed genome

### Usage
```bash
./runPipeline.sh <input_dir> <output_dir> <index_path> <splice_sites_file> <gtf_file>
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

4. **Mapping to Human Genome using Hisat2**

   - Script: `map.sh`
   - Usage: 
   ```bash 
   map.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
   ```

5. **Conversion of SAM to BAM**

   - Script: `samToBam.sh`
   - Usage: 
   ```bash 
   samToBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

6. **Indexing BAM Files**

   - Script: `indexBam.sh`
   - Usage: 
   ```bash 
   indexBam.sh "$INPUT_DIR" "$OUTPUT_DIR"
   ```

7. **Counting Reads for Each Gene Feature using FeatureCounts**

   - Script: `featureCount.sh`
   - Usage: 
   ```bash 
   featureCount.sh "$INPUT_DIR" "$OUTPUT_DIR" "$GTF_FILE"
   ```

## PrProcessed Pipeline

### Usage
```bash
./runPipelinePreP.sh <input_dir> <output_dir> <index_path> <splice_sites_file> <gtf_file>
```
### Individual Steps

1. **Mapping to Human Genome using Hisat2**

   - Script: `map2.sh`
   - Usage: 
   ```bash 
   map.sh "$INPUT_DIR" "$OUTPUT_DIR" "$INDEX_PATH" "$SPLICE_SITES"
   ```

2. **Conversion of SAM to BAM**

   - Script: `samToBam.sh2`
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
    merge_featureCounts.py
    ```

### R Script for DSeq2 Analysis

   - Script: `DSeq2_analysis.R`
   - Usage: `Execute in an R environment`