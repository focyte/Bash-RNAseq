#!/bin/bash

# install the required tools
sudo apt-get -y install hisat2
sudo apt-get -y install samtools
sudo apt intall trim-galore
pip install cutadapt
sudo apt install fastqc
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
sudo apt-get install rna-star

# navigate to the folder
cd ./Documents/Human_RNASeq

# run fastqc on all files
find . -maxdepth 1 -name \*fastq.gz -exec fastqc {} \;

# list files to check that it has made the correct files
ls -ltr

# This did not work, therefore ran without --fastqc
find . -maxdepth 1 -name \*fastq.gz -exec trim_galore {} \;

#then run the report again
find . -maxdepth 1 -name \*trimmed.fq.gz -exec fastqc {} \;
