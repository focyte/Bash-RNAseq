# install the required tools
sudo apt-get -y install hisat2
sudo apt-get -y install samtools
sudo apt intall trim-galore
pip install cutadapt
sudo apt install fastqc
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz

# navigate to the folder
cd ./Documents/RNASeq_Course

# run fastqc on files
fastqc SRR453566_yeast_rnaseq.fq.gz

# list files to check that it has made the correct files
ls -ltr

# open results of fastqc analysis in chrome
google-chrome SRR453566_yeast_rnaseq_fastqc.html

# trim the reads, adding --fastqc performs another analysis of the reads
trim_galore –-fastqc SRR453566_yeast_rnaseq.fq.gz

# This did not work, therefore ran without --fastqc
trim_galore SRR453566_yeast_rnaseq.fq.gz

#then run the report again
fastqc SRR453566_yeast_rnaseq_trimmed.fq.gz

# open results in google chrome
google-chrome SRR453566_yeast_rnaseq_trimmed_fastqc.html
