import pandas as pd

file1 = pd.read_csv("C:\Users\focyt\OneDrive\Documents\Project20230906\HCC1954_siControl_1.fastq.gz_trimmed_featurecounts.txt", sep = '\t', skiprows=[1])
file2 = pd.read_csv("C:\Users\focyt\OneDrive\Documents\Project20230906\HCC1954_siControl_3.fastq.gz_trimmed_featurecounts.txt", sep = '\t', skiprows=[1])
file3 = pd.read_csv("C:\Users\focyt\OneDrive\Documents\Project20230906\HCC1954_siRNF40_1.fastq.gz_trimmed_featurecounts.txt", sep = '\t', skiprows=[1])
file4 = pd.read_csv("C:\Users\focyt\OneDrive\Documents\Project20230906\HCC1954_siRNF40_2.fastq.gz_trimmed_featurecounts.txt", sep = '\t', skiprows=[1])
file5 = pd.read_csv("C:\Users\focyt\OneDrive\Documents\Project20230906\HCC1954_siRNF40_3.fastq.gz_trimmed_featurecounts.txt", sep = '\t', skiprows=[1])

print(file1)


