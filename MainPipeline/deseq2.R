library(DESeq2)
library(ggplot2)
library(readr)

countData <- read_csv("mergedCounts.csv")
countData <- as.data.frame(countData)
rownames(countData) <- countData$Geneid
countData <- countData[, -1]


# Create a data frame to store the sample group information
# Replace 'Group1&2' your sample groups
colData <- data.frame(
  Group = factor(c("Group1", "Group1", "Group1", "Group2", "Group2", "Group2")),
  row.names = colnames(countData)
)

# Create a DESeqDataSet object (dds)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Group)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Differential expression analysis 
# Comparing Group1 to Group2:
results <- results(dds, contrast = c("Group", "Group2", "Group1"))
write.csv(results, file = "all_genes.csv")

# Extract significant differentially expressed genes based on a set threshold
significant_genes <- subset(results, padj < 0.05)

# View the top differentially expressed genes
top_genes <- head(significant_genes, n = 10)

# Print the results
print(top_genes)
write.csv(significant_genes, file = "significant_genes.csv")

# Generate the variance-stabilized transformed data
vsd <- vst(dds, blind=FALSE)

# Perform PCA on the transformed data
pcaData <- plotPCA(vsd, intgroup=c("Group"), returnData=TRUE)

# Calculate the percentage of variance for the axes
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the PCA plot
ggplot(pcaData, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) +
  scale_color_manual(
    values = c("Group1" = "blue", "Group2" = "red"), # Custom colors
    labels = c("Group1" = "RSM3", "Group2" = "Control") # Custom labels
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA plot of DESeq2 Transformed Data")
