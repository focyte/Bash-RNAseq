library(DESeq2)
library(ggplot2)
library(readr)
library(ggplot2)
library(ggrepel)

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

# Volcano Plot

# Annotate the results data from the DESeq2 results 
df1 = read.csv("all_genes.csv")
df1$diffexpressed <- "NO"
df1$diffexpressed[df1$log2FoldChange > 0.6 & df1$pvalue < 0.05] <- "UP"
df1$diffexpressed[df1$log2FoldChange < -0.6 & df1$pvalue < 0.05] <- "DOWN"

# Add top gene annotations
df1$delabel <- NA

# Order by p-value
df1 <- df1[order(df1$pvalue), ]
df1$delabel[df1$diffexpressed == "UP"][1:5] <- df1$X[df1$diffexpressed == "UP"][1:5]
df1$delabel[df1$diffexpressed == "DOWN"][1:5] <- df1$X[df1$diffexpressed == "DOWN"][1:5]

# Volcano plot with enhancements
p <- ggplot(data=df1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
  geom_point(alpha=0.6, size=2) +  # Adjust alpha and size
  geom_text_repel(aes(label=delabel), size=3, max.overlaps=10, segment.color="grey50") +  # Repel annotations
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype="dashed") +  # Threshold lines
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
  scale_colour_manual(values=c("DOWN"="blue", "UP"="red", "NO"="grey")) +
  theme_minimal(base_size = 14) +  # Adjust font size
  labs(
    title="Volcano Plot",
    subtitle="Differential Expression Analysis",
    x="Log2 Fold Change",
    y="-Log10 p-value",
    color="Expression"
  ) +
  theme(
    plot.title=element_text(face="bold", size=16),
    legend.position="top"  # Move legend
  )

print(p)

