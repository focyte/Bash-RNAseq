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

# Volcano Plot

# Annotate the results data from the DESeq2 results 
df1 = results
df1$diffexpressed <- "NO"
df1$diffexpressed[df1$log2FoldChange > 0.6 & df1$pvalue < 0.05] <- "UP"
df1$diffexpressed[df1$log2FoldChange < -0.6 & df1$pvalue < 0.05] <- "DOWN"

# Plot the data
p <- ggplot(data=df1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

df1$delabel <- NA
df1$delabel[df1$diffexpressed != "NO"] <- df1$X[df1$diffexpressed != "NO"]

ggplot(data=df1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal()
