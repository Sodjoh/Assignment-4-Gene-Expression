# script to perform differential gene expression analysis using DESeq2 and edgeR package with reference to Tong et al., 2020
#I couldn't find my way around NCBI GEO, so I settled for EMBL-EBI
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("apeglm", "Glimma", "Mus.musculus"))
install.packages("apeglm")
# load all libraries
library(DESeq2)
library(tidyverse)
library(edgeR)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(pheatmap)
library(apeglm)
#Now let us prepare the data
#download the raw count data from ebi
count_matrix <- read.delim("raw-counts.tsv")
head(count_matrix)
class (count_matrix)
dim(count_matrix)
head(count_matrix[, 1:4])
colnames(count_matrix)

#download experimental data
sample_info <- read.delim("experiment-design.tsv")
head(sample_info)
class(sample_info)
colnames(sample_info)

#DO A LITTLE BIT OF FILTERING AND EXPLORATION
#Ensure the variables are in the right matrix form, i.e the genes in rows and samples as columns
head(count_matrix)
rownames(count_matrix) <- count_matrix$Gene.ID
head(count_matrix)
#So we can get ready to take the other columns that are not needed off
# Remove non-sample columns (Gene.ID and Gene.Name)
genes <- count_matrix[, c("Gene.ID", "Gene.Name")]
newcount <- count_matrix[, -c(1, 2)]
head(newcount)

head(sample_info)
rownames(sample_info) <- sample_info$Run
head(sample_info)

#only keep columns of interest
sample_info <- sample_info[, c("Sample.Characteristic.genotype.", "Sample.Characteristic.individual."), drop=FALSE]
sample_info
#so rename the columns of interest
colnames(sample_info) <- c("Genotype", "abbreviation")
sample_info

#trying to rename the long genotype
sample_info$Genotype[sample_info$Genotype =='Ptf1aCre/+; LSL-KrasG12D/+; mir-802+/+'] <- 'mir802_normal'
sample_info$Genotype[sample_info$Genotype == 'Ptf1aCre/+; LSL-KrasG12D/+; mir-802fl/fl'] <- 'mir802_floxed'
sample_info


#Turn genotype into a factor
#Ensuring sample names match with row names, turning everything into a factor
sample_info <- data.frame(
  sampleName = c("ERR5873116","ERR5873117","ERR5873118","ERR5873119",
                 "ERR5873120","ERR5873121","ERR5873122","ERR5873123"),
  genotype  = c(
    "mir802_normal",
    "mir802_normal",
    "mir802_normal",
    "mir802_normal",
    "mir802_floxed",
    "mir802_floxed",
    "mir802_floxed",
    "mir802_floxed"
  ),
  condition = factor(c(rep("WT",4), rep("KO",4)), levels = c("WT","KO")),
  stringsAsFactors = FALSE,
  row.names = c("ERR5873116","ERR5873117","ERR5873118","ERR5873119",
                "ERR5873120","ERR5873121","ERR5873122","ERR5873123")
)

# making sure the row names in sample_info matches to column names in count_matrix
all(colnames(newcount) %in% rownames(sample_info))

# are they in the same order?
all(colnames(newcount) == rownames(sample_info))
#Plot a simple graph showing substantial differences in library size that actually needs normalization
library_sizes <- colSums(newcount)
df_lib <- data.frame(Sample = names(library_sizes),LibrarySize = library_sizes)
ggplot(df_lib, aes(x = Sample, y = LibrarySize)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "Library Sizes Across Samples", x = "Sample ID", y = "Total Read Counts")

#Now, let us proceed to why we are here
# 1. Create the DESeqDataSet object
# genotype defines the comparison groups (wildtype vs knockout)
# confirming my variables are the exact way I want it
ncol(newcount) == nrow(sample_info)
ncol(newcount)
nrow(sample_info)
colnames(newcount)
rownames(sample_info)

# with respect to the paper, I will majorly focus on computational speed, accuracy, precision, and reliability
#I am adding a time start of the DESeq2 analysis to check speed according to tong et al 2020
time_start_deseq2 <- Sys.time()
dds <- DESeqDataSetFromMatrix(countData = newcount, colData=sample_info, design=~genotype)
#Keep only genes with at least 10 reads total
dds <- dds[rowSums(counts(dds))>10,]

# The below function performs:
# a) Normalization (Median of Ratios)
# b) Dispersion estimation
# c) Differential expression testing (Wald test
dds <- DESeq(dds)

#this gives the expected table output including log2foldchange, pvalue and padj
res <- results(dds, contrast =c("genotype", "mir802_floxed", "mir802_normal"), alpha=0.05)
res
res_shrink <- lfcShrink(dds, coef=2, type="apeglm")
res_shrink
#Adjusting the pvalues to less than 0.05
deseq2_adj <- subset(res, padj < 0.05)
view(deseq2_adj)

#this is the time stamp
time_end_deseq2 <- Sys.time()
time_deseq2 <- time_end_deseq2 - time_start_deseq2

#a quick digress to confirm my dds analysis with what I have from the original data base
#a brief comparison of the normal result in the data base with my analysis
res_df <- as.data.frame(res)    #this converts the object to a dataframe
head(res_df)
head(genes)
res_df <- merge(res_df, genes, by='row.names')
head(res_df)

#Here I want to check the outcome of nmy analyis with the database result for the log2fold change
#I have selected three genes to check, just to confirm the log2fold change value
genes_to_check <- c("Nlrplc-ps", "Sultlel", "Cxcl13")
res_df[res_df$Gene.Name %in% genes_to_check, ]    #The value obtained for log2fold change are very similar with a difference of +- 0.1
# Time the start of the edgeR analysis for your speed comparison
time_start_edger <- Sys.time()
#I performed filtering and then ran the TMM normalization, followed by dispersion estimation and the statistical tests.
# Ensure newcount columns are in the same order as sampleinfo rownames
counts_matrix <- newcount[, rownames(sample_info)]

# Create edgeR DGEList object
dge <- DGEList(counts = counts_matrix, group = sample_info$genotype)

# Filter out low-count genes
keep <- filterByExpr(dge, group = sample_info$genotype)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM normalization)
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~genotype, data = sample_info)
colnames(design)
design

#Estimate dispersion
dge <- estimateDisp(dge, design)
plotBCV(dge)

# Generate logCPM values (use prior.count to stabilize log)
logCPM <- cpm(dge, log = TRUE, prior.count = 2)
# Fit negative binomial GLM
fit <- glmFit(dge, design)

# Perform likelihood ratio test for KO vs WT
lrt <- glmLRT(fit, coef = 2)  # coef = 2 corresponds to condition KO, intercept vs KO

# Get top differentially expressed genes
edgeR_results <- topTags(lrt, n = Inf)$table
head(edgeR_results)
#Extract the final list of Differentially Expressed Genes (DEGs)
# Filter for genes with adjusted p-value < 0.05 (FDR is used by edgeR)
edgeR_adj<- subset(edgeR_results, FDR < 0.05)
edgeR_df <- edgeR_results
#Now that i have my data explored and sorted, I am going to start with some visualizations
#1 Calculate normalized counts
norm_counts <- counts(dds, normalized=TRUE)
# Take a random subset of genes for plotting clarity
set.seed(123)
subset_genes <- sample(rownames(norm_counts), 5000)
df_norm <- data.frame(
  Raw = log2(rowMeans(newcount[subset_genes, ]) + 1),
  Norm = log2(rowMeans(norm_counts[subset_genes, ]) + 1))
ggplot(df_norm, aes(x = Raw, y = Norm)) +
  geom_point(alpha=0.3, color="#E69F00") +  # orange — colorblind safe
  theme_minimal() +
  labs(title = "Raw vs Normalized Gene Counts", x = "Raw log2(Counts + 1)", y = "Normalized log2(Counts + 1)")

#2. MA plots of DESeq2 and edgeR
# MA Plot
res_df <- as.data.frame(res) %>%
  # Filter out rows with NA values (common for genes with zero counts)
  na.omit() %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"))

#Create the ggplot (MA Plot)
ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +  scale_x_log10() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  labs(
    x = "Mean of Normalized Counts (log10 Scale)",
    y = expression(log[2]~"Fold Change"),
    title = "MA Plot: Log2 Fold Change vs. Mean Expression") + theme_minimal()

plotMD(lrt, main = "edgeR: KnockOut vs Wildtype", ylim = c(-5, 5))

#3 Another important aspect is the sample clustering after normalization
# MDS Plot and PCA
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size=4) +
  scale_color_manual(values=c("mir802_normal"="green", "mir802_floxed"="steelblue")) + 
  theme_minimal() +
  labs(title = "PCA of Samples After rlog Normalization",
    x = paste0("PC1: ", percentVar[1], "% Variance"),
    y = paste0("PC2: ", percentVar[2], "% Variance")
  )
plotMDS(dge, labels = sample_info$genotype, col = as.numeric(sample_info$genotype), main="edgeR MDS Plot (Leading LogFC Dimensions)")

# Trying to see if i adjust my stat data, I can have similar answer for deseq and edger
common_genes <- intersect(
  rownames(subset(res_df, padj < 0.05)),
  rownames(subset(edgeR_df, FDR < 0.05))
)
length(common_genes)
head(common_genes)
#Comparing how Median Ratio vs TMM changes normalization, log2 fold-changes and which genes are called significant without adjusting
common <- intersect(rownames(res_df), rownames(edgeR_df))
length(common)
head(common)

#4Volcano report for DESeq2 and EdgeR
res_df <- as.data.frame(res)
res_df <- res_df %>%
  na.omit() %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha=0.5, size=1.2) +
  scale_color_manual(values=c(
    "Up"="red",    
    "Down"="blue", 
    "NS"="grey60"
  )) +
  theme_minimal() +
  labs(title = "DESeq2 Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted p-value")

#Volcano plot edgeR
edgeR_df <- as.data.frame(edgeR_results)
edgeR_df <- edgeR_df %>%
  mutate(
    Significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"))
ggplot(edgeR_df, aes(x = logFC, y = -log10(FDR), color = Significance)) +
  geom_point(alpha=0.5, size=1.2) +
  scale_color_manual(values=c(
    "Up"="red",
    "Down"="blue",
    "NS"="grey60"
  )) +
  theme_minimal() +
  labs(
    title = "edgeR Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 FDR")

#Record time and summarize
time_end_edger <- Sys.time()
time_edger <- time_end_edger - time_start_edger

# Get the gene names (IDs) for significant results
deseq2_genes <- rownames(deseq2_adj)
edgeR_genes <- rownames(edgeR_adj)
#5 Venn diagram
# Calculate the overlap
overlap_genes <- intersect(deseq2_genes, edgeR_genes)
view(overlap_genes)
# Create a list of the sets for the VennDiagram function
gene_sets <- list(
  "DESeq2 DEGs" = deseq2_genes,
  "edgeR DEGs" = edgeR_genes)

deseq2_genes <- rownames(subset(res_df, padj < 0.05))
edgeR_genes <- rownames(subset(edgeR_results, FDR < 0.05))

venn.diagram(
  x = list(DESeq2 = deseq2_genes, edgeR  = edgeR_genes),
  filename = "Venn_DESeq2_edgeR.png",
  fill = c("blue", "green"),
  alpha = 0.6,
  cex = 1.6,
  cat.cex = 1.4,
  lwd = 2)