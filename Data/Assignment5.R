# script to perform differential gene expression analysis using DESeq2 and edgeR package with reference to Tong et al., 2020
#I couldn't find my way around NCBI GEO, so I settled for EMBL-EBI
install.packages("pheatmap")

# load all libraries
library(DESeq2)
library(tidyverse)
library(edgeR)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(pheatmap)
#Now let us prepare the data
#download the raw count data from ebi

count_matrix <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10411/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(count_matrix)
class (count_matrix)
dim(count_matrix)
head(count_matrix[, 1:4])
colnames(count_matrix)

#download experimental data
sample_info <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10411/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(sample_info)
class(sample_info)
colnames(sample_info)


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
sample_info$Genotype <- factor(sample_info$Genotype, levels = c("mir802_normal", "mir802_floxed"))
sample_info$Genotype


# making sure the row names in sample_info matches to column names in count_matrix
all(colnames(newcount) %in% rownames(sample_info))

# are they in the same order?
all(colnames(newcount) == rownames(sample_info))






sampleTable <- data.frame(
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
# 1. Create the DESeqDataSet object
# genotype defines the comparison groups (wildtype vs knockout)
ncol(newcount) == nrow(sample_info)
ncol(newcount)
nrow(sample_info)
colnames(count_matrix)
rownames(sample_info)
#I am adding a time start of the DESeq2 analysis to check speed according to tong et al 2020
time_start_deseq2 <- Sys.time()
dds <- DESeqDataSetFromMatrix(countData = newcount, colData=sample_info, design=~Genotype)
#Keep only genes with at least 10 reads total
dds <- dds[rowSums(counts(dds))>10,]

# The below function performs:
# a) Normalization (Median of Ratios)
# b) Dispersion estimation
# c) Differential expression testing (Wald test
dds <- DESeq(dds)

res <- results(dds, contrast =c("Genotype", "mir802_floxed", "mir802_normal"), alpha=1e-5)
res
view(res)
deseq2_sig <- subset(res, padj < 0.05)
view(deseq2_sig)

time_end_deseq2 <- Sys.time()
time_deseq2 <- time_end_deseq2 - time_start_deseq2


#a brief comparison of the normal result in the data base with my analysis
res_df <- as.data.frame(res)    #this converts the object to a dataframe
head(res_df)
head(genes)
res_df <- merge(res_df, genes, by='row.names')
head(res_df)

genes_to_check <- c("Nlrplc-ps", "Catspere1", "Gm13736")
res_df[res_df$Gene.Name %in% genes_to_check, ]

#Some visualization
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


# Time the start of the edgeR analysis for your speed comparison
time_start_edger <- Sys.time()
library(edgeR)

#I performed filtering and then ran the TMM normalization, followed by dispersion estimation and the statistical tests.
# Ensure countMatrix columns are in the same order as sampleTable rownames
count_matrix <- count_matrix[, rownames(sampleTable)]

# Create edgeR DGEList object
dge <- DGEList(counts = count_matrix, group = sampleTable$genotype)

# Filter out low-count genes
keep <- filterByExpr(dge, group = sampleTable$genotype)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM normalization)
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~genotype, data = sampleTable)
colnames(design)
design

#Estimate dispersion
dge <- estimateDisp(dge, design)
plotBCV(dge)

# Fit negative binomial GLM
fit <- glmFit(dge, design)

# Perform likelihood ratio test for KO vs WT
lrt <- glmLRT(fit, coef = 2)  # coef = 2 corresponds to condition KO, intercept vs KO

# Get top differentially expressed genes
edgeR_results <- topTags(lrt, n = Inf)$table
head(edgeR_results)
#Extract the final list of Differentially Expressed Genes (DEGs)
# Filter for genes with adjusted p-value < 0.05 (FDR is used by edgeR)
edgeR_sig<- subset(edgeR_results, FDR < 0.05)


head(res_edger_sig)

# MA Plot
plotMD(lrt, main = "edgeR: KnockOut vs Wildtype", ylim = c(-5, 5))

# MDS Plot (similar to PCA)
plotMDS(dge, labels = sampleTable$genotype, col = as.numeric(sampleTable$genotype))

common_genes <- intersect(
  rownames(subset(res_df, padj < 0.05)),
  rownames(subset(edgeR_results, FDR < 0.05))
)
length(common_genes)
head(common_genes)

edgeR_df <- edgeR_results
edgeR_df$gene <- rownames(edgeR_df)
ggplot(edgeR_df, aes(x = logFC, y = -log10(FDR))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  ggtitle("Volcano plot (edgeR glmLRT)") +
  xlab("log2 Fold Change (KO vs WT)") + ylab("-log10(FDR)")

library(ggplot2)

# edgeR
edgeR_df <- as.data.frame(edgeR_results)
ggplot(edgeR_df, aes(x=logFC, y=-log10(FDR))) +
  geom_point(alpha=0.4) +
  theme_minimal() +
  ggtitle("Volcano Plot - edgeR")

# DESeq2
res_df <- as.data.frame(res)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) +
  theme_minimal() +
  ggtitle("Volcano Plot - DESeq2")


topGenes <- head(rownames(edgeR_results), 50)

mat <- assay(rld)[topGenes, ]
pheatmap(mat, scale="row", 
         main="Top 50 Differentially Expressed Genes in edgeR")

topGenes2 <- head(rownames(res_df), 50)

mat <- assay(rld)[topGenes, ]
pheatmap(mat, scale="row", 
         main="Top 50 Differentially Expressed Genes in deseq")



# rlog transformation for visualization
rld <- rlog(dds, blind=FALSE)

# PCA plot
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples (rlog normalized)") +
  theme_minimal()

#Record time and summarize
time_end_edger <- Sys.time()
time_edger <- time_end_edger - time_start_edger


# Get the gene names (IDs) for significant results
deseq2_genes <- rownames(deseq2_sig)
edgeR_genes <- rownames(edgeR_sig)

# Calculate the overlap
overlap_genes <- intersect(deseq2_genes, edgeR_genes)
view(overlap_genes)

# Create a list of the sets for the VennDiagram function
gene_sets <- list(
  "DESeq2 DEGs" = deseq2_genes,
  "edgeR DEGs" = edgeR_genes
)
# Generate and save the Venn Diagram
venn.diagram(
  x = gene_sets,
  category.names = names(gene_sets), # Use the descriptive names for the labels
  filename = 'DGE_Overlap_Venn_Diagram.png', # Name of the output file
  output = TRUE, # Saves the file
  
  # --- Aesthetics for the circles ---
  lwd = 2,
  lty = 'solid',
  fill = c("skyblue", "pink"), # Colors for the two circles
  alpha = 0.6,
  
  # --- Aesthetics for the numbers/labels ---
  cex = 1.5, # Size of the numbers
  cat.cex = 1.2, # Size of the category labels (DESeq2 DEGs, edgeR DEGs)
  cat.pos = c(-27, 27), # Position of the category labels (adjust for best fit)
  ext.percent = 5, # Extends the labels slightly outside the circles
  margin = 0.05 # Margin around the plot
)
