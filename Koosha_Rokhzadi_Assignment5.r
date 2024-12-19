#*************************************************************************************************************************************
# Differential Gene Expression Analysis: edgeR vs DESeq2
# Course: Software Tools (BINF*6210)
# Koosha Rokhzadi  
# Assignment: 5
# the data I will be using comes from the airway dataset, which is part of the airway Bioconductor package.
# This dataset represents RNA-Seq experiments in airway smooth muscle cells and was described in:
# Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Kher SS, Ladd-Acosta C, Yang IV, Bonner D, Soto-Quiros ME, et al.
# (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates
# Cytokine Function in Airway Smooth Muscle Cells. 
#************************************** Install and Load Necessary Libraries ********************************************************

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(
#   c("edgeR", "DESeq2", "airway", "VennDiagram", "ggplot2", 
#     "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"),
#   force = TRUE
# )
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("edgeR", force = TRUE)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("limma", force = TRUE)


library(edgeR)           # For edgeR RNA-Seq differential expression analysis.
library(DESeq2)          # For DESeq2 RNA-Seq differential expression analysis.
library(airway)          # Dataset used for this analysis.
library(VennDiagram)     # To create Venn diagrams to compare results.
library(ggplot2)         # For high-quality plots.
library(EnhancedVolcano) # To create volcano plots for visualization of DEGs.
library(clusterProfiler) # For GO enrichment analysis of gene lists.
library(org.Hs.eg.db)    # Human annotation database for GO analysis.
library(pheatmap)        # For heatmap visualization.

#******************************************* Loading airway Data ********************************************************************************
# The following codes loads the airway dataset, a (RangedSummarizedExperiment) object containing RNA-Seq count data and metadata.
# This data provides the input data (raw counts and metadata) for downstream differential expression analysis.
data(airway)
airway$dex <- relevel(airway$dex, ref = "untrt") # Set "untreated" as the reference condition for comparison which Ensures that statistical comparisons are made relative to the untreated condition.
counts <- assay(airway) # Extract the counts matrix for differential analysis This .Prepares the data for analysis by separating counts and metadata for convenience.
metadata <- colData(airway) # Extract sample metadata for grouping and experimental conditions.

#******************************************* Missing Data Handling ********************************************************************************
# Verify if there are any missing values in the dataset
missing_counts <- sum(is.na(counts))  # Check for missing values in the counts matrix
missing_metadata <- sum(is.na(metadata))  # Check for missing values in metadata

if (missing_counts == 0 & missing_metadata == 0) {
  print("No missing values were detected in the dataset.")
} else {
  print(paste("Missing values detected: Counts -", missing_counts, ", Metadata -", missing_metadata))
}
# *Observation: No missing values were detected during data exploration, so no imputation was required.

#************************************* Data Exploration and Quality Control ***********************************************************************
# Summarize count data
# The following codes ensures there are no anomalies in the data (e.g., too many zeros or negative values).
summary(counts) # Summarize counts data to understand its distribution.
# *Observation: The summary shows most genes have zero counts (median = 0), typical of RNA-Seq data.
# Variability in mean and max values highlights the need for normalization and filtering.*

dim(counts) # Get the dimensions of the counts matrix (genes x samples). Dimensions help verify the number of genes and samples

# Data preparation
# Calculate library sizes (total counts for each sample) for visualization
library_sizes <- colSums(counts)  # Calculate the total number of reads per sample. Library sizes help identify samples with insufficient sequencing depth.
sample_names <- colnames(counts)  # Extract sample names.

#*********************************************** Potential Biases in Data ************************************************************************
# Compute summary statistics for library sizes
library_sizes_summary <- summary(library_sizes)
print(library_sizes_summary)

# Highlight sequencing depth variability
cat("Observation: Library sizes vary between samples, ranging from",
    min(library_sizes), "to", max(library_sizes),
    "reads. This variability can introduce biases in differential expression analysis,
     which will be addressed through normalization methods (TMM in edgeR and size factors in DESeq2).")

# *Observation: Observation: Library sizes vary between samples, ranging from 15163415 to 30818215 reads.
# This variability can introduce biases in differential expression analysis,
# which will be addressed through normalization methods (TMM in edgeR and size factors in DESeq2).*

#******************************************** Create a barplot of library sizes ******************************************************************

# Create a data frame for library size visualization
# This code organizes the data for visualization with ggplot2
library_sizes_df <- data.frame(Sample = sample_names, LibrarySize = library_sizes)
# *Observation: The library_sizes_df data frame contains two columns. 1:Sample: Lists the sample names (e.g., SRR1039508, SRR1039509). 2:Displays the total number of raw counts (reads) for each sample.

# Create a barplot of library sizes using ggplot2
# Barplot ensures library sizes are roughly equal across samples, indicating no outliers or poorly sequenced samples.
ggplot(library_sizes_df, aes(x = Sample, y = LibrarySize)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
  geom_text(aes(label = scales::comma(LibrarySize)), vjust = -0.5, size = 3.5) +  # Add numeric labels on bars
  theme_minimal() +  # Apply minimal theme for cleaner appearance
  labs(
    title = "Library Sizes for RNA-Seq Samples",
    x = "Sample Name",
    y = "Library Size (Total Counts)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and style title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  scale_y_continuous(labels = scales::comma)  # Format y-axis numbers with commas

# *Observation: This barplot displays the total library sizes (raw read counts) for each RNA-Seq sample,
# ranging from ~15M (SRR1039513) to ~31M (SRR1039517), highlighting variability in sequencing depth across samples.
# Such variation is common in RNA-Seq experiments and underscores the importance of normalization methods,
# like TMM (edgeR) or size factor estimation (DESeq2), to adjust for differences and enable reliable comparisons between samples in downstream analyses.*

#************************************** Log-transform and Filter Lowly Expressed Genes ********************************************************
# Log-transform the counts for visualization of distribution
# This code examines data distribution across samples and ensures no extreme variability that could affect downstream analysis.
log_counts <- log2(counts + 1)  # Apply log2 transformation to stabilize variance and avoid zeros
boxplot(log_counts, las = 2, col = "lightblue", 
        main = "Log-transformed Counts Distribution",
        xlab = "Samples", ylab = "Log2(Counts)")
# *Observation: This boxplot shows the distribution of log2-transformed RNA-Seq counts for each sample.
# The log transformation reduces variability and highlights trends across samples, making data more comparable.
# Each sample exhibits a similar distribution, with most counts concentrated below a log2 value of ~5.
# Outliers (points above whiskers) represent highly expressed genes. The uniformity across samples indicates
# consistent sequencing and data quality, preparing the dataset for normalization and differential analysis.*

#************************************************ Filter Lowly Expressed Genes ****************************************************************

# This filter retains only informative genes, reducing noise and improving statistical power.
keep <- rowSums(counts >= 10) >= 2  # Retain genes expressed in at least 2 samples with >=10 counts
filtered_counts <- counts[keep, ] # Subset the filtered counts matrix
# *Observation: This table shows the raw RNA-Seq read counts for each gene (rows, identified by their Ensembl IDs)
# across samples (columns labeled by sample IDs).
dim(filtered_counts) # Check the dimensions of the filtered matrix

#******************************************* Differential Expression Analysis with edgeR and DESeq2 ********************************************
# Perform edgeR analysis
# This analysis will prepare normalized data for statistical modeling and identifies differentially
# expressed genes (DEGs) in the dataset and also Focuses on statistically significant DEGs.
dge <- DGEList(counts = filtered_counts, group = metadata$dex) # Create DGEList object for edgeR
dge <- calcNormFactors(dge)             # Apply TMM normalization to account for library size differences
design <- model.matrix(~ metadata$dex)  # Create a design matrix for the differential analysis
dge <- estimateDisp(dge, design)        # Estimate dispersion for modeling
fit <- glmQLFit(dge, design)            # Fit a quasi-likelihood model
edgeR_results <- glmQLFTest(fit)         # Perform the statistical test for DEGs
edgeR_table <- topTags(edgeR_results, n = Inf)$table  #Extract all genes and associated statistics

# *Observation: The table (edgeR_table) shows the **edgeR results** for differential expression, including log2 fold changes (`logFC`),
# overall expression levels (`logCPM`), and statistical significance (`FDR`). Genes with FDR < 0.05 are significant,
# with positive `logFC` indicating upregulation and negative `logFC` indicating downregulation in treated samples.
# These results highlight key differentially expressed genes for further analysis.*


# Filter significant genes from edgeR results
edgeR_significant <- edgeR_table[edgeR_table$FDR < 0.05, ]  # Keep genes with FDR < 0.05.


# Perform DESeq2 analysis
# Normalizes counts and performs differential expression analysis using DESeq2
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = metadata,
                              design = ~ dex)  # Create DESeq2 dataset.
dds <- DESeq(dds)   # Perform DESeq2 normalization and differential analysis
deseq_results <- results(dds, alpha = 0.05)  # Extract results with significance threshold of 0.05 (Results table with log fold changes, p-values, and FDR)
deseq_significant <- as.data.frame(deseq_results[!is.na(deseq_results$padj) & 
                                                   deseq_results$padj < 0.05, ])  # Filters significant DEGs (FDR < 0.05). # Identifies DEGs using DESeq2 for comparison with edgeR results
# *Observation: This table shows the **DESeq2 results** for significant genes, with columns for mean expression (`baseMean`), 
# log2 fold change (`log2FoldChange`), standard error (`lfcSE`), and statistical significance (`pvalue` and `padj`).
# Genes with `padj` < 0.05 are significant, where positive `log2FoldChange` indicates upregulation and negative values
# indicate downregulation. These results identify key genes for further biological interpretation or validation.*

#************************************************** Compare Results: Venn Diagram ********************************************************************
# Compare results using a Venn diagram
# Overlap in significant genes
overlap <- intersect(rownames(edgeR_significant), rownames(deseq_significant)) # Find overlapping genes.

# Venn Diagram
# Creates a Venn diagram showing overlapping significant genes between edgeR and DESeq2.
# This compares results from the two methods to assess consistency
venn.diagram(
  x = list(
    edgeR = rownames(edgeR_significant), 
    DESeq2 = rownames(deseq_significant)
  ),
  filename = "venn_diagram.png",  
  category.names = c("edgeR Significant Genes", "DESeq2 Significant Genes"), # Descriptive category names
  fill = c("blue", "red"),               # Colors for the sets
  alpha = 0.7,                           # Adjust transparency
  cex = 3,                               # Increase font size for counts
  cat.cex = 2.5,                         # Increase font size for category names
  cat.fontface = "bold",                 # Bold category labels
  cat.col = c("blue", "red"),            # Match category label colors to the sets
  cat.pos = c(-20, 20),                  # Adjust label positions
  cat.dist = c(0.06, 0.06),              # Adjust distance of category labels from circles
  margin = 0.1,                          # Add margin for spacing
  lwd = 3,                               # Increase line width for circle outlines
  col = c("blue", "red"),                # Outline colors for the circles
  main = "Comparison of Significant Genes Between edgeR and DESeq2", # Add a descriptive title
  main.cex = 3,                          # Increase font size for the title
  main.fontface = "bold",                # Bold the title
  height = 4000,                         # Increase height of the image (in pixels)
  width = 4000,                          # Increase width of the image (in pixels)
  resolution = 300                       # Set resolution to 300 dpi for high quality
)
# *Observation: The Venn diagram compares significant genes identified by edgeR and DESeq2, showing 2,028 shared genes,
# 758 unique to DESeq2, and 28 unique to edgeR. This highlights strong overlap, indicating robustness in shared genes.*

#********************************************* Visualize Results: Volcano Plots *******************************************************************
# Volcano Plot for edgeR
# The volcano plot highlighting significant genes from edgeR.
# The chosen fold change cutoff (logFC > |1.5|) highlights biologically meaningful changes,
# while the p-value cutoff (FDR < 0.05) ensures statistical significance after correcting for multiple testing.
EnhancedVolcano(edgeR_table,
                lab = rownames(edgeR_table),
                x = 'logFC',
                y = 'FDR',
                title = 'Volcano Plot (edgeR)',  # Title 
                subtitle = NULL,              
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = c(ifelse(edgeR_table$FDR < 0.05, 3, 1)),
                labSize = 3,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'))

# *Observation: This volcano plot shows the differential expression results from edgeR, where each point represents a gene.
# The x-axis indicates log2 fold change (upregulated on the right, downregulated on the left), and the y-axis represents the -log10 p-value,
# highlighting statistical significance. Red points denote genes with significant p-values (FDR < 0.05) and large fold changes,
# while blue and grey points represent less significant or non-significant genes. This visualization identifies key genes for
# further analysis based on fold change and significance (if possible 🙆).*

# Volcano Plot for DESeq2
# The volcano plot highlighting significant genes from DESeq2
# Similarly, log2 fold changes greater than |1.5| capture biologically relevant changes,
# and adjusted p-value (padj) < 0.05 ensures significance under stringent criteria.

EnhancedVolcano(deseq_results,
                lab = rownames(deseq_results),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano Plot (DESeq2)',  # title
                subtitle = NULL,                
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = c(ifelse(deseq_results$padj < 0.05, 3, 1)),
                labSize = 3,
                col = c('grey30', 'purple', 'gold', 'darkorange'))

# *Observation: This volcano plot visualizes the differential expression results from DESeq2. Each point represents a gene,
# with the x-axis showing the log2 fold change (positive for upregulated and negative for downregulated genes)
# and the y-axis representing -log10 p-values. Yellow points highlight genes with significant p-values (adjusted p-value < 0.05)
# and large fold changes, while grey and purple points represent non-significant genes or those with only large fold changes.
# This plot identifies key differentially expressed genes for further analysis based on both statistical significance and biological relevance (if possible 🙆).*

#*************************************************** Gene Ontology (GO) Enrichment ***********************************************************************
# Purpose: To derive biological insights by annotating significant DEGs identified by edgeR and DESeq2
# with enriched Gene Ontology terms for biological processes.
# This performs Gene Ontology enrichment analysis on significant genes.
annotate_genes <- function(gene_ids) {
  enrichGO(
    gene = gene_ids,               # The input genes for enrichment analysis
    OrgDb = org.Hs.eg.db,          # The annotation database for Homo sapiens (human)
    keyType = "ENSEMBL",           # Specifies the type of gene IDs (here, ENSEMBL IDs)
    ont = "BP",                    # Focuses on the Biological Process (BP) ontology
    pAdjustMethod = "BH",          # Adjusts p-values using the Benjamini-Hochberg method
    pvalueCutoff = 0.05            # Filters results to show only significant GO terms (adjusted p-value < 0.05)
  )
}

# Annotate edgeR significant genes
ego_edgeR <- annotate_genes(rownames(edgeR_significant))
print(head(ego_edgeR))
# *Observation: This table presents the Gene Ontology (GO) enrichment results for edgeR significant genes, highlighting enriched
# biological processes such as "connective tissue development" and "response to peptide hormone.
# " Each term shows the proportion of significant genes (GeneRatio), enrichment strength (FoldEnrichment),
# and statistical significance (p.adjust < 0.05). These results reveal key pathways and processes associated with the experimental conditions,
# offering insights into the biological relevance of the differentially expressed genes.*

# Annotate DESeq2 significant genes
ego_DESeq2 <- annotate_genes(rownames(deseq_significant))
print(head(ego_DESeq2))
# *Observation: This table summarizes the Gene Ontology (GO) enrichment results for significant genes identified by DESeq2,
# highlighting biological processes associated with differentially expressed genes. Key enriched terms
# include "connective tissue development" (80 genes), "extracellular matrix organization" (88 genes), and "response to peptide hormone" (99 genes).
# Each term is associated with the proportion of significant genes (GeneRatio), enrichment strength (FoldEnrichment), and adjusted p-values (p.adjust < 0.05),
# confirming their statistical significance. These results provide insights into the functional roles of DESeq2-identified genes,
# emphasizing pathways related to structural organization and hormone responses under experimental conditions.*

# Define the function for Gene Ontology (GO) enrichment
annotate_genes <- function(gene_ids) {
  enrichGO(
    gene = gene_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",  # Benjamini-Hochberg for multiple testing correction
    pvalueCutoff = 0.05  # Only significant terms
  )
}

# Perform GO enrichment for edgeR significant genes
if (!is.null(edgeR_significant) && nrow(edgeR_significant) > 0) {
  ego_edgeR <- annotate_genes(rownames(edgeR_significant))
  print("Top GO terms for edgeR significant genes:")
  print(head(ego_edgeR))
  
  # Visualize results if enrichment is found
  if (!is.null(ego_edgeR) && nrow(as.data.frame(ego_edgeR)) > 0) {
    barplot(ego_edgeR, showCategory = 10, title = "Top 10 GO Terms (edgeR - Biological Process)")
    dotplot(ego_edgeR, showCategory = 10, title = "GO Enrichment Dotplot (edgeR - Biological Process)")
  } else {
    print("No significant GO terms found for edgeR significant genes.")
  }
}

# *Observation: This *dot plot* visualizes the *Gene Ontology (GO) enrichment analysis results* for biological processes identified using edgeR.
# The x-axis represents the **GeneRatio** (proportion of significant genes associated with each process), and the y-axis lists enriched GO terms
# such as "response to peptide hormone," "axonogenesis," and "connective tissue development." Dot size corresponds to the number of genes (`Count`)
# nvolved in each process, while the color gradient represents the adjusted p-value (`p.adjust`), with red indicating more significant enrichment.
# This plot highlights processes related to tissue development, signaling, and cellular organization as key biological pathways impacted in the dataset.*

# Perform GO enrichment for DESeq2 significant genes
if (!is.null(deseq_significant) && nrow(deseq_significant) > 0) {
  ego_DESeq2 <- annotate_genes(rownames(deseq_significant))
  print("Top GO terms for DESeq2 significant genes:")
  print(head(ego_DESeq2))
  
  # Visualize results if enrichment is found
  if (!is.null(ego_DESeq2) && nrow(as.data.frame(ego_DESeq2)) > 0) {
    barplot(ego_DESeq2, showCategory = 10, title = "Top 10 GO Terms (DESeq2 - Biological Process)")
    dotplot(ego_DESeq2, showCategory = 10, title = "GO Enrichment Dotplot (DESeq2 - Biological Process)")
  } else {
    print("No significant GO terms found for DESeq2 significant genes.")
  }
}

# *Observation: This dot plot illustrates the results of *Gene Ontology (GO) enrichment analysis* for biological processes based on DESeq2 results.
# The x-axis represents the *GeneRatio* (proportion of significant genes associated with each biological process), while the y-axis lists enriched GO terms,
# such as "ossification," "response to peptide hormone," and "cartilage development." Dot sizes correspond to the number of genes (`Count`) involved in each process,
# and the color gradient indicates the adjusted p-value (`p.adjust`), with darker red denoting higher significance.
# This visualization emphasizes critical biological processes related to tissue development, extracellular structure organization, and signaling pathways influenced in the dataset.*

#************************************************ Heatmap for Both edgeR and DESeq2 ***************************************************************************
# Extract significant genes from edgeR results
edgeR_DEGs <- edgeR_significant[edgeR_significant$FDR < 0.05, ]  # Filter significant DEGs
dge_normalized <- cpm(dge, log = TRUE)  # Normalize counts using log-transformed CPM
heatmap_data_edgeR <- dge_normalized[rownames(edgeR_DEGs), ]  # Subset normalized data for significant genes

# Create a clustered heatmap for significant genes from edgeR
pheatmap(
  heatmap_data_edgeR,
  cluster_rows = TRUE,  # Cluster rows (genes) based on expression patterns
  cluster_cols = TRUE,  # Cluster columns (samples) based on similarity
  main = "Clustered Heatmap of Significant Genes (edgeR)",
  fontsize_row = 6,  # Adjust font size for rows (gene names)
  fontsize_col = 8,  # Adjust font size for columns (samples)
  color = colorRampPalette(c("blue", "white", "red"))(50)  # Define color gradient
)

# Explanation:
# The heatmap clusters significant genes (rows) and samples (columns) based on normalized expression levels. 
# This visualization helps identify patterns of co-expressed genes and explore functional groupings or sample similarities.

# *Observation: The heatmap visualizes the normalized expression levels of significant differentially expressed genes (DEGs) identified using edgeR,
# with rows representing genes and columns representing samples. A clear color gradient shows high expression in red and low expression in blue,
# highlighting transcriptional variability across samples. Clustering reveals co-expressed gene groups and sample similarities,
# likely reflecting biological differences between conditions (e.g., treated vs. untreated). Distinct expression patterns,
# such as upregulation or downregulation in specific samples, indicate treatment-related transcriptional changes. Overall,
# the clustering suggests functional relationships among DEGs and consistent biological responses across replicates.

# Clustering DEGs by Functional Categories: DESeq2
# Extract significant genes from DESeq2 results
DESeq2_DEGs <- deseq_significant[!is.na(deseq_significant$padj) & deseq_significant$padj < 0.05, ]  # Filter significant DEGs
dds_normalized <- rlog(dds)  # Apply regularized log transformation for normalized counts
heatmap_data_DESeq2 <- assay(dds_normalized)[rownames(DESeq2_DEGs), ]  # Subset normalized data for significant genes

# Create a clustered heatmap for significant genes from DESeq2
pheatmap(
  heatmap_data_DESeq2,
  cluster_rows = TRUE,  # Cluster rows (genes) based on expression patterns
  cluster_cols = TRUE,  # Cluster columns (samples) based on similarity
  main = "Clustered Heatmap of Significant Genes (DESeq2)",
  fontsize_row = 6,  # Adjust font size for rows (gene names)
  fontsize_col = 8,  # Adjust font size for columns (samples)
  color = colorRampPalette(c("purple", "white", "orange"))(50)  # Define color gradient
)

# Explanation:
# The heatmap clusters significant genes (rows) and samples (columns) based on normalized expression levels from DESeq2.
# It visualizes patterns of co-expression among significant DEGs, offering insights into functional groupings and biological relevance.
# *Observation: The heatmap visualizes the normalized expression levels of significant differentially expressed genes (DEGs) identified using DESeq2,
# with rows representing genes and columns representing RNA-Seq samples. The color gradient ranges from purple (low expression) to orange (high expression), 
# indicating transcriptional variability across samples. Clustering reveals distinct co-expression patterns among genes and similarities between samples,
# likely reflecting differences between experimental conditions (e.g., treated vs. untreated). Genes with similar expression profiles are grouped together,
# suggesting shared biological functions, while sample clustering highlights consistent responses across replicates.
# This visualization provides insights into both the biological processes underlying DEGs and the treatment-related transcriptional changes.

#**************************************************** Execution Time Comparison *****************************************************************************************
# Execution Time Comparison
# Purpose: Evaluate computational efficiency of edgeR and DESeq2.
start_edgeR <- Sys.time()
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
glmQLFTest(fit)
end_edgeR <- Sys.time()

start_DESeq2 <- Sys.time()
dds <- DESeq(dds)
results(dds)
end_DESeq2 <- Sys.time()

print(paste("Execution time for edgeR:", end_edgeR - start_edgeR))
#"Execution time for edgeR: 16.9364078044891"
print(paste("Execution time for DESeq2:", end_DESeq2 - start_DESeq2))
#"Execution time for DESeq2: 10.1060450077057"

# Save Outputs
# Saves the significant genes from edgeR and DESeq2 to CSV files
write.csv(edgeR_significant, "edgeR_significant_genes.csv")
write.csv(deseq_significant, "DESeq2_significant_genes.csv")

# Explanation: The execution time comparison highlights the computational efficiency of edgeR and DESeq2,
# providing practical insights into their performance in differential gene expression (DGE) analysis.
# DESeq2 completed the analysis in approximately 10.1 seconds, faster than edgeR’s 16.9 seconds,
# suggesting that DESeq2 may be more suitable for larger datasets or scenarios requiring quicker turnaround times.
# While edgeR’s slightly longer runtime could be attributed to its robust handling of small or variable datasets,
# this comparison demonstrates the importance of balancing computational efficiency with analytical needs.
# Including execution time in the evaluation ensures a comprehensive assessment of these tools,
# addressing both biological relevance and practical usability in RNA-Seq workflows.


print(" Analysis completed and Outputs saved (: ")  


