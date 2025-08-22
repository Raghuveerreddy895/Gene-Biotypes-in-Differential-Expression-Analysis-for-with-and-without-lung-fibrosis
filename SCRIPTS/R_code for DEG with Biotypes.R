# ================================
# RNA-seq DESeq2 Workflow with lncRNA extraction
# ================================

# Install once (skip if already installed)
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano", "biomaRt", "AnnotationDbi"))
# install.packages(c("ggplot2", "RColorBrewer"))

# Libraries
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(biomaRt)

# ================================
# Load Data
# ================================
counts <- read.csv("C:/Users/raghu/Downloads/Other_folder/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv",
                   row.names = 1)
coldata <- read.csv("C:/Users/raghu/Downloads/Other_folder/meta.csv", row.names = 1)

# Condition (Normal vs Disease)
coldata$condition <- ifelse(grepl("Normal", coldata$Sample_characteristics_ch1),
                            "Normal", "Disease")
coldata$condition <- factor(coldata$condition)

# Reorder coldata to match counts
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(coldata)))

# ================================
# DESeq2 Analysis
# ================================
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]   # filter low counts
dds <- DESeq(dds)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
res_df$ensembl_gene_id <- sub("\\.\\d+$", "", res_df$ensembl_gene_id)  # remove version suffix

# ================================
# Annotation (Ensembl biomaRt)
# ================================
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = res_df$ensembl_gene_id,
  mart = mart
)

res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id", all.x = TRUE)

# ================================
# Classification
# ================================
cnts <- counts(dds)

# 1. Upregulated
upregulated <- res_annotated[res_annotated$log2FoldChange > 1 & res_annotated$padj < 0.05, ]
upregulated$Status <- "Upregulated"

# 2. Downregulated
downregulated <- res_annotated[res_annotated$log2FoldChange < -1 & res_annotated$padj < 0.05, ]
downregulated$Status <- "Downregulated"


# 4. Not Significant
expressed_genes <- rownames(cnts)[rowSums(cnts) > 0]
not_sig_genes <- setdiff(expressed_genes, c(upregulated$ensembl_gene_id,
                                            downregulated$ensembl_gene_id,
                                            unexpressed$ensembl_gene_id))
not_significant <- res_annotated[res_annotated$ensembl_gene_id %in% not_sig_genes, ]
not_significant$Status <- "Not Significant"

# Combine
all_genes_combined <- rbind(upregulated, downregulated, unexpressed, not_significant)
write.csv(all_genes_combined, "DEGs_up_down_unexpressed_notsig.csv", row.names = FALSE)

# ================================
# Combine all 4 categories
# ================================
all_genes_combined <- rbind(upregulated, downregulated, unexpressed, not_significant)

# Add Regulated_Biotype column (Status + gene_biotype)
all_genes_combined$Regulated_Biotype <- paste(all_genes_combined$Status,
                                              all_genes_combined$gene_biotype,
                                              sep = "_")

# Save
write.csv(all_genes_combined, "DEGs_up_down_unexpressed_notsig.csv", row.names = FALSE)

# ================================
# Extract lncRNAs only
# ================================
lnc_types <- c("lncRNA", "antisense", "processed_transcript",
               "sense_intronic", "sense_overlapping",
               "3prime_overlapping_ncRNA", "macro_lncRNA",
               "bidirectional_promoter_lncRNA")

lnc_all <- subset(all_genes_combined, gene_biotype %in% lnc_types)

# Save lncRNA file
write.csv(lnc_all, "lncRNA_up_down_unexpressed_notsig.csv", row.names = FALSE)

# Also split into separate regulated lncRNA files (optional)
lnc_up <- subset(lnc_all, Status == "Upregulated")
lnc_down <- subset(lnc_all, Status == "Downregulated")
lnc_unexp <- subset(lnc_all, Status == "Unexpressed")
lnc_notsig <- subset(lnc_all, Status == "Not Significant")

write.csv(lnc_up, "lncRNA_upregulated.csv", row.names = FALSE)
write.csv(lnc_down, "lncRNA_downregulated.csv", row.names = FALSE)
write.csv(lnc_unexp, "lncRNA_unexpressed.csv", row.names = FALSE)
write.csv(lnc_notsig, "lncRNA_not_significant.csv", row.names = FALSE)


# ================================
# Extract lncRNAs
# ================================
lnc_types <- c("lncRNA", "antisense", "processed_transcript",
               "sense_intronic", "sense_overlapping",
               "3prime_overlapping_ncRNA", "macro_lncRNA",
               "bidirectional_promoter_lncRNA")

lnc_all <- subset(all_genes_combined, gene_biotype %in% lnc_types)
write.csv(lnc_all, "lncRNA_up_down_unexpressed_notsig.csv", row.names = FALSE)

# ================================
# Plots
# ================================
EnhancedVolcano(res_annotated,
                lab = res_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano Plot',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.5)

vsd <- vst(dds, blind = FALSE)
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:50]
pheatmap(assay(vsd)[select, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         annotation_col = coldata)

plotPCA(vsd, intgroup = "condition")

