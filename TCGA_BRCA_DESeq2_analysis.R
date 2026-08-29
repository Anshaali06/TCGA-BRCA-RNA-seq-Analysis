
# ============================================================
# 1. INSTALL REQUIRED PACKAGES
# ============================================================

install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks",
  "DESeq2",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

install.packages(c(
  "ggplot2",
  "ggrepel"
))

BiocManager::install(c(
  "ComplexHeatmap",
  "circlize"
))


# ============================================================
# 2. LOAD REQUIRED LIBRARIES
# ============================================================

library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)


# ============================================================
# 3. QUERY TCGA-BRCA DATA
# ============================================================

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

results <- getResults(query)

nrow(results)

table(results$sample_type)


# ============================================================
# 4. DOWNLOAD TCGA-BRCA DATA
# ============================================================

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 30
)


# ============================================================
# 5. PREPARE DOWNLOADED DATA
# ============================================================

data <- GDCprepare(query)


# ============================================================
# 6. EXTRACT COUNT MATRIX AND SAMPLE INFORMATION
# ============================================================

counts <- assay(
  data,
  "unstranded"
)

coldata <- as.data.frame(
  colData(data)
)


# ============================================================
# 7. DEFINE TUMOUR AND NORMAL CONDITIONS
# ============================================================

coldata$condition <- ifelse(
  coldata$sample_type == "Primary Tumor",
  "tumour",
  ifelse(
    coldata$sample_type == "Solid Tissue Normal",
    "normal",
    NA
  )
)

table(
  coldata$condition,
  useNA = "ifany"
)


# ============================================================
# 8. KEEP ONLY TUMOUR AND NORMAL SAMPLES
# ============================================================

keep <- !is.na(
  coldata$condition
)

counts <- counts[
  ,
  keep
]

coldata <- coldata[
  keep,
  ,
  drop = FALSE
]


# ============================================================
# 9. CREATE DESEQ2 DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)


# ============================================================
# 10. FILTER LOW-COUNT GENES
# ============================================================

dds <- dds[
  rowSums(counts(dds)) >= 10,
  ]

nrow(dds)


# ============================================================
# 11. RUN DESEQ2 ANALYSIS
# ============================================================

dds <- DESeq(dds)


# ============================================================
# 12. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

res <- results(
  dds,
  contrast = c(
    "condition",
    "tumour",
    "normal"
  )
)

res_df <- as.data.frame(res)

res_df$gene <- rownames(res_df)

res_df <- res_df[
  order(res_df$padj),
  ]

head(res_df)


# ============================================================
# 13. FILTER SIGNIFICANT DIFFERENTIALLY EXPRESSED GENES
# ============================================================

sig_genes <- res_df[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  abs(res_df$log2FoldChange) > 1,
]

nrow(sig_genes)


# ============================================================
# 14. IDENTIFY UPREGULATED GENES
# ============================================================

upregulated <- sig_genes[
  sig_genes$log2FoldChange > 1,
  ]

nrow(upregulated)


# ============================================================
# 15. IDENTIFY DOWNREGULATED GENES
# ============================================================

downregulated <- sig_genes[
  sig_genes$log2FoldChange < -1,
  ]

nrow(downregulated)
