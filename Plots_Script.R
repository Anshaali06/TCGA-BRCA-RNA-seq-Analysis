# ============================================================
#  PREPARE DATA FOR HEATMAP
# ============================================================

# Take the top 50 significant genes
top_genes <- head(
  sig_genes[order(sig_genes$padj), ],
  50
)

# Normalized counts
norm_counts <- DESeq2::counts(
  dds,
  normalized = TRUE
)

# Select top genes
heatmap_data <- norm_counts[
  rownames(top_genes),
]

# Log2 transformation
heatmap_data <- log2(
  heatmap_data + 1
)

# Row-wise Z-score
heatmap_data <- t(
  scale(
    t(heatmap_data)
  )
)

# Remove rows containing NA values
heatmap_data <- heatmap_data[
  complete.cases(heatmap_data),
]

# Sample group annotation
annotation_col <- data.frame(
  Condition = coldata$condition
)

rownames(annotation_col) <- colnames(
  heatmap_data
)

# Draw heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 50 Differentially Expressed Genes"
)

# ============================================================
# Save Heatmap as PNG
# ============================================================

png(
  "TCGA_BRCA_DEG_heatmap.png",
  width = 1800,
  height = 1600,
  res = 200
)

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

#check where R saved it:

getwd()


# ============================================================
#  Alternative: ComplexHeatmap
# ============================================================

ha <- HeatmapAnnotation(
  Condition = annotation_col$Condition
)

Heatmap(
  heatmap_data,
  name = "Z-score",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "TCGA-BRCA Top 50 DEGs"
)


# Save Complex Heatmap as PNG

png(
  "TCGA_BRCA_ComplexHeatmap.png",
  width = 2000,
  height = 1800,
  res = 200
)

Heatmap(
  heatmap_data,
  name = "Z-score",
  top_annotation = HeatmapAnnotation(
    Condition = annotation_col$Condition
  ),
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "TCGA-BRCA Samples",
  row_title = "Top 50 Differentially Expressed Genes"
)

dev.off()


# ============================================================
#  Volcano Plot
# ============================================================

#  Labelling Significance Categories

res_df$Significance <- "Not Significant"

res_df$Significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange > 1
] <- "Upregulated"

res_df$Significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange < -1
] <- "Downregulated"

res_df$minus_log10_padj <- -log10(res_df$padj)


# ============================================================
#  Basic Plot (ggplot2)
# ============================================================

library(ggplot2)

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "TCGA-BRCA Differential Gene Expression",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

# ============================================================
# Save Basic Volcano Plot
# ============================================================

png(
  "TCGA_BRCA_Volcano_Plot.png",
  width = 2000,
  height = 1600,
  res = 200
)

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "TCGA-BRCA Differential Gene Expression",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

dev.off()

# ============================================================
#  Label the Top 10 Genes
# ============================================================

library(ggrepel)

label_genes <- res_df[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  abs(res_df$log2FoldChange) > 1,
]

# Select the 10 most significant
label_genes <- head(
  label_genes[order(label_genes$padj), ],
  10
)


# ============================================================
# Now draw the labelled volcano plot
# ============================================================

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(
      label = gene
    ),
    size = 3,
    max.overlaps = 10
  ) +
  labs(
    title = "TCGA-BRCA Differential Gene Expression",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )


# ============================================================
# Save Labelled Volcano Plot
# ============================================================

png(
  "TCGA_BRCA_Volcano_Plot_Labelled.png",
  width = 2000,
  height = 1600,
  res = 200
)

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 10
  ) +
  labs(
    title = "TCGA-BRCA Differential Gene Expression",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

dev.off()

#Check where it was saved:
getwd()





























