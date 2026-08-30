# Load required libraries

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# Prepare gene IDs by removing Ensembl version numbers

up_gene_ids <- sub("\\..*", "", upregulated$gene)
down_gene_ids <- sub("\\..*", "", downregulated$gene)


# Convert upregulated Ensembl IDs to Entrez IDs

up_entrez <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Check mapped upregulated genes

head(up_entrez)

nrow(up_entrez)

length(up_gene_ids)

(nrow(up_entrez) / length(up_gene_ids)) * 100


# Convert downregulated Ensembl IDs to Entrez IDs

down_entrez <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Check mapped downregulated genes

head(down_entrez)

nrow(down_entrez)

length(down_gene_ids)

(nrow(down_entrez) / length(down_gene_ids)) * 100


# Remove duplicate Entrez IDs

up_entrez <- up_entrez[
  !duplicated(up_entrez$ENTREZID),
]

down_entrez <- down_entrez[
  !duplicated(down_entrez$ENTREZID),
]


# Check unique mapped genes

nrow(up_entrez)

nrow(down_entrez)


# GO Biological Process enrichment for upregulated genes

ego_up_BP <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO BP results for upregulated genes

head(
  as.data.frame(ego_up_BP)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO BP bar plot for upregulated genes

go_up_BP_plot <- barplot(
  ego_up_BP,
  showCategory = 15,
  title = "GO Biological Process - Upregulated Genes"
)

go_up_BP_plot


# Save GO BP plot for upregulated genes

ggsave(
  "GO_BP_Upregulated.png",
  plot = go_up_BP_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)


# GO Biological Process enrichment for downregulated genes

ego_down_BP <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO BP results for downregulated genes

head(
  as.data.frame(ego_down_BP)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO BP bar plot for downregulated genes

go_down_BP_plot <- barplot(
  ego_down_BP,
  showCategory = 15,
  title = "GO Biological Process - Downregulated Genes"
)

go_down_BP_plot


# Save GO BP plot for downregulated genes

ggsave(
  "GO_BP_Downregulated.png",
  plot = go_down_BP_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)
