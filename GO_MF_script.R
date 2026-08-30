# GO Molecular Function enrichment for upregulated genes

ego_up_MF <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO MF results for upregulated genes

head(
  as.data.frame(ego_up_MF)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO MF bar plot for upregulated genes

go_up_MF_plot <- barplot(
  ego_up_MF,
  showCategory = 15,
  title = "GO Molecular Function - Upregulated Genes"
)

go_up_MF_plot


# Save GO MF plot for upregulated genes

ggsave(
  "GO_MF_Upregulated.png",
  plot = go_up_MF_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)


# GO Molecular Function enrichment for downregulated genes

ego_down_MF <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO MF results for downregulated genes

head(
  as.data.frame(ego_down_MF)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO MF bar plot for downregulated genes

go_down_MF_plot <- barplot(
  ego_down_MF,
  showCategory = 15,
  title = "GO Molecular Function - Downregulated Genes"
)

go_down_MF_plot


# Save GO MF plot for downregulated genes

ggsave(
  "GO_MF_Downregulated.png",
  plot = go_down_MF_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)
