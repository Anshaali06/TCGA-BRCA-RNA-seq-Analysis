# GO Cellular Component enrichment for upregulated genes

ego_up_CC <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO CC results for upregulated genes

head(
  as.data.frame(ego_up_CC)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO CC bar plot for upregulated genes

go_up_CC_plot <- barplot(
  ego_up_CC,
  showCategory = 15,
  title = "GO Cellular Component - Upregulated Genes"
)

go_up_CC_plot


# Save GO CC plot for upregulated genes

ggsave(
  "GO_CC_Upregulated.png",
  plot = go_up_CC_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)


# GO Cellular Component enrichment for downregulated genes

ego_down_CC <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Check GO CC results for downregulated genes

head(
  as.data.frame(ego_down_CC)[, c(
    "ID",
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )]
)


# Create GO CC bar plot for downregulated genes

go_down_CC_plot <- barplot(
  ego_down_CC,
  showCategory = 15,
  title = "GO Cellular Component - Downregulated Genes"
)

go_down_CC_plot


# Save GO CC plot for downregulated genes

ggsave(
  "GO_CC_Downregulated.png",
  plot = go_down_CC_plot,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300
)
