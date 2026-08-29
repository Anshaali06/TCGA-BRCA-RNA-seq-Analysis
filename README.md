# TCGA-BRCA RNA-seq Differential Gene Expression Analysis

## 📌 Project Overview

This project demonstrates a beginner-friendly RNA-seq differential gene expression analysis using publicly available **TCGA-BRCA (Breast Invasive Carcinoma)** data.

The analysis is performed in **R** using Bioconductor and other commonly used bioinformatics packages. The workflow is documented step-by-step so that beginners can understand not only the code, but also the purpose of each analysis step.

---

## 🎯 Project Objectives

The main objectives of this project are to:

1. Obtain TCGA-BRCA RNA-seq gene expression data.
2. Select **Primary Tumor** and **Solid Tissue Normal** samples.
3. Prepare the gene expression count matrix.
4. Filter low-count genes.
5. Perform differential expression analysis using **DESeq2**.
6. Identify significantly upregulated and downregulated genes.
7. Visualize differential expression using a **Volcano Plot** and **Heatmap**.
8. Perform **Gene Ontology (GO) enrichment analysis** to identify associated biological processes.

---

## 🧬 Dataset

**Project:** TCGA-BRCA  
**Cancer Type:** Breast Invasive Carcinoma  
**Data Type:** RNA-seq gene expression counts  
**Data Source:** The Cancer Genome Atlas (TCGA)

The data were accessed using the `TCGAbiolinks` R package.

For the differential expression analysis, the following sample types were selected:

- Primary Tumor
- Solid Tissue Normal

Metastatic samples were excluded from the tumour-vs-normal comparison.

---

## 🛠️ Tools & Packages

The analysis was performed using **R/RStudio** and the following Bioconductor and R packages:

- **TCGAbiolinks**
- **SummarizedExperiment**
- **DESeq2**
- **dplyr**
- **ggplot2**
- **EnhancedVolcano**
- **pheatmap**
- **ComplexHeatmap**
- **clusterProfiler**
- **org.Hs.eg.db**
- **enrichplot**
- **circlize**
- **grid**

## 🔬 Analysis Workflow

```text
TCGA-BRCA Data
       ↓
Sample Selection
       ↓
Primary Tumor vs Solid Tissue Normal
       ↓
Count Matrix Preparation
       ↓
Low-Count Gene Filtering
       ↓
DESeq2 Differential Expression Analysis
       ↓
Significant DEGs
       ↓
 ┌───────────────┐
 ↓               ↓
Upregulated   Downregulated
 ↓               ↓
 └───────┬───────┘
         ↓
Visualization
 ↓               ↓
Volcano Plot   Heatmap
         ↓
GO Enrichment Analysis
         ↓
Biological Interpretation
```

## 📊 Results

The differential expression analysis was performed to compare gene expression between **Primary Tumor** and **Solid Tissue Normal** samples from TCGA-BRCA.

### Differentially Expressed Genes

Using an adjusted p-value cutoff of **< 0.05** and an absolute log2 fold-change cutoff of **> 1**, significant differentially expressed genes (DEGs) were identified.

The significant DEGs were further classified into:

- **Upregulated genes:** log2FoldChange > 1
- **Downregulated genes:** log2FoldChange < -1

### Volcano Plot

The volcano plot provides a visual representation of the differential expression results.

- The **x-axis** represents log2 fold change.
- The **y-axis** represents adjusted p-value.
- Genes on the right show increased expression in tumour samples.
- Genes on the left show decreased expression in tumour samples.
- Significant genes are highlighted based on the selected statistical and fold-change thresholds.


### Heatmap

The heatmap displays the expression patterns of selected differentially expressed genes across the analysed samples.

It helps visualize similarities and differences in gene expression patterns between **tumour and normal samples**.


### Gene Ontology Enrichment

GO enrichment analysis was performed separately for the upregulated and downregulated genes.

The analysis explores three major Gene Ontology categories:

- **Biological Process (BP)**
- **Molecular Function (MF)**
- **Cellular Component (CC)**

The enriched terms provide insight into the biological functions and processes associated with the identified DEGs.

## 🧾 Conclusion

This project demonstrates a complete workflow for **RNA-seq differential gene expression analysis using TCGA-BRCA data in R**.

The workflow covers data retrieval, sample selection, count matrix preparation, low-count gene filtering, differential expression analysis using DESeq2, identification of upregulated and downregulated genes, visualization, and Gene Ontology enrichment analysis.

Overall, this project provides a practical introduction to how computational approaches can be used to analyse gene expression data and explore its biological significance in cancer research.

---

## 💡 Final Note

Bioinformatics is not something that can be learned by simply reading about it. It becomes easier through **hands-on practice, working with real datasets, understanding each step, and learning from mistakes**.

This guide was created with the aim of making an RNA-seq workflow easier for beginners to follow and understand.

I hope this project helps students take their first step into **RNA-seq analysis, genomics, and bioinformatics**, and encourages them to explore beyond this basic workflow.

> **Learn the code. Understand the biology. Question the results. Keep exploring.**

---

## 👩‍💻 About the Author
### Ansha Ali

**Master’s Graduate | Bioinformatics | Genomics | NGS**

I am a Master's graduate with skills in **Bioinformatics, Genomics, NGS, and computational analysis of biological data**.

This project demonstrates my practical experience in **RNA-seq analysis using R**, covering the workflow from TCGA data acquisition and preprocessing to differential gene expression analysis, visualization, and Gene Ontology enrichment.

I created this documentation to present the workflow in a **clear, structured, and reproducible format**, making it useful for students and beginners who want to learn and practice RNA-seq analysis.







