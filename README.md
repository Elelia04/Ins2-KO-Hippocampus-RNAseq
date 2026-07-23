# Ins2-KO-Hippocampus-RNAseq

Differential expression analysis of hippocampal RNA-seq data from female *Ins2⁻/⁻* knockout mice (GEO accession [GSE305902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305902)), replicating and extending the transcriptomic analysis from:

> Baehring SK, O'Leary TP, Cen HH, et al. (2026). Loss of brain insulin production impairs learning and memory in female mice. *Metabologia* 2:4. https://doi.org/10.1007/s44357-026-00004-8

📊 **[View the full rendered analysis (HTML report)](https://elelia04.github.io/Ins2-KO-Hippocampus-RNAseq/)**

## Dataset

- **Organism**: *Mus musculus*
- **Tissue**: hippocampus
- **Design**: 23 samples, *Ins2⁺/⁺* (WT) vs *Ins2⁻/⁻* (KO), female mice, 12 months of age

## Pipeline overview

The analysis (`Progetto_bioI.R`) follows these steps:

1. **Data retrieval** — metadata and raw count matrix downloaded from GEO via `GEOquery`
2. **Metadata cleaning** — genotype extraction from sample titles, standardized `Snumber` sample IDs
3. **DESeq2 setup** — `DESeqDataSetFromMatrix`, low-count gene filtering (`rowSums(counts) > 1`), `design = ~ Genotype`
4. **Outlier detection (iterative)**:
   - PCA on all 23 samples reveals one extreme outlier (Sample 8)
   - PCA after removing Sample 8 reveals two further outliers (Samples 6, 12, 17), consistent with tissue-contamination signal (Pmch, a hypothalamic marker gene, among top PCA loadings)
   - Final cleaned dataset: 19 samples, 4 outliers removed — matching the outlier set reported in the original paper's methods
5. **Differential expression** — DESeq2 Wald test (`Genotype: KO vs WT`) on both the full (23 samples) and cleaned (19 samples) datasets, compared side-by-side
6. **Ccnd1 validation** — direct check of *Ccnd1* (the top downregulated gene reported in the original paper) to validate the pipeline against published results
7. **Visualization**:
   - Sample-to-sample distance heatmaps (23 samples / 22 samples / 19 samples cleaned)
   - MA-plots (dirty vs clean dataset, with *Ccnd1* highlighted)
   - Top 50 DE genes heatmaps (VST and rlog, row-scaled)
   - Volcano plot of the 19 significant DEGs (padj < 0.05)
8. **Functional enrichment**:
   - GO enrichment (Biological Process, Molecular Function) via `clusterProfiler::enrichGO`, both with strict BH correction and as an exploratory nominal p-value analysis (given the small number of DEGs)
   - Reactome pathway enrichment via `ReactomePA::enrichPathway` (strict BH and exploratory)
   - GSEA with Hallmark gene sets (`msigdbr` + `clusterProfiler::GSEA`), including NES barplot, top 5 enrichment plots, and a heatmap of the top pathway's core enrichment genes
9. **Summary statistics** — DEG counts at padj < 0.05 / < 0.10, up/down split, significant GO categories, and export of full and significant DEG tables (CSV/TXT)

## Key findings

- **4 outlier samples** identified independently via PCA/hierarchical clustering, matching those excluded in the original study
- **19 significant DEGs** (padj < 0.05) in the cleaned dataset
- ***Ccnd1*** confirmed as the top downregulated gene (log2FC ≈ ‑1.16, padj ≈ 8×10⁻⁴⁴), consistent with the original paper
- **GSEA (Hallmark)**: oxidative phosphorylation is the top upregulated pathway in *Ins2⁻/⁻* mice, consistent with the GO-based GSEA results (aerobic respiration, oxidative phosphorylation) reported in the original study
- GO/Reactome enrichment under strict BH correction yields limited/no significant categories, reflecting the small number of DEGs — addressed with a complementary exploratory (nominal p-value) analysis

## Requirements

```r
install.packages(c("ggplot2", "ggrepel", "pheatmap", "stringr", "dplyr", "gplots", "msigdbr"))
BiocManager::install(c("GEOquery", "DESeq2", "AnnotationDbi", "org.Mm.eg.db",
                        "clusterProfiler", "ReactomePA", "enrichplot"))
```

## Outputs

- `All_Genes_DE_Results.csv` — full DESeq2 results table (cleaned dataset)
- `Significant_DEGs_Results.txt` — 19 significant DEGs (padj < 0.05)

## Reference

Baehring SK, O'Leary TP, Cen HH, et al. (2026). Loss of brain insulin production impairs learning and memory in female mice. *Metabologia* 2:4. https://doi.org/10.1007/s44357-026-00004-8
