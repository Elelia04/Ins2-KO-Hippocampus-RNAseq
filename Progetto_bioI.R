#DIFFERENTIAL EXPRESSION ANALYSIS (GSE305902)

#Required libraries 

library(GEOquery)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)   
library(clusterProfiler)
library(ReactomePA)
library(stringr)
library(msigdbr)
library(gplots)
library(dplyr)
library(enrichplot)

#Phenotypical matrix extraction 
geo <- GEOquery::getGEO(GEO = "GSE305902", GSEMatrix = TRUE)
colData <- Biobase::pData(geo[[1]])

#Raw counts extraction
GEOquery::getGEOSuppFiles(GEO = "GSE305902") 



############################ DATA CLEANING #####################################

#Cleaning colData

metadata_clean <- colData[, c("geo_accession", "title")]

metadata_clean$Genotype <- ifelse(grepl("WT", metadata_clean$title), "WT", "KO")
metadata_clean$Genotype <- as.factor(metadata_clean$Genotype)

sample_numbers <- gsub("WT_|KO_", "", metadata_clean$title)

metadata_clean$Snumber <- paste0("Snumber", sample_numbers)
rownames(metadata_clean) <- metadata_clean$Snumber

#Raw counts 
counts_path <- "GSE305902/GSE305902_jjohnson_hippocampus_counts.txt.gz"
counts <- read.delim(counts_path, row.names = 1, header = TRUE)


#Re-ordering counts: colnames(counts) must be in the same identical order of rownames(metadata_clean) for DESeq2
counts <- counts[, rownames(metadata_clean)]

#Decimal rounding 
counts <- round(counts)
all(colnames(counts) == rownames(metadata_clean))



##################################### DESeq ####################################

#Creating DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata_clean,
  design = ~ Genotype
)


#Filtering Low Expressed Genes 
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]


#Run DESeq2 pipeline to obtain dds object containing: size factors, dispersion estimates, linear model fits and hypothesis testing 
dds <- DESeq(dds)


############################# DIAGNOSTIC PLOTS #################################

# PCA 1 (23 samples)

#Data transformation for PCA: making variance constant for all genes
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c("Genotype", "title", "Snumber"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Grabbing sample number from labels
pcaData$JustNumber <- gsub("[^0-9]", "", pcaData$title)

ggplot(pcaData, aes(PC1, PC2, color = Genotype, label = JustNumber)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf,       
    box.padding = 0.3,        
    point.padding = 0.3,      
    segment.color = "grey60", 
    segment.size = 0.3        
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Initial PCA: 23 samples")


#PCA 2: sample 8 removal

dds_no8 <- dds[, colnames(dds) != "Snumber8"]
dds_no8$Genotype <- droplevels(dds_no8$Genotype) 

vsd_no8 <- vst(dds_no8, blind = FALSE)
pcaData_no8 <- plotPCA(vsd_no8, intgroup = c("Genotype", "title", "Snumber"), returnData = TRUE)
percentVar_no8 <- round(100 * attr(pcaData_no8, "percentVar"))

pcaData_no8$JustNumber <- gsub("[^0-9]", "", pcaData_no8$title)

ggplot(pcaData_no8, aes(PC1, PC2, color = Genotype, label = JustNumber)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  xlab(paste0("PC1: ", percentVar_no8[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_no8[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Step 2: Sample 8 Removal")


## PCA 3: removal of all the outliers (4 samples)

outliers <- c("Snumber8", "Snumber6", "Snumber12", "Snumber17")
dds_clean <- dds[, !(colnames(dds) %in% outliers)]
dds_clean$Genotype <- droplevels(dds_clean$Genotype)

vsd_clean <- vst(dds_clean, blind = FALSE)
pcaData_clean <- plotPCA(vsd_clean, intgroup = c("Genotype", "title", "Snumber"), returnData = TRUE)
percentVar_clean <- round(100 * attr(pcaData_clean, "percentVar"))

pcaData_clean$JustNumber <- gsub("[^0-9]", "", pcaData_clean$title)

ggplot(pcaData_clean, aes(PC1, PC2, color = Genotype, label = JustNumber)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  xlab(paste0("PC1: ", percentVar_clean[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_clean[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Step 3: 4 Outliers Removed")



#HEATMAPS 

#23 samples

#Calculating Euclidean distances between samples considering gene expression
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

just_numbers <- gsub("[^0-9]", "", vsd$title)
sample_names <- paste(vsd$Genotype, just_numbers, sep="-")

rownames(sampleDistMatrix) <- sample_names
colnames(sampleDistMatrix) <- sample_names

annotation_col <- data.frame(Genotype = vsd$Genotype)
rownames(annotation_col) <- sample_names

colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         annotation_col = annotation_col,
         main = "Initial Sample-to-Sample Distances (23 samples)")


##HEATMAP WITH CLEANED DATASET

sampleDists_clean <- dist(t(assay(vsd_clean)))
sampleDistMatrix_clean <- as.matrix(sampleDists_clean)

just_numbers_clean <- gsub("[^0-9]", "", vsd_clean$title)
sample_names_clean <- paste(vsd_clean$Genotype, just_numbers_clean, sep="-")

rownames(sampleDistMatrix_clean) <- sample_names_clean
colnames(sampleDistMatrix_clean) <- sample_names_clean

annotation_col_clean <- data.frame(Genotype = vsd_clean$Genotype)
rownames(annotation_col_clean) <- sample_names_clean

#Cleaned heatmap
pheatmap(sampleDistMatrix_clean,
         clustering_distance_rows = sampleDists_clean,
         clustering_distance_cols = sampleDists_clean,
         col = colors, 
         annotation_col = annotation_col_clean,
         main = "Cleaned Sample-to-Sample Distances (19 samples)")


## HEATMAP All samples but no 8

sampleDists_no8 <- dist(t(assay(vsd_no8)))
sampleDistMatrix_no8 <- as.matrix(sampleDists_no8)

just_numbers_no8 <- gsub("[^0-9]", "", vsd_no8$title)
sample_names_no8 <- paste(vsd_no8$Genotype, just_numbers_no8, sep="-")

rownames(sampleDistMatrix_no8) <- sample_names_no8
colnames(sampleDistMatrix_no8) <- sample_names_no8

annotation_col_no8 <- data.frame(Genotype = vsd_no8$Genotype)
rownames(annotation_col_no8) <- sample_names_no8

pheatmap(sampleDistMatrix_no8,
         clustering_distance_rows = sampleDists_no8,
         clustering_distance_cols = sampleDists_no8,
         col = colors, 
         annotation_col = annotation_col_no8,
         main = "Sample-to-Sample Distances (Sample 8 Removed, 22 samples)")



## Confronting Before and After Differentially Expressed Genes

#23 samples
dds <- DESeq(dds)

res_dirty <- results(dds, contrast = c("Genotype", "KO", "WT"), alpha = 0.05)
summary(res_dirty)


#Cleaned up dataset (removal of 4 outliers)
#Filter for low expressed genes (>1)

keep_clean <- rowSums(counts(dds_clean)) > 1
dds_clean <- dds_clean[keep_clean, ]

dds_clean <- DESeq(dds_clean)

res_clean <- results(dds_clean, contrast = c("Genotype", "KO", "WT"), alpha = 0.05)

summary(res_clean)




## Checking Ccnd1 behaviour

res_clean$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(res_clean),
  column = "SYMBOL",
  keytype = "ENSEMBL", 
  multiVals = "first"
)

ccnd1_row <- which(res_clean$symbol == "Ccnd1")

if(length(ccnd1_row) > 0) {
  print(res_clean[ccnd1_row, ])
} else {
  cat("No Ccnd1 gene.\n")
}




# Generate the MA-plot

#23 samples 
plotMA(res_dirty, main = "MA Plot: KO vs WT (Initial Dataset, 23 samples)")

#Cleaned dataset MA plot

plotMA(res_clean, main = "MA Plot: KO vs WT (4 Outliers Removed, 19 samples)")

## Ccnd1 in MA plot

plotMA(res_clean, main = "MA Plot: KO vs WT (Ccnd1)")

ccnd1_id <- "ENSMUSG00000070348"
ccnd1_data <- res_clean[ccnd1_id, ]

points(
  x = ccnd1_data$baseMean, 
  y = ccnd1_data$log2FoldChange, 
  col = "red", 
  lwd = 2, 
  cex = 1.8, 
  pch = 1
)

text(
  x = ccnd1_data$baseMean, 
  y = ccnd1_data$log2FoldChange, 
  labels = "Ccnd1", 
  pos = 4, 
  col = "red", 
  font = 2, 
  cex = 0.9
)


################# HEATMAPS 50 TOP EXPRESSED GENES ##############################

#Running rlog on cleaned dataset
rld_clean <- rlog(dds_clean, blind = FALSE)

#Extracting top 50 genes sorted by their adjusted p-value from clean results
top50_genes <- head(order(res_clean$padj), 50)

just_numbers_clean <- gsub("[^0-9]", "", vsd_clean$title)
sample_names_clean <- paste(vsd_clean$Genotype, just_numbers_clean, sep="-")

annotation_col_clean <- data.frame(Genotype = vsd_clean$Genotype)
rownames(annotation_col_clean) <- sample_names_clean

# VST MATRIX
vsd_clean <- vst(dds_clean, blind = FALSE)

vst_matrix <- assay(vsd_clean)[top50_genes, ]
colnames(vst_matrix) <- sample_names_clean

pheatmap(vst_matrix, 
         cluster_rows = TRUE, #clustering genes with similar gene expression
         show_rownames = FALSE,
         cluster_cols = TRUE, 
         annotation_col = annotation_col_clean,
         scale = "row", #applies Z-score on rows
         main = "Top 50 DE Genes, VST (Clean Dataset, Row Scaled)")

# RLOG MATRIX
rlog_matrix <- assay(rld_clean)[top50_genes, ]
colnames(rlog_matrix) <- sample_names_clean

pheatmap(rlog_matrix, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         cluster_cols = TRUE, 
         annotation_col = annotation_col_clean,
         scale = "row",
         main = "Top 50 DE Genes, RLOG (Clean Dataset, Row Scaled)")



## TOP 50 GENES FOR PPT

top50_data <- res_clean[top50_genes, ]

top50_df <- as.data.frame(top50_data)
top50_summary <- top50_df[, c("baseMean", "log2FoldChange", "padj")]

top50_summary$baseMean <- round(top50_summary$baseMean, 1)
top50_summary$log2FoldChange <- round(top50_summary$log2FoldChange, 2)

top50_summary$padj <- formatC(top50_summary$padj, format = "e", digits = 2)

head(top50_summary, 30)



## VOLCANO PLOT

#Isolating 19 significant genes
degs_19 <- as.data.frame(res_clean[which(res_clean$padj < 0.05), ])

#Adding gene symbol
degs_19$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = rownames(degs_19), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
degs_19$Gene_Symbol <- ifelse(is.na(degs_19$Gene_Symbol), rownames(degs_19), degs_19$Gene_Symbol)

degs_19$Gene_Symbol[rownames(degs_19) == "ENSMUSG00000096999"] <- "Gm26793"

#-log10 of padj for vertical scaling
degs_19$log_sig <- -log10(degs_19$padj)

ggplot(degs_19, aes(x = log2FoldChange, y = log_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = log2FoldChange > 0), size = 5, alpha = 0.8) +
  geom_text(aes(label = Gene_Symbol), 
            vjust = -1, 
            fontface = "italic", 
            size = 3.5, 
            color = "black", 
            check_overlap = TRUE) +
  scale_color_manual(values = c("TRUE" = "#d95f02", "FALSE" = "#7570b3"), guide = "none") +
  labs(
    title = "Statistical Confidence of the 19 DEGs",
    x = "Log2 Fold Change (KO vs WT)",
    y = "-Log10 Adjusted P-value"
  ) +
  ylim(0, 50) + 
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )


# UP REGULATED AND DOWN REGULATED DEGs

#Upregulated 
sum(degs_19$log2FoldChange > 0)

#Downregulated
sum(degs_19$log2FoldChange < 0)



########################## GO Visualization ####################################

res_df <- as.data.frame(res_clean)

#Mapping IDs
res_df$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(res_df),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = rownames(res_df),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

#Definition of Universe and Gene Lists
gene_universe <- res_df$symbol[!is.na(res_df$padj)] 
entrez_universe <- res_df$entrez[!is.na(res_df$padj)]


#Filtring table to extract significative genes (padj < 0.05) and separating them into up and down genes
up_genes_strict <- res_df$symbol[res_df$log2FoldChange > 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]
down_genes_strict <- res_df$symbol[res_df$log2FoldChange < 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]

#Compacting thw two lists into one
gene_list_strict <- list(UP = na.omit(up_genes_strict), DOWN = na.omit(down_genes_strict))


#Phase I, BH

#Biological Process (BP)
go_bp <- compareCluster(
  geneClusters = gene_list_strict, 
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db, #Database for Mus Musculus
  keyType      = "SYMBOL", #Specifies that Gene IDs are official
  ont          = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
  universe     = gene_universe
)


if(!is.null(go_bp)) {
  print(
    dotplot(go_bp, showCategory = 10, title = "GO Enrichment: BP (Strict BH Correction)") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) + 
      theme(axis.text.y = element_text(size = 8, face = "plain")) 
  )
}

#Molecular Function (MF), Bh
go_mf <- compareCluster(
  geneClusters = gene_list_strict,
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe     = gene_universe
)

if(!is.null(go_mf)) {
  print(
    dotplot(go_mf, showCategory = 10, title = "GO Enrichment: MF (Strict BH Correction)") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
      theme(axis.text.y = element_text(size = 8))
  )
}


#Phase 2: explorative approach (hypotesis)

#Biological Process (BP)
go_bp_none <- compareCluster(
  geneClusters = gene_list_strict, 
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  universe     = gene_universe
)

print(
  dotplot(go_bp_none, showCategory = 10, title = "GO BP (Exploratory - Nominal p-value)") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
    theme(axis.text.y = element_text(size = 8))
)

#Molecular Function (MF)
go_mf_none <- compareCluster(
  geneClusters = gene_list_strict, 
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  universe     = gene_universe
)

print(
  dotplot(go_mf_none, showCategory = 10, title = "GO MF (Exploratory, Nominal p-value)") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
    theme(axis.text.y = element_text(size = 8))
)



#######################  PATHWAY ENRICHMENT ANALYSIS WITH REACTOME #############


#Entrez IDs for the 19 significant core DEGs (padj < 0.05)
up_entrez_strict <- res_df$entrez[res_df$log2FoldChange > 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]
down_entrez_strict <- res_df$entrez[res_df$log2FoldChange < 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]

#Strict list for comparative analysis
gene_list_reactome_strict <- list(
  UP = na.omit(up_entrez_strict), 
  DOWN = na.omit(down_entrez_strict)
)

#Performing the comparative pathway enrichment analysis with Reactome
reactome_strict <- compareCluster(
  geneClusters  = gene_list_reactome_strict,
  fun           = "enrichPathway",
  organism      = "mouse",     
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",        
  universe      = entrez_universe
)

if (!is.null(reactome_strict) && nrow(as.data.frame(reactome_strict)) > 0) {
  print(
    dotplot(reactome_strict, showCategory = 10, title = "Reactome Pathway Enrichment (Strict BH)") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
      theme(axis.text.y = element_text(size = 8))
  )
} else {
  cat("No Reactome pathways reached statistical significance under strict BH correction.\n")
}


## CHECKING LOST SIGNIFICAT GENES WHEN MAPPING 
#How many DEGs lost while converting

length(up_entrez_strict); length(na.omit(up_entrez_strict))
length(down_entrez_strict); length(na.omit(down_entrez_strict))
length(entrez_universe); length(na.omit(entrez_universe))

sum(is.na(res_df$entrez[res_df$log2FoldChange > 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]))
sum(is.na(res_df$entrez[res_df$log2FoldChange < 0 & res_df$padj < 0.05 & !is.na(res_df$padj)]))


reactome_up_raw <- enrichPathway(gene = na.omit(up_entrez_strict),
                                 organism = "mouse",
                                 universe = entrez_universe,
                                 pvalueCutoff = 1,   
                                 pAdjustMethod = "BH")

if (!is.null(reactome_up_raw)) { head(as.data.frame(reactome_up_raw)) }

reactome_down_raw <- enrichPathway(gene = na.omit(down_entrez_strict),
                                   organism = "mouse",
                                   universe = entrez_universe,
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH")

if (!is.null(reactome_down_raw)) { head(as.data.frame(reactome_down_raw)) }


# REACTOME EXPLORATIVE 

reactome_explorative <- compareCluster(
  geneClusters  = gene_list_reactome_strict,
  fun           = "enrichPathway",
  organism      = "mouse",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "none",   
  universe      = entrez_universe
)

if (!is.null(reactome_explorative) && nrow(as.data.frame(reactome_explorative)) > 0) {
  print(
    dotplot(reactome_explorative, showCategory = 10, 
            title = "Reactome Pathway Enrichment (Exploratory)") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
      theme(axis.text.y = element_text(size = 8))
  )
} else {
  cat("No Reactome pathway with nominal p-value\n")
}




######################### GSEA #################################################

#Hallmark GSEA
h_gene_sets <- msigdbr(species = "mouse", category = "H")

#Table with 2 columns with pathway name and Ensembl codes
msigdbr_t2g <- as.data.frame(dplyr::distinct(h_gene_sets, gs_name, ensembl_gene))

#Ranked list of all genes
res_gsea_df <- as.data.frame(res_clean)
res_gsea_df <- res_gsea_df[!is.na(res_gsea_df$stat), ]
rank <- res_gsea_df$stat
names(rank) <- rownames(res_gsea_df)
rank <- sort(rank, decreasing = TRUE) #Ranking from most upregulated to less 

gsea_res <- GSEA(
  geneList      = rank,
  TERM2GENE     = msigdbr_t2g,
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_res_df <- as.data.frame(gsea_res)
gsea_res_df$Description <- gsub("HALLMARK_", "", gsea_res_df$Description)
gsea_res_df$Description <- gsub("_", " ", gsea_res_df$Description)
gsea_res_df$Regulation <- ifelse(gsea_res_df$NES > 0, "Up-regulated (Activated)", "Down-regulated (Suppressed)")
gsea_res_df <- gsea_res_df[order(gsea_res_df$NES, decreasing = TRUE), ]

# Bar plot NES
ggplot(gsea_res_df, aes(x = reorder(Description, NES), y = NES, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Up-regulated (Activated)" = "#d95f02", "Down-regulated (Suppressed)" = "#1f78b4")) +
  labs(title = "Hallmark Gene Sets Ranked by NES", x = "Hallmark Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")

#Enrichment plots, top 5
gsea_res_ordered <- gsea_res_df[order(gsea_res_df$pvalue), ]
top5_plots <- head(gsea_res_ordered$ID, 5)

for (pathway_id in top5_plots) {
  clean_title <- gsub("_", " ", gsub("HALLMARK_", "", pathway_id))
  print(gseaplot2(gsea_res, geneSetID = pathway_id, 
                  title = paste("Enrichment Plot:", clean_title), 
                  pvalue_table = TRUE))
}

#Heatmap
top_pathway_name <- gsea_res_ordered$Description[1]
core_ensembl <- unlist(strsplit(gsea_res_ordered$core_enrichment[1], "/"))

counts_matrix <- assay(vsd_clean)
heatmap_mat <- counts_matrix[rownames(counts_matrix) %in% core_ensembl, , drop = FALSE]

gene_symbols <- mapIds(org.Mm.eg.db, keys = rownames(heatmap_mat), 
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols[is.na(gene_symbols)] <- rownames(heatmap_mat)[is.na(gene_symbols)]
rownames(heatmap_mat) <- make.unique(gene_symbols)

annotation_col <- data.frame(
  Genotype = colData(dds_clean)[colnames(heatmap_mat), "Genotype"]
)
rownames(annotation_col) <- colnames(heatmap_mat)

pheatmap(
  heatmap_mat,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = TRUE,
  show_rownames = FALSE,
  fontsize_row = 6,
  annotation_col = annotation_col,
  main = paste("Core Enrichment:", top_pathway_name),
  color = colorRampPalette(c("#1f78b4", "white", "#e31a1c"))(50)
)



############################### QUESTIONS ######################################

res_clean_df <- as.data.frame(res_clean)

#Question 1
genes_padj_005 <- sum(res_clean_df$padj < 0.05, na.rm = TRUE)
genes_padj_01  <- sum(res_clean_df$padj < 0.10, na.rm = TRUE)

#Question 2
genes_up_005   <- sum(res_clean_df$padj < 0.05 & res_clean_df$log2FoldChange > 0, na.rm = TRUE)
genes_down_005 <- sum(res_clean_df$padj < 0.05 & res_clean_df$log2FoldChange < 0, na.rm = TRUE)

#Question 3

go_df <- as.data.frame(go_bp)

num_categorie <- sum(go_df$p.adjust < 0.05, na.rm = TRUE)

cat("Significative categories:", num_categorie, "\n")


#Question 4
geni_top <- go_df$Count[1]
nome_top <- go_df$Description[1]

cat("Category:", nome_top, "| Number of genes:", geni_top, "\n")

######################## Csv format for the table ##############################

write.csv(as.data.frame(res_clean), 
          file = "All_Genes_DE_Results.csv", 
          row.names = TRUE)

# 19 Significant genes (DEGs with padj < 0.05)
degs_signif <- res_clean[which(res_clean$padj < 0.05), ]

# 19 Significant Genes Table (TXT)
write.table(as.data.frame(degs_signif), 
            file = "Significant_DEGs_Results.txt", 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA)
