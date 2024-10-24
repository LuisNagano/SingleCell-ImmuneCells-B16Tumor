# Single-Cell RNA-Seq Analysis of Infiltrating Immune Cells in B16 Tumor

## Table of Contents

- [Overview](#overview)
- [Package Installation](#package-installation)
- [Creating the Seurat Object](#creating-the-seurat-object)
- [Cell Filtering](#cell-filtering)
- [Data Normalization and Clustering](#data-normalization-and-clustering)
- [Cell Type Identification](#cell-type-identification)
- [Subsetting and Reclustering](#subsetting-and-reclustering)
- [Seurat Native Visualizations](#seurat-native-visualizations)
  - [1. FeaturePlot](#1-featureplot)
  - [2. VlnPlot](#2-vlnplot)
  - [3. RidgePlot](#3-ridgeplot)
  - [4. DotPlot](#4-dotplot)
  - [5. Heatmap](#5-heatmap)
  - [6. FeatureScatter](#6-featurescatter)
- [Custom Visualizations](#custom-visualizations)
  - [1. Cell Type Proportion Plot](#1-cell-type-proportion-plot)
  - [2. UMI Count Plot](#2-umi-count-plot)
  - [3. Gene Signature Plot](#3-gene-signature-plot)
  - [4. Cell Identification by Barcode](#4-cell-identification-by-barcode)
  - [5. Percentage of Positive Cells by Cluster](#5-percentage-of-positive-cells-by-cluster)
  - [6. Median Gene Splitting](#6-median-gene-splitting)
- [Group-wise Analysis](#group-wise-analysis)
  - [1. Group-wise Cell Proportion Plot](#1-group-wise-cell-proportion-plot)
  - [2. Heatmap of Sample and Cell Type Similarity](#2-heatmap-of-sample-and-cell-type-similarity)
  - [3. Median Splitting Across Groups](#3-median-splitting-across-groups)
  - [4. Differential Gene Expression Analysis Between Groups](#4-differential-gene-expression-analysis-between-groups)
  - [5. Gene Set Enrichment Analysis (GSEA) Across Groups](#5-gene-set-enrichment-analysis-gsea-across-groups)
- [Conclusions](#conclusions)

## Overview

This document provides a step-by-step analysis of a publicly available single-cell RNA-Seq dataset as described by Ishizuka et al. (2019). We utilized the R package `Seurat` to process and analyze the data, which were generated using the 10X Genomics platform from leukocytes infiltrating B16 mouse tumors. The article can be accessed at: [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30559380).

The goal of this analysis is to offer a practical guide for researchers with basic R knowledge, facilitating the understanding and analysis of single-cell RNA-Seq experiments. We also include custom visualizations for a deeper exploration of the data.

A rendered version of this R notebook is available at: [Single-Cell RNA-Seq Analysis](https://erilu.github.io/single-cell-rnaseq-analysis/).

---

## Package Installation

```r
library(Seurat)
library(tidyverse)
library(future)
```

This analysis was conducted using `Seurat` version 3.0.2. The `future` package is employed for parallel processing, speeding up time-consuming steps. `tidyverse` provides efficient tools for data manipulation and visualization.

## Creating the Seurat Object

Data can be obtained from: [GEO - GSE110746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110746). Download the `GSE110746_barcodes.tsv.gz`, `GSE110746_genes.tsv.gz`, and `GSE110746_matrix.mtx.gz` files. These files contain the cell barcode matrices, gene names, and expression values, respectively, generated by 10X Genomics' Cell Ranger software.

```r
# Load raw data
raw_data <- Read10X(data.dir = "raw_data/haining")
```

After reading the data, we create a Seurat object:

```r
b16 <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = 200, project = "b16", names.delim = "-", names.field = 2)
```

- `min.cells`: Filters genes expressed in at least 3 cells.
- `min.features`: Filters cells expressing at least 200 genes.

This initial filtering helps remove rare genes and low-quality cells.

Check the object's metadata:

```r
head(b16@meta.data)
```

```r
unique(b16@meta.data$orig.ident)
```

```r
b16
```

## Cell Filtering

We add mitochondrial and ribosomal gene information to the metadata and visualize their distribution.

```r
# Identify mitochondrial and ribosomal genes
grep("^mt-", rownames(b16@assays$RNA@data), value = TRUE)
grep("Rps|Rpl|Mrpl|Mrps", rownames(b16@assays$RNA@data), value = TRUE)
```

```r
# Calculate percentages
b16[["percent.mito"]] <- PercentageFeatureSet(object = b16, pattern = "^mt-")
b16[["percent.ribo"]] <- PercentageFeatureSet(object = b16, pattern = "Rps|Rpl|Mrpl|Mrps")

# Visualize percentages
plot1 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "percent.mito") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot2 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot3 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "percent.ribo") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
```

![Mitochondrial and Ribosomal Percentage Plots](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/percent-mito-ribo-nfeature-graphs.jpeg)

We identify cells with high percentages of mitochondrial and ribosomal genes, possibly indicating compromised cells.

```r
b16 <- subset(x = b16, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40)
```

Check the reduction in the number of cells:

```r
b16
```

## Data Normalization and Clustering

Proceed with normalization, variable gene identification, and data scaling.

```r
b16 <- NormalizeData(b16, normalization.method = "LogNormalize", scale.factor = 10000)
b16 <- FindVariableFeatures(b16, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(b16)
plan("multiprocess", workers = 4)
b16 <- ScaleData(b16, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multiprocess", workers = 1)
```

Perform Principal Component Analysis (PCA) for dimensionality reduction:

```r
b16 <- RunPCA(b16, features = VariableFeatures(object = b16), ndims.print = 1:2)
```

Evaluate the number of PCs to retain using the elbow plot:

```r
ElbowPlot(object = b16, ndims = 50)
```

![Elbow Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/scree-plot.jpeg)

We opt for 40 PCs for subsequent analyses and proceed with clustering and t-SNE visualization:

```r
b16 <- FindNeighbors(b16, dims = 1:40)
b16 <- FindClusters(b16, resolution = c(0.2, 0.4, 0.6))
b16 <- RunTSNE(b16, dims = 1:40)

# Save processed object
# save(b16, file = "b16_all_v3.Robj")
```

We explore different resolutions to determine the granularity of clusters:

```r
r02 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.2") + labs(title = "Resolution 0.2")
r04 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.4") + labs(title = "Resolution 0.4")
r06 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.6") + labs(title = "Resolution 0.6")

CombinePlots(plots = list(r02, r04, r06), ncol = 3)
```

![Resolution Comparison](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/compare-resolutions.jpeg)

We choose a resolution of 0.6 to balance specificity and biological relevance.

## Cell Type Identification

We identify specific marker genes for each cluster:

```r
# Identify markers
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Save markers
# save(b16.markers, file = "b16_all_markers_res06_v3.Robj")
```

Visualize the top markers:

```r
top_markers_all <- b16.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top_markers_all
```

We assign cell types based on the identified markers:

```r
# Define new cluster identities
new.cluster.ids <- c(
  "0- TAM 1",
  "1- Monocyte",
  "2- CD8 T cell",
  "3- TAM 2",
  "4- Melanocyte",
  "5- TAM 3",
  "6- cDC2",
  "7- Melanocyte",
  "8- cDC1",
  "9- Migratory DC",
  "10- NK",
  "11- TAM 4",
  "12- Treg",
  "13- Melanocyte",
  "14- Fibroblast",
  "15- pDC",
  "16- TAM 5"
)

b16.labeled <- b16
names(new.cluster.ids) <- levels(b16.labeled)
b16.labeled <- RenameIdents(b16.labeled, new.cluster.ids)

DimPlot(b16.labeled, label = TRUE, reduction = "tsne")
```

![Labeled t-SNE Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/labeled-tsne-plot.jpeg)

## Subsetting and Reclustering

We focus on immune cells by excluding non-immune cells such as melanocytes and fibroblasts, verifying the expression of CD45 (Ptprc):

```r
FeaturePlot(b16.labeled, "Ptprc")
```

![FeaturePlot of Ptprc](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/ptprc-featureplot.jpeg)

```r
VlnPlot(b16.labeled, "Ptprc", pt.size = 0.5) + NoLegend()
```

![Violin Plot of Ptprc](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/ptprc-vlnplot.jpeg)

We exclude clusters with low CD45 expression:

```r
# Remove non-immune clusters (4, 7, 13, 14)
b16 <- subset(b16, idents = c(0:16)[-c(5,8,14,15)])
DimPlot(b16, label = TRUE, reduction = "tsne")
```

![Non-Immune Cells Removed](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/remove-non-immune-tsne.jpeg)

We reprocess the filtered data:

```r
# Normalize, scale, and recluster
b16 <- NormalizeData(b16, normalization.method = "LogNormalize", scale.factor = 10000)
b16 <- FindVariableFeatures(b16, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(b16)
plan("multiprocess", workers = 4)
b16 <- ScaleData(b16, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multiprocess", workers = 1)
b16 <- RunPCA(b16, features = VariableFeatures(object = b16))
b16 <- FindNeighbors(b16, dims = 1:40)
b16 <- FindClusters(b16, resolution = 0.6)
b16 <- RunTSNE(b16, dims = 1:40)

# Save object
# save(b16, file = "b16_immune_v3.Robj")
```

Visualize the clusters again:

```r
DimPlot(b16, label = TRUE, reduction = "tsne")
```

![Reclustered t-SNE Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/tsne-reclustered-b16-immune.jpeg)

We re-identify markers for the immune subset:

```r
# Identify markers in the immune subset
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Save markers
# save(b16.markers, file = "b16_immune_markers_res06_v3.Robj")
```

Assign cell types again:

```r
new.cluster.ids <- c(
  "0- Classical Monocyte",
  "1- TAM 1",
  "2- TAM 2",
  "3- CD8 T cell",
  "4- TAM 3",
  "5- cDC2",
  "6- Patrolling Monocyte",
  "7- Treg",
  "8- cDC1",
  "9- NK cell",
  "10- TAM 4",
  "11- Migratory DC",
  "12- pDC",
  "13- TAM 5"
)

b16.labeled <- b16
names(new.cluster.ids) <- levels(b16.labeled)
b16.labeled <- RenameIdents(b16.labeled, new.cluster.ids)

DimPlot(b16.labeled, label = TRUE, reduction = "tsne")
```

![Labeled Immune Subset t-SNE Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/labeled-tsne-b16-immune.jpeg)

Reorder clusters for better visualization:

```r
# Reorder clusters
reordered.new.cluster.ids <- c(
  "0- Classical Monocyte",
  "6- Patrolling Monocyte",
  "1- TAM 1",
  "2- TAM 2",
  "4- TAM 3",
  "10- TAM 4",
  "13- TAM 5",
  "8- cDC1",
  "5- cDC2",
  "11- Migratory DC",
  "12- pDC",
  "3- CD8 T cell",
  "7- Treg",
  "9- NK cell"
)

b16.labeled@active.ident <- factor(b16.labeled@active.ident, levels = reordered.new.cluster.ids)

DimPlot(b16.labeled, label = TRUE, reduction = "tsne")
```

![Reordered Clusters](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/reordered-tsne-b16-immune.jpeg)

We aggregate clusters for an overview:

```r
combined.cluster.ids <- c(
  "Monocyte",
  "TAM",
  "TAM",
  "CD8",
  "TAM",
  "cDC2",
  "Monocyte",
  "Treg",
  "cDC1",
  "NK",
  "TAM",
  "Migratory DC",
  "pDC",
  "TAM"
)

b16.combined <- b16
names(combined.cluster.ids) <- levels(b16.combined)
b16.combined <- RenameIdents(b16.combined, combined.cluster.ids)

# Reorder for myeloid and lymphoid hierarchy
b16.combined@active.ident <- factor(b16.combined@active.ident, levels = c("TAM", "Monocyte", "cDC1", "cDC2", "pDC", "Migratory DC", "CD8", "Treg", "NK"))
# Save aggregated identities
b16.combined[["res06_aggregated"]] <- b16.combined@active.ident
```

Visualize the aggregated clusters:

```r
DimPlot(b16.combined, label = TRUE, reduction = "tsne")
```

![Aggregated t-SNE Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/aggregated-tsne-b16.jpeg)

We identify markers for the aggregated clusters:

```r
# Identify aggregated markers
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Save markers
# save(b16.markers, file = "b16_aggregate_markers_res06_v3.Robj")
```

```r
top_markers_aggregate <- b16.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top_markers_aggregate
```

## Seurat Native Visualizations

### 1. FeaturePlot

`FeaturePlot()` allows visualization of specific gene expression within the t-SNE reduced space.

```r
fp <- FeaturePlot(b16.combined, "Ccr2", pt.size = 0.6)
tp <- DimPlot(b16.combined, pt.size = 0.6, label = TRUE, reduction = "tsne")
CombinePlots(plots = list(tp, fp))
```

![Ccr2 FeaturePlot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/featureplot-ccr2.jpeg)

### 2. VlnPlot

`VlnPlot()` generates violin plots to quantify gene expression per cluster.

```r
VlnPlot(b16.combined, "Ccr2", pt.size = 0.1)
```

![Violin Plot of Ccr2](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/vlnplot-ccr2.jpeg)

### 3. RidgePlot

`RidgePlot()` offers an alternative density-based visualization of gene expression.

```r
RidgePlot(b16.combined, "Ccr2")
```

![RidgePlot of Ccr2](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/ridgeplot-ccr2.jpeg)

### 4. DotPlot

`DotPlot()` facilitates the comparison of multiple genes across clusters.

```r
cluster_specific_genes <- c("C1qa", "Apoe", "Ly6c2", "Ccr2", "Xcr1", "Clec9a", "Cd209a", "Clec4a4", "Siglech", "Ccr9", "Ccr7", "Fscn1", "Cd8a", "Cd3e", "Foxp3", "Il2ra", "Ncr1", "Klra4")

DotPlot(b16.combined, features = rev(cluster_specific_genes), cols = "RdBu", plot.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

![DotPlot of Specific Genes](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/dotplot-cell-specific-markers.jpeg)

### 5. Heatmap

`DoHeatmap()` visualizes gene expression across individual cells within clusters.

```r
top_b16_markers <- b16.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DoHeatmap(b16.combined, cells = WhichCells(b16.combined, downsample = 50, seed = 1), features = top_b16_markers$gene, size = 3.5, angle = 15, hjust = 0)
```

![Heatmap of Specific Markers](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/heatmap-cell-specific-markers.jpeg)

### 6. FeatureScatter

`FeatureScatter()` assesses the correlation between two genes within a cluster.

```r
FeatureScatter(cd8, "Tox", "Gzmb") +
  geom_smooth(method = "lm", se = FALSE)
```

![FeatureScatter of Tox vs Gzmb](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/featurescatter-tox-gzmb.jpeg)

## Custom Visualizations

### 1. Cell Type Proportion Plot

Visualize the percentage of cells each cluster represents in the total.

```r
plot_proportions <- function(seuratobj){
  proporcao <- table(seuratobj@active.ident) / length(seuratobj@active.ident) * 100
  proporcao <- as.data.frame(proporcao)
  ggplot(proporcao, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", aes(fill = Var1)) +
    theme_classic() +
    labs(x = NULL, y = "Total Cell Percentage", title = "Cell Type Distribution") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    )
}

plot_proportions(b16.combined)
```

![Cell Type Proportion](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/plot-celltype-proportions.jpeg)

### 2. UMI Count Plot

Assess the distribution of UMIs per cell.

```r
umi_plot <- function(seuratobj) {
  seuratobj[["tsne1"]] <- Embeddings(object = seuratobj, reduction = "tsne")[, "tSNE_1"]
  seuratobj[["tsne2"]] <- Embeddings(object = seuratobj, reduction = "tsne")[, "tSNE_2"]
  ggplot(seuratobj@meta.data, aes(x = tsne1, y = tsne2)) +
    geom_point(aes(color = log(nCount_RNA))) +
    scale_color_gradientn(colors = viridis::magma(10))
}

umi_plot(b16.combined)
```

![UMI Count Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/umi-ncountrna-plot.jpeg)

Summary of statistics:

```r
summarize_metadata <- function(seuratobject, group.by = "RNA_snn_res.0.6") {
  seuratobject@meta.data %>% 
    group_by_at(group.by) %>% 
    summarize(
      avg_UMI = mean(nCount_RNA),
      avg_nFeature = mean(nFeature_RNA),
      avg_percent.mito = mean(percent.mito),
      avg_percent.ribo = mean(percent.ribo)
    ) %>%
    arrange(desc(avg_UMI))
}

summarize_metadata(b16.combined, group.by = "RNA_snn_res.0.6")
```

### 3. Gene Signature Plot

Overlay gene signature scores on the t-SNE plot.

```r
signature_genes <- c("Havcr2", "Tigit", "Lag3","Pdcd1","Tox")

signature_plot <- function(seuratobj, signature_genes){
  sig_data <- as.data.frame(Embeddings(object = seuratobj, reduction = "tsne"))
  sig_data$signature <- rowMeans(t(seuratobj@assays$RNA@scale.data[signature_genes, ]))
  ggplot(sig_data, aes(tSNE_1, tSNE_2)) +
    geom_point(aes(color = signature)) +
    scale_color_gradientn(colors = topo.colors(8)) +
    theme_classic()
}

signature_plot(b16.combined, signature_genes)
```

![Gene Signature Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/signature-plot.jpeg)

### 4. Cell Identification by Barcode

Isolate and visualize specific cell subpopulations.

```r
# Extract CD8 cells
cd8 <- subset(b16.combined, idents = "CD8")

# Reclustering CD8
cd8 <- FindNeighbors(cd8, dims = 1:40)
cd8 <- FindClusters(cd8, resolution = 0.2)

DimPlot(cd8, label = TRUE, reduction = "tsne")
```

![CD8 Subclustered t-SNE Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/tsne-plot-cd8-subclustered.jpeg)

Extract barcodes of a specific subcluster:

```r
cd8_cluster1_barcodes <- WhichCells(cd8, idents = 1)
```

Highlight cells on the t-SNE:

```r
ta <- DimPlot(b16, label = TRUE, cells.highlight = cd8_cluster1_barcodes, reduction = "tsne") + labs(title = "Using cells.highlight")

b16@meta.data$cd8_cluster1 <- ifelse(rownames(b16@meta.data) %in% cd8_cluster1_barcodes, TRUE, FALSE)
tb <- DimPlot(b16, group.by = "cd8_cluster1", reduction = "tsne") + labs(title = "Adding to meta.data and using group.by")

CombinePlots(plots = list(ta, tb))
```

![Highlighted Cells on t-SNE](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/highlight-cells-barcode-tsne-plot.jpeg)

Assign a new identity to isolated cells:

```r
levels(b16.combined@active.ident) <- c(levels(b16.combined@active.ident), "unknown")
b16.combined@active.ident[top_right_barcodes] <- "unknown"
```

Check new identities:

```r
table(b16.combined@active.ident)
```

![Identity Distribution](#)

Visualize on t-SNE:

```r
DimPlot(b16.combined, label = TRUE, reduction = "tsne")
```

![Unknown Cluster](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/tsne-plot-unknown-labeled.jpeg)

Perform differential expression analysis:

```r
unknown_markers <- FindMarkers(b16.combined, "unknown")

head(unknown_markers[order(unknown_markers$avg_logFC, decreasing = TRUE)[1:10], ])
```

### 5. Percentage of Positive Cells by Cluster

Calculate the fraction of cells expressing a specific gene in each cluster.

```r
get_percent_positive <- function(seuratobj, marker, idents = "all"){
  seuratobj@meta.data$markerpos <- ifelse(seuratobj@assays$RNA@data[marker, ] > 0, 1, 0)
  seuratobj@meta.data$active.ident <- as.character(seuratobj@active.ident)
  if (idents != "all")
    seuratobj@meta.data$active.ident[!(seuratobj@meta.data$active.ident %in% idents)] <- "Other"
  summary_data <- seuratobj@meta.data %>% 
    group_by(active.ident) %>% 
    summarize(
      proportion_cd45 = n() / nrow(seuratobj@meta.data) * 100,
      proportion_markerpos = mean(markerpos) * 100
    )
  summary_data
}

get_percent_positive(b16.combined, "Ccr2")
```

```r
plot_percent_positive <- function(seuratobj, marker, idents = "all") {
  summary_data <- get_percent_positive(seuratobj, marker, idents)
  ggplot(summary_data, aes(x = active.ident, y = proportion_markerpos)) +
    geom_bar(aes(fill = active.ident), stat = "identity") +
    labs(y = paste("Percent Positive for", marker),
         x = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_percent_positive(b16.combined, "Ccr2")
```

![Percentage of Positive Cells](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/percent-positive-plot.jpeg)

### 6. Median Gene Splitting

Classify cells based on above or below median expression for a specific gene.

```r
cd8 <- subset(b16.combined, idents = "CD8")

median_split <- function(seuratobj, gene, method = "plot") {
  gene_expression <- seuratobj@assays$RNA@data[gene, ]
  median_split <- factor(ifelse(gene_expression >= median(gene_expression), paste(gene, "high", sep = "-"), paste(gene, "low", sep = "-")))
  split_name <- paste(gene, "split", sep = ".")
  seuratobj[[split_name]] <- median_split
  if (method == "plot") {
    DimPlot(seuratobj, group.by = split_name, reduction = "tsne")
  } else if (method == "differential genes") {
    seuratobj <- SetIdent(seuratobj, value = split_name)
    de_results <- FindMarkers(seuratobj, paste(gene, "high", sep = "-"), paste(gene, "low", sep = "-"), min.pct = 0.25, logfc.threshold = 0.25)
    de_results <- de_results[order(de_results$avg_logFC, decreasing = TRUE), ]
    # Export results if needed
    # write.csv(de_results, file = paste0(gene, "_high_vs_low_de_results.csv"))
    de_results
  } else {
    print("Invalid method")
  }
}

median_split(cd8, "Pdcd1")
```

![Median Split](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/median-split-by-gene-tsne.jpeg)

Analyze differentially expressed genes:

```r
cd8_pd1_de_genes <- median_split(cd8, "Pdcd1", "differential genes")
head(cd8_pd1_de_genes)
```

## Group-wise Analysis

### 1. Group-wise Cell Proportion Plot

Compare the distribution of cell types between ADAR1 KO and Control samples.

```r
get_group_proportions <- function(seuratobj, group.by = "active.ident") {
  if (group.by == "active.ident") {
    seuratobj[["active.ident"]] <- seuratobj@active.ident
  }
  total_populations <- seuratobj@meta.data %>% group_by(orig.ident) %>% summarize(total_pop = n())
  count_populations <- seuratobj@meta.data %>% group_by_at(vars(group.by, "orig.ident")) %>% summarize(n = n())
  count_populations <- left_join(count_populations, total_populations, by = "orig.ident")
  count_populations <- count_populations %>% mutate(proportion = n / total_pop)
  count_populations
}

get_group_proportions(b16.combined)
```

Visualize:

```r
plot_group_proportions <- function(seuratobj, graph.type = "dodge") {
  count_populations <- get_group_proportions(seuratobj)
  if (graph.type == "dodge") {
    ggplot(count_populations, aes(x = active.ident, y = proportion)) +
      geom_bar(aes(fill = orig.ident), stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else if (graph.type == "stacked") {
    ggplot(count_populations, aes(x = orig.ident, y = proportion)) +
      geom_bar(aes(fill = active.ident), stat = "identity", position = "fill")
  } else {
    print("Invalid graph type")
  }
}

plot_group_proportions(b16.combined, graph.type = "dodge")
```

![Group-wise Cell Proportion - Dodge](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/group-celltype-proportion-dodge.jpeg)

Stacked bar plot:

```r
plot_group_proportions(b16.combined, graph.type = "stacked")
```

![Group-wise Cell Proportion - Stacked](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/group-celltype-proportion-stacked.jpeg)

### 2. Heatmap of Sample and Cell Type Similarity

Assess similarity in cell composition between samples.

```r
plot_heatmap_proportions <- function(seuratobj, graph.type = "by.cell") {
  count_populations <- get_group_proportions(seuratobj)
  reformatted <- count_populations %>% reshape2::dcast(formula = orig.ident ~ active.ident, value.var = "proportion")
  heatmap_matrix <- as.matrix(reformatted[, -1])
  rownames(heatmap_matrix) <- reformatted[, 1]
  if (graph.type == "by.cell") {
    pheatmap::pheatmap(heatmap_matrix, scale = "column")
  } else if (graph.type == "by.sample") {
    pheatmap::pheatmap(cor(t(heatmap_matrix)))
  } else {
    print("Invalid graph type")
  }
}

plot_heatmap_proportions(b16.combined, graph.type = "by.sample")
```

![Heatmap by Sample](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/heatmap-by-sample-proportions.jpeg)

```r
plot_heatmap_proportions(b16.combined, graph.type = "by.cell")
```

![Heatmap by Cell Type](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/heatmap-by-celltype-proportions.jpeg)

### 3. Median Splitting Across Groups

Classify cells based on gene expression and compare proportions across groups.

```r
proportions_median_split <- function(seuratobj, gene) {
  gene_expression <- seuratobj@assays$RNA@data[gene, ]
  median_split <- factor(ifelse(gene_expression >= median(gene_expression), paste(gene, "high", sep = "-"), paste(gene, "low", sep = "-")))
  split_name <- paste(gene, "split", sep = ".")
  seuratobj[[split_name]] <- median_split
  
  count_populations <- get_group_proportions(seuratobj, group.by = split_name)
  
  ggplot(count_populations, aes(x = orig.ident, y = proportion)) +
    geom_bar(stat = "identity", aes(fill = split_name), position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Proportion within Group")
}

proportions_median_split(b16.combined, "Ccr2")
```

![Median Split - Ccr2](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/median-split-by-group.jpeg)

### 4. Differential Gene Expression Analysis Between Groups

Identify differentially expressed genes between ADAR1 KO and Control samples.

```r
de_groups <- c("CD8", "Monocyte")

get_DE_genes <- function(seuratobj, de_groups) {
  seuratobj[["group"]] <- str_split_fixed(seuratobj@meta.data$orig.ident, "_", 2)[, 1]
  
  for (item in de_groups) {
    de_subset <- subset(seuratobj, idents = item)
    de_subset <- SetIdent(de_subset, value = "group")
    de_results <- FindMarkers(de_subset, "ADAR", "Control", min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(de_results[order(de_results$avg_logFC, decreasing = TRUE), ], file = paste0(item, "_ADAR_vs_Control_DE.csv"))
  }
}

get_DE_genes(b16.combined, de_groups)
```

```r
cd8_de <- read.csv("output/CD8_ADAR_vs_Control_DE.csv", header = TRUE, row.names = 1)
head(cd8_de)
```

### 5. Gene Set Enrichment Analysis (GSEA) Across Groups

Evaluate the enrichment of biological pathways between groups using the `fgsea` package.

```r
library(fgsea)
```

```r
# Load pathways from MSigDB
hallmark_pathway <- gmtPathways("pathway_files/h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))
```

Prepare the ranked list:

```r
prepare_ranked_list <- function(ranked_list) { 
  if (sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(. ~ Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = TRUE), ]
  }
  ranked_list <- na.omit(ranked_list)
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

cd8_ranked_list <- prepare_ranked_list(cd8_ranked_list)
head(cd8_ranked_list)
```

Run GSEA:

```r
fgsea_results <- fgsea(
  pathways = hallmark_pathway,
  stats = cd8_ranked_list,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

fgsea_results %>% arrange(desc(NES)) %>% select(pathway, padj, NES) %>% head()
```

Plot enrichment for specific pathways:

```r
plot_enrichment <- function(pathway, ranked_list) {
  plotEnrichment(hallmark_pathway[[pathway]], ranked_list) +
    labs(title = pathway)
}

# Interferon Gamma Response
plot_enrichment("HALLMARK_INTERFERON_GAMMA_RESPONSE", cd8_ranked_list)

# IL2 STAT5 Signaling
plot_enrichment("HALLMARK_IL2_STAT5_SIGNALING", cd8_ranked_list)
```

![Interferon Gamma Response Enrichment](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/gsea-ifn-gamma-response-curve.jpeg)
![IL2 STAT5 Signaling Enrichment](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/gsea-il2-stat5-signaling-curve.jpeg)

Visualize results with a waterfall plot:

```r
waterfall_plot <- function(fsgea_results, graph_title) {
  fsgea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_", 2)[, 2]) %>%
    ggplot(aes(reorder(short_name, NES), NES)) +
      geom_bar(stat = "identity", aes(fill = padj < 0.05)) +
      coord_flip() +
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title) +
      theme(
        axis.text.y = element_text(size = 7), 
        plot.title = element_text(hjust = 1)
      )
}

waterfall_plot(fgsea_results, "Pathway Enrichment in CD8 T Cells ADAR1 KO vs Control")
```

![GSEA Waterfall Plot](https://github.com/LuisNagano/SingleCell-ImmuneCells-B16Tumor/blob/main/figures/gsea-waterfall-plot.jpeg)

---

## Conclusions

The analysis provided significant insights into the tumor microenvironment (TME) of the B16 melanoma model and the impact of ADAR1 knockout (ADAR1 KO):

- **Myeloid Dominance**: The B16 model is predominantly infiltrated by myeloid-origin cells, especially various macrophage subtypes (C1q⁺, APOE⁺, Ms4a7⁺).
- **Immune Diversity**: Multiple immune cell types are present, including various dendritic cell subtypes (cDC1, cDC2, pDC, Migratory DCs) and distinct T cell populations (CD8⁺ T cells, T-regulatory cells).
- **Effects of ADAR1 KO**: ADAR1 KO samples show an increased frequency of CD8⁺ T cells, a reduction in Arg1⁺ macrophages, and a slight increase in T-regulatory cells. Classical monocytes and certain TAM subtypes (TAM 2, TAM 3) are less abundant in ADAR1 KO tumors.
- **Enhanced Interferon Response**: GSEA indicates significant enrichment of the interferon gamma response pathways in CD8⁺ T cells from ADAR1 KO samples, suggesting a more robust cytotoxic inflammatory response.

These findings provide a fundamental understanding of the immune landscape within B16 tumors and highlight the modulatory role of ADAR1 in the composition and functional states of immune cell populations. Further investigations and comparative analyses with other tumor models may deepen these biological insights.

## Contact

For any questions, feedback, or collaboration opportunities, please contact:

**Luis Fernando Nagano**

- **Email**: nagano.luis@gmail.com
- **LinkedIn**: [Luis Fernando Nagano](www.linkedin.com/in/luis-fernando-nagano-7585b82a8)
- **GitHub**: [luisfernando](https://github.com/LuisNagano)

