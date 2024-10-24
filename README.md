# Análise de RNA-Seq de Células Imunes Infiltrantes no Tumor B16

**Autor:** Erick Lu

**Data:** 12 de fevereiro de 2020

## Sumário

- [Visão Geral](#visão-geral)
- [Instalação de Pacotes](#instalação-de-pacotes)
- [Criação do Objeto Seurat](#criação-do-objeto-seurat)
- [Filtragem de Células](#filtragem-de-células)
- [Normalização de Dados e Agrupamento](#normalização-de-dados-e-agrupamento)
- [Identificação de Tipos Celulares](#identificação-de-tipos-celulares)
- [Subconjuntos e Reagrupamento](#subconjuntos-e-reagrupamento)
- [Visualizações Nativas do Seurat](#visualizações-nativas-do-seurat)
  - [1. FeaturePlot](#1-featureplot)
  - [2. VlnPlot](#2-vlnplot)
  - [3. RidgePlot](#3-ridgeplot)
  - [4. DotPlot](#4-dotplot)
  - [5. Heatmap](#5-heatmap)
  - [6. FeatureScatter](#6-featurescatter)
- [Visualizações Personalizadas](#visualizações-personalizadas)
  - [1. Gráfico de Proporção de Células](#1-gráfico-de-proporção-de-células)
  - [2. Gráfico de Contagem de UMI](#2-gráfico-de-contagem-de-umi)
  - [3. Gráfico de Assinatura Genética](#3-gráfico-de-assinatura-genética)
  - [4. Identificação de Células por Barcode](#4-identificação-de-células-por-barcode)
  - [5. Percentual de Células Positivas por Cluster](#5-percentual-de-células-positivas-por-cluster)
  - [6. Divisão por Mediana de Genes](#6-divisão-por-mediana-de-genes)
- [Análise por Grupo](#análise-por-grupo)
  - [1. Gráfico de Proporção por Grupo](#1-gráfico-de-proporção-por-grupo)
  - [2. Heatmap de Similaridade entre Amostras e Tipos Celulares](#2-heatmap-de-similaridade-entre-amostras-e-tipos-celulares)
  - [3. Divisão por Mediana entre Grupos](#3-divisão-por-mediana-entre-grupos)
  - [4. Análise de Expressão Diferencial entre Grupos](#4-análise-de-expressão-diferencial-entre-grupos)
  - [5. Análise de Enriquecimento de Conjuntos de Genes (GSEA) entre Grupos](#5-análise-de-enriquecimento-de-conjuntos-de-genes-gsea-entre-grupos)
- [Conclusões](#conclusões)

## Visão Geral

Este documento descreve uma análise passo a passo de um conjunto de dados de RNA-Seq de célula única disponível publicamente, conforme descrito por Ishizuka et al. (2019). Utilizamos o pacote R `Seurat` para processar e analisar os dados, que foram gerados utilizando a plataforma 10X Genomics a partir de leucócitos infiltrantes em tumores de camundongos B16. O artigo pode ser acessado em: [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30559380).

O objetivo desta análise é fornecer um guia prático para pesquisadores com conhecimentos básicos em R, facilitando a compreensão e análise de experimentos de RNA-Seq de célula única. Incluímos também visualizações personalizadas para uma exploração mais aprofundada dos dados.

A versão renderizada deste notebook em R está disponível em: [Análise de RNA-Seq de Célula Única](https://erilu.github.io/single-cell-rnaseq-analysis/).

---

## Instalação de Pacotes

```r
library(Seurat)
library(tidyverse)
library(future)
```

Esta análise foi conduzida utilizando a versão 3.0.2 do `Seurat`. O pacote `future` é empregado para processamento paralelo, acelerando etapas que demandam mais tempo. O `tidyverse` fornece ferramentas eficientes para manipulação e visualização de dados.

## Criação do Objeto Seurat

Os dados podem ser obtidos em: [GEO - GSE110746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110746). Baixe os arquivos `GSE110746_barcodes.tsv.gz`, `GSE110746_genes.tsv.gz` e `GSE110746_matrix.mtx.gz`. Estes arquivos contêm as matrizes de códigos de barras celulares, nomes de genes e valores de expressão, respectivamente, gerados pelo software Cell Ranger da 10X Genomics.

```r
# Carregar dados brutos
raw_data <- Read10X(data.dir = "raw_data/haining")
```

Após a leitura dos dados, criamos um objeto Seurat:

```r
b16 <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = 200, project = "b16", names.delim = "-", names.field = 2)
```

- `min.cells`: Filtra genes expressos em pelo menos 3 células.
- `min.features`: Filtra células que expressam pelo menos 200 genes.

Este filtro inicial ajuda a remover genes raros e células de baixa qualidade.

Verifique os metadados do objeto:

```r
head(b16@meta.data)
```

```r
unique(b16@meta.data$orig.ident)
```

```r
b16
```

## Filtragem de Células

Adicionamos informações sobre genes mitocondriais e ribossomais aos metadados e visualizamos a distribuição.

```r
# Identificar genes mitocondriais e ribossomais
grep("^mt-", rownames(b16@assays$RNA@data), value = TRUE)
grep("Rps|Rpl|Mrpl|Mrps", rownames(b16@assays$RNA@data), value = TRUE)
```

```r
# Calcular porcentagens
b16[["percent.mito"]] <- PercentageFeatureSet(object = b16, pattern = "^mt-")
b16[["percent.ribo"]] <- PercentageFeatureSet(object = b16, pattern = "Rps|Rpl|Mrpl|Mrps")

# Visualizar porcentagens
plot1 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "percent.mito") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot2 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot3 <- FeatureScatter(object = b16, feature1 = "nCount_RNA", feature2 = "percent.ribo") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
```

![Gráficos de Porcentagem de Mitocôndria e Ribossomos](images/figure-markdown_github/percent-mito-ribo-nfeature-graphs-1.jpeg)

Identificamos células com altas porcentagens de genes mitocondriais e ribossomais, possivelmente células comprometidas.

```r
b16 <- subset(x = b16, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40)
```

Verifique a redução no número de células:

```r
b16
```

## Normalização de Dados e Agrupamento

Procedemos com a normalização, identificação de genes variáveis e escalonamento dos dados.

```r
b16 <- NormalizeData(b16, normalization.method = "LogNormalize", scale.factor = 10000)
b16 <- FindVariableFeatures(b16, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(b16)
plan("multiprocess", workers = 4)
b16 <- ScaleData(b16, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multiprocess", workers = 1)
```

Realizamos a Análise de Componentes Principais (PCA) para redução de dimensionalidade:

```r
b16 <- RunPCA(b16, features = VariableFeatures(object = b16), ndims.print = 1:2)
```

Avalie o número de PCs a reter utilizando o gráfico de cotovelo:

```r
ElbowPlot(object = b16, ndims = 50)
```

![Gráfico de Cotovelo](images/figure-markdown_github/scree-plot-1.jpeg)

Optamos por 40 PCs para as análises subsequentes e procedemos com o agrupamento e visualização via t-SNE:

```r
b16 <- FindNeighbors(b16, dims = 1:40)
b16 <- FindClusters(b16, resolution = c(0.2, 0.4, 0.6))
b16 <- RunTSNE(b16, dims = 1:40)

# Salvar objeto processado
# save(b16, file = "b16_all_v3.Robj")
```

Exploramos diferentes resoluções para determinar a granularidade dos clusters:

```r
r02 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.2") + labs(title = "Resolução 0.2")
r04 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.4") + labs(title = "Resolução 0.4")
r06 <- DimPlot(b16, label = TRUE, reduction = "tsne", group.by = "RNA_snn_res.0.6") + labs(title = "Resolução 0.6")

CombinePlots(plots = list(r02, r04, r06), ncol = 3)
```

![Comparação de Resoluções](images/figure-markdown_github/compare-resolutions-b16-all-1.jpeg)

Optamos pela resolução 0.6 para equilibrar especificidade e relevância biológica.

## Identificação de Tipos Celulares

Identificamos genes marcadores específicos para cada cluster:

```r
# Identificação de marcadores
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Salvar marcadores
# save(b16.markers, file = "b16_all_markers_res06_v3.Robj")
```

Visualizamos os principais marcadores:

```r
top_markers_all <- b16.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top_markers_all
```

Atribuímos tipos celulares com base nos marcadores identificados:

```r
# Definição de novos identificadores de clusters
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

![Gráfico T-SNE com Tipos Celulares](images/figure-markdown_github/labeled-tsne-plot-b16-all-1.jpeg)

## Subconjuntos e Reagrupamento

Focamos nas células imunes excluindo células não imunes como melanócitos e fibroblastos, verificando a expressão de CD45 (Ptprc):

```r
FeaturePlot(b16.labeled, "Ptprc")
```

![FeaturePlot de Ptprc](images/figure-markdown_github/ptprc-featureplot-b16-all-1.jpeg)

```r
VlnPlot(b16.labeled, "Ptprc", pt.size = 0.5) + NoLegend()
```

![Violin Plot de Ptprc](images/figure-markdown_github/ptprc-vlnplot-b16-all-1.jpeg)

Excluímos os clusters com baixa expressão de CD45:

```r
# Remover clusters não imunes (4, 7, 13, 14)
b16 <- subset(b16, idents = c(0:16)[-c(5,8,14,15)])
DimPlot(b16, label = TRUE, reduction = "tsne")
```

![Remoção de Células Não Imunes](images/figure-markdown_github/remove-non-immune-tsne-1.jpeg)

Reprocessamos os dados filtrados:

```r
# Normalizar, escalar e re-agrupar
b16 <- NormalizeData(b16, normalization.method = "LogNormalize", scale.factor = 10000)
b16 <- FindVariableFeatures(b16, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(b16)
plan("multiprocess", workers = 4)
b16 <- ScaleData(b16, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multiprocess", workers = 1)
b16 <- RunPCA(b16, features = VariableFeatures(object = b16))
b16 <- FindNeighbors(b16, dims = 1:40)
b16 <- FindClusters(b16, resolution = c(0.6))
b16 <- RunTSNE(b16, dims = 1:40)

# Salvar objeto
# save(b16, file = "b16_immune_v3.Robj")
```

Visualizamos novamente os clusters:

```r
DimPlot(b16, label = TRUE, reduction = "tsne")
```

![Gráfico T-SNE Reagrupado](images/figure-markdown_github/tsne-reclustered-b16-immune-1.jpeg)

Reidentificamos os marcadores para o subconjunto imune:

```r
# Identificação de marcadores no subconjunto imune
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Salvar marcadores
# save(b16.markers, file = "b16_immune_markers_res06_v3.Robj")
```

Atribuímos novamente os tipos celulares:

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

![Gráfico T-SNE Subconjunto Imune](images/figure-markdown_github/labeled-tsne-b16-immune-1.jpeg)

Reordenamos os clusters para melhor visualização:

```r
# Reordenar clusters
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

![Clusters Reordenados](images/figure-markdown_github/reordered-tsne-b16-immune-1.jpeg)

Agrupamos os clusters para uma visão geral:

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

# Reordenar para hierarquia mieloide e linfóide
b16.combined@active.ident <- factor(b16.combined@active.ident, levels = c("TAM", "Monocyte", "cDC1", "cDC2", "pDC", "Migratory DC", "CD8", "Treg", "NK"))
# Salvar identidades agregadas
b16.combined[["res06_aggregated"]] <- b16.combined@active.ident
```

Visualizamos os clusters agregados:

```r
DimPlot(b16.combined, label = TRUE, reduction = "tsne")
```

![Gráfico T-SNE Agrupado](images/figure-markdown_github/aggregated-tsne-b16-combined-1.jpeg)

Identificamos marcadores para os clusters agregados:

```r
# Identificação de marcadores agregados
plan("multiprocess", workers = 4)
b16.markers <- FindAllMarkers(b16.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plan("multiprocess", workers = 1)

# Salvar marcadores
# save(b16.markers, file = "b16_aggregate_markers_res06_v3.Robj")
```

```r
top_markers_aggregate <- b16.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top_markers_aggregate
```

## Visualizações Nativas do Seurat

### 1. FeaturePlot

`FeaturePlot()` permite visualizar a expressão de genes específicos no espaço reduzido do t-SNE.

```r
fp <- FeaturePlot(b16.combined, "Ccr2", pt.size = 0.6)
tp <- DimPlot(b16.combined, pt.size = 0.6, label = TRUE, reduction = "tsne")
CombinePlots(plots = list(tp, fp))
```

![FeaturePlot de Ccr2](images/figure-markdown_github/featureplot-ccr2-with-tsne-1.jpeg)

### 2. VlnPlot

`VlnPlot()` gera gráficos de violino para quantificar a expressão gênica por cluster.

```r
VlnPlot(b16.combined, "Ccr2", pt.size = 0.1)
```

![Violin Plot de Ccr2](images/figure-markdown_github/vlnplot-ccr2-1.jpeg)

### 3. RidgePlot

`RidgePlot()` oferece uma visualização alternativa baseada em densidade da expressão gênica.

```r
RidgePlot(b16.combined, "Ccr2")
```

![RidgePlot de Ccr2](images/figure-markdown_github/ridgeplot-ccr2-1.jpeg)

### 4. DotPlot

`DotPlot()` facilita a comparação de múltiplos genes através dos clusters.

```r
cluster_specific_genes <- c("C1qa", "Apoe", "Ly6c2", "Ccr2", "Xcr1", "Clec9a", "Cd209a", "Clec4a4", "Siglech", "Ccr9", "Ccr7", "Fscn1", "Cd8a", "Cd3e", "Foxp3", "Il2ra", "Ncr1", "Klra4")

DotPlot(b16.combined, features = rev(cluster_specific_genes), cols = "RdBu", plot.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

![DotPlot de Genes Específicos](images/figure-markdown_github/dotplot-cell-specific-markers-1.jpeg)

### 5. Heatmap

`DoHeatmap()` visualiza a expressão gênica em células individuais dentro dos clusters.

```r
top_b16_markers <- b16.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DoHeatmap(b16.combined, cells = WhichCells(b16.combined, downsample = 50, seed = 1), features = top_b16_markers$gene, size = 3.5, angle = 15, hjust = 0)
```

![Heatmap de Marcadores Específicos](images/figure-markdown_github/heatmap-cell-specific-markers-1.jpeg)

### 6. FeatureScatter

`FeatureScatter()` avalia a correlação entre dois genes dentro de um cluster.

```r
FeatureScatter(cd8, "Tox", "Gzmb") +
  geom_smooth(method = "lm", se = FALSE)
```

![FeatureScatter de Tox vs Gzmb](images/figure-markdown_github/featurescatter-tox-gzmb-1.jpeg)

## Visualizações Personalizadas

### 1. Gráfico de Proporção de Células

Visualize a porcentagem de células que cada cluster representa no total.

```r
plot_proportions <- function(seuratobj){
  proporcao <- table(seuratobj@active.ident) / length(seuratobj@active.ident) * 100
  proporcao <- as.data.frame(proporcao)
  ggplot(proporcao, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", aes(fill = Var1)) +
    theme_classic() +
    labs(x = NULL, y = "Porcentagem Total de Células", title = "Distribuição por Tipo Celular") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    )
}

plot_proportions(b16.combined)
```

![Proporção de Tipos Celulares](images/figure-markdown_github/plot-celltype-proportions-1.jpeg)

### 2. Gráfico de Contagem de UMI

Avalie a distribuição de UMIs por célula.

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

![Gráfico de UMI](images/figure-markdown_github/umi-ncountrna-plot-1.jpeg)

Resumo das estatísticas:

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

### 3. Gráfico de Assinatura Genética

Sobreponha pontuações de assinaturas genéticas no t-SNE.

```r
signature_genes <- c("Havcr2", "Tigit", "Lag3","Pdcd1","Tox")

signature_plot <- function(seuratobject, signature_genes){
  sig_data <- as.data.frame(Embeddings(object = seuratobject, reduction = "tsne"))
  sig_data$signature <- rowMeans(t(seuratobject@assays$RNA@scale.data[signature_genes, ]))
  ggplot(sig_data, aes(tSNE_1, tSNE_2)) +
    geom_point(aes(color = signature)) +
    scale_color_gradientn(colors = topo.colors(8)) +
    theme_classic()
}

signature_plot(b16.combined, signature_genes)
```

![Gráfico de Assinatura Genética](images/figure-markdown_github/signature-plot-1.jpeg)

### 4. Identificação de Células por Barcode

Isolar e visualizar subpopulações específicas de células.

```r
# Extrair células CD8
cd8 <- subset(b16.combined, idents = "CD8")

# Reagrupar CD8
cd8 <- FindNeighbors(cd8, dims = 1:40)
cd8 <- FindClusters(cd8, resolution = 0.2)

DimPlot(cd8, label = TRUE, reduction = "tsne")
```

![Subgrupo de CD8](images/figure-markdown_github/tsne-plot-cd8-subclustered-1.jpeg)

Extrair barcodes de um subgrupo específico:

```r
cd8_cluster1_barcodes <- WhichCells(cd8, idents = 1)
```

Realçar células no t-SNE:

```r
ta <- DimPlot(b16, label = TRUE, cells.highlight = cd8_cluster1_barcodes, reduction = "tsne") + labs(title = "Usando cells.highlight")

b16@meta.data$cd8_cluster1 <- ifelse(rownames(b16@meta.data) %in% cd8_cluster1_barcodes, TRUE, FALSE)
tb <- DimPlot(b16, group.by = "cd8_cluster1", reduction = "tsne") + labs(title = "Adicionando a meta.data e usando group.by")

CombinePlots(plots = list(ta, tb))
```

![Células Destacadas no t-SNE](images/figure-markdown_github/highlight-cells-barcode-tsne-plot-1.jpeg)

Atribuir nova identidade a células isoladas:

```r
levels(b16.combined@active.ident) <- c(levels(b16.combined@active.ident), "unknown")
b16.combined@active.ident[top_right_barcodes] <- "unknown"
```

Verificar novas identidades:

```r
table(b16.combined@active.ident)
```

![Distribuição de Identidades](#)

Visualizar no t-SNE:

```r
DimPlot(b16.combined, label = TRUE, reduction = "tsne")
```

![Cluster "Unknown"](images/figure-markdown_github/tsne-plot-unknown-labeled-1.jpeg)

Analisar expressão diferencial:

```r
unknown_markers <- FindMarkers(b16.combined, "unknown")
# save(unknown_markers, file = "unknown_cluster_markers.Robj")
```

```r
unknown_markers[order(unknown_markers$avg_logFC, decreasing = TRUE)[1:10],]
```

### 5. Percentual de Células Positivas por Cluster

Calcule a fração de células que expressam um gene específico em cada cluster.

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
         x = "Tipo Celular") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_percent_positive(b16.combined, "Ccr2")
```

![Percentual de Células Positivas](images/figure-markdown_github/percent-positive-plot-1.jpeg)

### 6. Divisão por Mediana de Genes

Classifique células com base na expressão acima ou abaixo da mediana para um gene específico.

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
    # Exportar resultados se necessário
    # write.csv(de_results, file = paste0(gene, "_high_vs_low_de_results.csv"))
    de_results
  } else {
    print("Método inválido")
  }
}

median_split(cd8, "Pdcd1")
```

![Divisão por Mediana](images/figure-markdown_github/median-split-by-gene-tsne-1.jpeg)

Analise genes diferencialmente expressos:

```r
cd8_pd1_de_genes <- median_split(cd8, "Pdcd1", "differential genes")
head(cd8_pd1_de_genes)
```

## Análise por Grupo

### 1. Gráfico de Proporção por Grupo

Compare a distribuição de tipos celulares entre amostras ADAR1 KO e Controle.

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
    print("Tipo de gráfico inválido")
  }
}

plot_group_proportions(b16.combined, graph.type = "dodge")
```

![Proporção por Grupo - Dodge](images/figure-markdown_github/group-celltype-proportion-plot-dodge-1.jpeg)

Gráfico empilhado:

```r
plot_group_proportions(b16.combined, graph.type = "stacked")
```

![Proporção por Grupo - Empilhado](images/figure-markdown_github/group-celltype-proportion-plot-stacked-1.jpeg)

### 2. Heatmap de Similaridade entre Amostras e Tipos Celulares

Avalie a similaridade na composição celular entre amostras.

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
    print("Tipo de gráfico inválido")
  }
}

plot_heatmap_proportions(b16.combined, graph.type = "by.sample")
```

![Heatmap por Amostra](images/figure-markdown_github/heatmap-by-sample-proportions-1.jpeg)

```r
plot_heatmap_proportions(b16.combined, graph.type = "by.cell")
```

![Heatmap por Tipo Celular](images/figure-markdown_github/heatmap-by-celltype-proportions-1.jpeg)

### 3. Divisão por Mediana entre Grupos

Classifique células com base na expressão de um gene e compare entre grupos.

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
    labs(x = "", y = "Proporção dentro do Grupo")
}

proportions_median_split(b16.combined, "Ccr2")
```

![Divisão por Mediana - Ccr2](images/figure-markdown_github/median-split-by-group-1.jpeg)

### 4. Análise de Expressão Diferencial entre Grupos

Identifique genes diferencialmente expressos entre amostras ADAR1 KO e Controle.

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

### 5. Análise de Enriquecimento de Conjuntos de Genes (GSEA) entre Grupos

Avalie o enriquecimento de vias biológicas entre grupos usando o pacote `fgsea`.

```r
library(fgsea)
```

```r
# Carregar vias do MSigDB
hallmark_pathway <- gmtPathways("pathway_files/h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))
```

Preparar a lista ranqueada:

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

Executar GSEA:

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

Plotar enriquecimento para vias específicas:

```r
plot_enrichment <- function(pathway, ranked_list) {
  plotEnrichment(hallmark_pathway[[pathway]], ranked_list) +
    labs(title = pathway)
}

# Resposta à interferon gama
plot_enrichment("HALLMARK_INTERFERON_GAMMA_RESPONSE", cd8_ranked_list)

# Sinalização IL2 STAT5
plot_enrichment("HALLMARK_IL2_STAT5_SIGNALING", cd8_ranked_list)
```

![Enriquecimento de Interferon Gama](images/figure-markdown_github/gsea-ifn-gamma-response-curve-1.jpeg)
![Enriquecimento de IL2 STAT5](images/figure-markdown_github/gsea-il2-stat5-signaling-curve-1.jpeg)

Visualizar resultados com um gráfico de cascata:

```r
waterfall_plot <- function(fsgea_results, graph_title) {
  fsgea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_", 2)[, 2]) %>%
    ggplot(aes(reorder(short_name, NES), NES)) +
      geom_bar(stat = "identity", aes(fill = padj < 0.05)) +
      coord_flip() +
      labs(x = "Vias Hallmark", y = "Normalized Enrichment Score", title = graph_title) +
      theme(
        axis.text.y = element_text(size = 7), 
        plot.title = element_text(hjust = 1)
      )
}

waterfall_plot(fgsea_results, "Enriquecimento de Vias em CD8 T Cells ADAR1 KO vs Controle")
```

![Waterfall Plot de GSEA](images/figure-markdown_github/gsea-waterfall-plot-1.jpeg)

---

## Conclusões

A análise revelou insights significativos sobre o microambiente tumoral (TME) do modelo de melanoma B16 e o impacto da knockout de ADAR1 (ADAR1 KO):

- **Domínio Mielóide**: O modelo B16 é predominantemente infiltrado por células de origem mielóide, especialmente diversos subtipos de macrófagos (C1q⁺, APOE⁺, Ms4a7⁺).
- **Diversidade Imune**: Diversos tipos de células imunes estão presentes, incluindo múltiplos subtipos de células dendríticas (cDC1, cDC2, pDC, Migratory DCs) e distintas populações de células T (CD8⁺ T cells, T-regulatory cells).
- **Efeitos do ADAR1 KO**: Amostras ADAR1 KO apresentam aumento na frequência de células T CD8⁺, redução em macrófagos Arg1⁺ e ligeiro aumento em células T-regulatory. Monócitos clássicos e certos subtipos de TAM (TAM 2, TAM 3) são menos abundantes em tumores ADAR1 KO.
- **Resposta de Interferon Aumentada**: GSEA indica um enriquecimento significativo das vias de resposta ao interferon gama em células T CD8⁺ de amostras ADAR1 KO, sugerindo uma resposta inflamatória citotóxica mais robusta.

Esses achados fornecem uma compreensão fundamental do panorama imune dentro dos tumores B16 e destacam o papel modulador de ADAR1 na composição e nos estados funcionais das populações celulares imunes. Investigações adicionais e análises comparativas com outros modelos de tumor podem aprofundar ainda mais esses insights biológicos.

---

Espero que este guia facilite a análise de seus próprios experimentos de RNA-Seq de célula única, fornecendo tanto uma estrutura metodológica quanto visualizações personalizadas para uma exploração mais detalhada dos dados.
