setwd("~/Documents/Seurat")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(SeuratObject)

#Setup the Seurat Object
SRR_obj <- Read10X_h5(filename = '/Users/achisha_saikia/Downloads/SRR13932927_filtered_feature_bc_matrix.h5')
SRR<- CreateSeuratObject(counts = SRR_obj, project = "SRR13932927",min.cells=3,min.features=20)
SRR

#QC and selecting cells for further analysis
# operator can add columns to object metadata. This is a great place to stash QC stats
SRR[["percent.mt"]] <- PercentageFeatureSet(SRR, pattern = "Mt-*")

## Visualize QC metrics as a violin plot
VlnPlot(SRR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SRR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
SRR <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
SRR<- NormalizeData(SRR)

#Identification of highly variable features (feature selection)
SRR <- FindVariableFeatures(SRR, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SRR), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SRR)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(SRR)
SRR <- ScaleData(SRR, features = all.genes)
SRR <- ScaleData(SRR, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
SRR <- RunPCA(SRR, features = VariableFeatures(object = SRR))

# Examine and visualize PCA results a few different ways
print(SRR[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SRR, dims = 1:2, reduction = "pca")
DimPlot(SRR, reduction = "pca")

DimHeatmap(SRR, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(SRR, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
SRR <- JackStraw(SRR, num.replicate = 100)
SRR <- ScoreJackStraw(SRR, dims = 1:20)

JackStrawPlot(SRR, dims = 1:15)
ElbowPlot(SRR)

#Cluster the cells
SRR <- FindNeighbors(SRR, dims = 1:10)
SRR <- FindClusters(SRR, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(SRR), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
SRR <- RunUMAP(SRR, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(SRR, reduction ="umap")

#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(SRR, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(SRR, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
SRR.markers <- FindAllMarkers(SRR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SRR.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(SRR, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(SRR, features = c("Gdf15", "Cxcl1"))

# you can plot raw counts as well
VlnPlot(SRR, features = c("Ifrd1", "Lcp1"), slot = "counts", log = TRUE)

FeaturePlot(SRR, features = c("Selp","Fcer1g","S100a6","Ctss","Coro1a","Lcp1","Hbb-bt","Tanc2","Adh7"))

SRR.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SRR, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                     "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
