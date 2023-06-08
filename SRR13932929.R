setwd("~/Documents/GBC/Seurat")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(SeuratObject)

#Setup the Seurat Object
SRR29_obj <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932929_filtered_feature_bc_matrix.h5')
SRR29<- CreateSeuratObject(counts = SRR29_obj, project = "SRR13932929",min.cells=3,min.features=20)
SRR29

#QC and selecting cells for further analysis
# operator can add columns to object metadata. This is a great place to stash QC stats
SRR29[["percent.mt"]] <- PercentageFeatureSet(SRR29, pattern = "Mt-*")

## Visualize QC metrics as a violin plot
VlnPlot(SRR29, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SRR29, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR29, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
SRR29 <- NormalizeData(SRR29, normalization.method = "LogNormalize", scale.factor = 10000)
SRR29<- NormalizeData(SRR29)

#Identification of highly variable features (feature selection)
SRR29 <- FindVariableFeatures(SRR29, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SRR29), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SRR29)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(SRR29)
SRR29 <- ScaleData(SRR29, features = all.genes)
SRR29 <- ScaleData(SRR29, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
SRR29 <- RunPCA(SRR29, features = VariableFeatures(object = SRR29))

# Examine and visualize PCA results a few different ways
print(SRR29[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SRR29, dims = 1:2, reduction = "pca")
DimPlot(SRR29, reduction = "pca")

DimHeatmap(SRR29, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(SRR29, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
SRR29 <- JackStraw(SRR29, num.replicate = 100)
SRR29 <- ScoreJackStraw(SRR29, dims = 1:20)

JackStrawPlot(SRR29, dims = 1:15)
ElbowPlot(SRR29)

#Cluster the cells
SRR29 <- FindNeighbors(SRR29, dims = 1:10)
SRR29 <- FindClusters(SRR29, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(SRR29), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
SRR29 <- RunUMAP(SRR29, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(SRR29, reduction ="umap")

#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(SRR29, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(SRR29, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
SRR29.markers <- FindAllMarkers(SRR29, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SRR29.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(SRR29, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(SRR29, features = c("Gdf15", "Cxcl1"))

# you can plot raw counts as well
VlnPlot(SRR29, features = c("Ifrd1", "Lcp1"), slot = "counts", log = TRUE)

FeaturePlot(SRR29, features = c("Selp","Fcer1g","S100a6","Ctss","Coro1a","Lcp1","Hbb-bt","Tanc2","Adh7"))

SRR29.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SRR29, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                     "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
