setwd("~/Documents/GBC/Seurat")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(SeuratObject)


#Setup the Seurat Object
SRR28_obj <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932928_filtered_feature_bc_matrix.h5')
SRR28<- CreateSeuratObject(counts = SRR28_obj, project = "SRR13932928",min.cells=3,min.features=20)
SRR28


#QC and selecting cells for further analysis
# operator can add columns to object metadata. This is a great place to stash QC stats
SRR28[["percent.mt"]] <- PercentageFeatureSet(SRR28, pattern = "Mt-*")

## Visualize QC metrics as a violin plot
VlnPlot(SRR28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SRR28, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR28, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
SRR28 <- NormalizeData(SRR28, normalization.method = "LogNormalize", scale.factor = 10000)
SRR28<- NormalizeData(SRR28)

#Identification of highly variable features (feature selection)
SRR28 <- FindVariableFeatures(SRR28, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SRR28), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SRR28)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(SRR28)
SRR28 <- ScaleData(SRR28, features = all.genes)
SRR28 <- ScaleData(SRR28, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
SRR28 <- RunPCA(SRR28, features = VariableFeatures(object = SRR28))

# Examine and visualize PCA results a few different ways
print(SRR28[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SRR28, dims = 1:2, reduction = "pca")
DimPlot(SRR28, reduction = "pca")

DimHeatmap(SRR28, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(SRR28, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
SRR28 <- JackStraw(SRR28, num.replicate = 100)
SRR28 <- ScoreJackStraw(SRR28, dims = 1:20)

JackStrawPlot(SRR28, dims = 1:15)
ElbowPlot(SRR28)

#Cluster the cells
SRR28 <- FindNeighbors(SRR28, dims = 1:10)
SRR28 <- FindClusters(SRR28, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(SRR28), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
SRR28 <- RunUMAP(SRR28, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(SRR28, reduction ="umap")


#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(SRR28, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(SRR28, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
SRR28.markers <- FindAllMarkers(SRR28, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SRR28.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(SRR28, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(SRR28, features = c("Gdf15", "Cxcl1"))

# you can plot raw counts as well
VlnPlot(SRR28, features = c("Ifrd1", "Lcp1"), slot = "counts", log = TRUE)

FeaturePlot(SRR28, features = c("Selp","Fcer1g","S100a6","Ctss","Coro1a","Lcp1","Hbb-bt","Tanc2","Adh7"))

SRR28.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SRR28, features = top10$gene) + NoLegend()


#Assigning cell type identity to clusters
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                     "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#saveRDS(pbmc, file = "../output/pbmc3k_final.rds")