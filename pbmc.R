setwd("~/Documents/GBC/Seurat")

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

saveRDS(pbmc, file = "pbmc_8EC.rds")


#Read H5 files for each one 
dr <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932927_filtered_feature_bc_matrix.h5')
dl <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932929_filtered_feature_bc_matrix.h5')
wr <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932928_filtered_feature_bc_matrix.h5')
wl <- Read10X_h5(filename = '/Users/achisha_saikia/Documents/GBC/Seurat/SRR13932930_filtered_feature_bc_matrix.h5')


#Creat an object for each one
dr.object <- CreateSeuratObject(counts = dr, project = "2D-R")
dl.object <- CreateSeuratObject(counts = dl, project = "2D-L")
wr.object <- CreateSeuratObject(counts = wr, project = "2W-R")
wl.object <- CreateSeuratObject(counts = wl, project = "2W-L")

#Combine all the objects
pbmc <- merge(dr.object, y = c(dl.object, wr.object, wl.object), add.cell.ids = c("2dr", "2dl", "2wr", "2wl"), project = "scRNAseq")
pbmc
head(colnames(pbmc))
tail(colnames(pbmc))
unique(sapply(X = strsplit(colnames(pbmc), split = "_"), FUN = "[", 1))
table(pbmc$orig.ident, pbmc@active.ident)

# Inspect metadata 
(pbmc@meta.data)
levels(pbmc@active.ident)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "Mt-*")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7600 & percent.mt < 10)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
CombinePlots(plots = list(plot1, plot2))

#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 7, pt.size = 2.5)+NoLegend()
DimPlot(pbmc, reduction = "umap", pt.size = 2.5)+NoLegend()
DimPlot(pbmc, reduction = "umap", pt.size = 2.5, label = FALSE, label.size = 7, split.by = "orig.ident", ncol=2)
DimPlot(pbmc, reduction = "umap", pt.size = 2.5, label = TRUE, label.size = 7, split.by = "orig.ident", ncol=2)+NoLegend()

#SAFE OBJECT
saveRDS(pbmc, file = "pbmc_8EC.rds")

######################################################################
######################################################################

#Finding differentially expressed features (cluster biomarkers)
#Find all markers of cluster 1
cluster1<- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1, n = 20) 
#Find all markers of cluster 0
cluster0<- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0, n = 20) 
#Find all markers of cluster 8
cluster8<- FindMarkers(pbmc, ident.1 = 8, min.pct = 0.25)
head(cluster8, n = 20) 
#Find all markers of cluster 3
cluster3<- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3, n = 20) 

#Plot certain genes ECs
VlnPlot(pbmc, features = c("Klf2", "Klf4", "Heg1"))
VlnPlot(pbmc, features = c("Cdh5", "Pecam1", "Icam2", "Cldn5"), ncol = 2)

#UMAP Gene plots ECs
FeaturePlot(pbmc, features = c("Cdh5", "Pecam1", "Icam2", "Cldn5"), ncol=2) 
FeaturePlot(pbmc, features = c("Klf2","Klf4", "Klk10", "Nos3"), ncol=2)

#Plot certain genes SMCs
VlnPlot(pbmc, features = c("Cnn1", "Myl9", "Speg", "Myh11"), ncol = 2)
VlnPlot(pbmc, features = c("Snai1", "Snai2", "Tagln", "Twist1", "Acta2", "Zeb1"))

#UMAP Gene plots SMCs
FeaturePlot(pbmc, features = c("Cnn1", "Myl9", "Speg","Myh11"), ncol=2)
VlnPlot(pbmc, features = c("Acta2"), ncol=2)

#Plot certain genes Fibroblast
VlnPlot(pbmc, features = c("Medag", "Tcf21", "Dcn", "Pdpn"), ncol=2)

#UMAP Gene plots Fibroblast
FeaturePlot(pbmc, features = c("Medag", "Tcf21", "Dcn", "Pdpn"), ncol=2)

#Plot certain genes Mon/Mac
VlnPlot(pbmc, features = c("C1qa", "C5ar1", "C1qb", "C1qc"), ncol=2)

#UMAP Gene plots Mon/Mac
FeaturePlot(pbmc, features = c("C1qa", "C5ar1", "C1qb", "C1qc"), ncol=2)

#Plot certain genes DC
VlnPlot(pbmc, features = c("Ccr7", "Mmp25", "Flt3", "Kit"), ncol=2)

#UMAP Gene plots DC
FeaturePlot(pbmc, features = c("Ccr7", "Mmp25", "Flt3", "Kit"), ncol=2)

#Plot certain genes T cells
VlnPlot(pbmc, features = c("Cd3e", "Itk","Cd3e", "Itk"), ncol=2)

#UMAP Gene plots T cells
FeaturePlot(pbmc, features = c("Cd3e", "Itk", "Cd3e", "Itk"), ncol=2)

#########################################################################
#########################################################################
pbmc1 <- pbmc
#Label the number of cluster with Cell Name  
pbmc <- RenameIdents(object = pbmc, '0' = 'E2', '1' = 'E8', '2' = 'E4', '3' = 'E1', '4' = 'SMCs', '5' = 'Mo1', '6' = 'E3', '7' = 'Fibro', '8' = 'E5', '9' = 'Mo2', '10' = 'Mo4', '11' = 'DC', '12' = 'T', '13' = 'E6', '14' = 'Mo3', '15' = 'E7')
pbmc1 <- pbmc
pbmc <- pbmc1
DimPlot(pbmc, pt.size = 2.5, label = TRUE, label.size = 5, ncol=2) + NoLegend() 
DimPlot(pbmc, pt.size = 2.5, label = TRUE, label.size = 5, split.by = "orig.ident", ncol=2) + NoLegend()

#Generate a cluster tree
pbmc_small<- BuildClusterTree(object = pbmc)
PlotClusterTree(object = pbmc_small)

#Organize clusters
levels(pbmc@active.ident)
cluster_names <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "SMCs", "Fibro", "Mo1", "Mo2", "Mo3", "Mo4", "DC", "T")
pbmc@active.ident <- factor(x = pbmc@active.ident, levels = cluster_names)
cluster_names <- levels(pbmc@active.ident) 
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

#Crear factors to have levels in orig.ident
my_levels <- c("2D-R", "2D-L", "2W-R", "2W-L")
pbmc@meta.data$orig.ident <- factor(x = pbmc@meta.data$orig.ident, levels = my_levels)


#Color clusters ith custome colors
kolore <- c("#E6194B", "#04FA1C", "#481930", "#29FFFF", "#B76FFF", "#D7026A", "#FF9999", "#FFA500", "#000073", "#F032E6", "#FFE011", "#FFE011", "#FFE011", "#FFE011", "#CC7000", '#7F7F7F')
DimPlot(pbmc,  pt.size = 2.5, label = FALSE) 
DimPlot(pbmc, reduction = "umap", pt.size = 2.5, label = FALSE) + scale_color_manual(values = kolore) + NoLegend()
ggsave("UMAP.tiff", units="in", width=7.02, height=4.62, dpi=1200)

DimPlot(pbmc, reduction = "umap", split.by = "orig.ident", pt.size = 2.5, label = FALSE, ncol = 2) + NoLegend()+ scale_color_manual(values = kolore) + NoLegend()
#ggsave("UMAP.tiff", units="in", width=7.26, height=4.76, dpi=1200)

#Generate a Dot-plot
mygenes <- c("Itk", "Mmp25", "Ccr7", "C5ar1", "C1qc", "C1qb", "C1qa", "Lum", "Dpep1", "Medag", "Speg", "Myh11", "Cnn1", "Tie1", "Icam2", "Cdh5", "Pecam1")

DotPlot(pbmc, 
        features = c("Itk", "Mmp25", "Ccr7","C5ar1", "C1qc", "C1qb", "C1qa", "Lum", "Dpep1", "Medag", "Speg", "Myh11", "Cnn1","Tie1","Icam2","Cdh5", "Pecam1")
)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Feature-plot of some key genes
p1<- FeaturePlot(pbmc, features = c("Pecam1"))+ NoLegend()
p2<- FeaturePlot(pbmc, features = c("Myh11"))+ NoLegend()
p3<- FeaturePlot(pbmc, features = c("Dpep1"))+ NoLegend()
p4<- FeaturePlot(pbmc, features = c("C1qa"))+ NoLegend()
p5<- FeaturePlot(pbmc, features = c("Ccr7"))+ NoLegend()
p6<- FeaturePlot(pbmc, features = c("Itk"))+ NoLegend()

(p1|p2|p3)/(p4|p5|p6)

#Differentia analysis E2 vs other
cluster.markers <- FindMarkers(pbmc, ident.1 = "E2", ident.2 = "E8", min.pct = 0.25)
head(cluster.markers)
write.csv(cluster.markers, "scRNAseq differential analysis E2vsE8.csv")

#Differential analysis all versus all
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, "scRNAseq differential analysis.csv")

#Make a HEATMAP of differential genes for each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()+ theme(axis.text=element_text(size=5))
ggsave("UMAP.tiff", units="in", width=19, height=12.5, dpi=1200)

#Combine all 4 Mo clustes in one
new.cluster.ids <- c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7", "EC8", "SMCs", 
                     "Fibro", "Mo", "Mo", "Mo", "Mo", "DC", "T")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

#Create a subset of Endothelial population
ec <- subset(pbmc, idents = c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7", "EC8"))

DimPlot(ec, reduction = "umap", label = TRUE, label.size = 7, pt.size = 2.5)+NoLegend()

DimPlot(ec, pt.size = 2.5, label = TRUE, label.size = 5, split.by = "orig.ident", ncol=2) + NoLegend()

