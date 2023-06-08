setwd("~/Documents/GBC/Seurat")
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(dplyr)

#install.packages("devtools")
#devtools::install_github("timoast/signac")


my.list<-list(SRR27,SRR28,SRR29,SRR30)
features <- SelectIntegrationFeatures(object.list = my.list)

rca.anchors <- FindIntegrationAnchors(object.list = my.list, anchor.features = features)
# this command creates an 'integrated' data assay
rca.combined <- IntegrateData(anchorset = rca.anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(rca.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
rca.combined <- ScaleData(rca.combined, verbose = FALSE)
rca.combined <- RunPCA(rca.combined, npcs = 30, verbose = FALSE)
rca.combined <- RunUMAP(rca.combined, reduction = "pca", dims = 1:30)
rca.combined <- FindNeighbors(rca.combined, reduction = "pca", dims = 1:30)
DimPlot(rca.combined,group.by="orig.ident")
rca.combined <- FindClusters(rca.combined, resolution = 0.2)
DimPlot(rca.combined)


#rca.combined@meta.data
# Visualization
#p1 <- DimPlot(rca.combined, reduction = "umap", group.by = "orig.ident")
#p2 <- DimPlot(rca.combined, reduction = "umap", label = TRUE, repel = TRUE)
#p1 + p2

#Find all markers
rca.marker.genes <- FindAllMarkers(rca.combined) 

# Feature plot - visualize feature expression in low-dimensional space

#Endothelial markers
FeaturePlot(rca.combined, features = c("Pecam1","Cdh5", "Icam2", "Tie1"), min.cutoff = "q10", max.cutoff = "q90")

#EndMT
FeaturePlot(rca.combined, features = c("Tagln","Cnn1", "Acta2", "Snai1"), min.cutoff = "q10", max.cutoff = "q90")

#EndHT - Didn't find EPRC, Tie2 -> Tek
FeaturePlot(rca.combined, features = c("Sox7","Sox17", "Gata2", "Kit", "Notch1","EPRC" , "Tek", "Bmp4"), min.cutoff = "q10", max.cutoff = "q90")

#ESC/EPC - CD157 -BST1???? / Sca1 -> Ly6a
FeaturePlot(rca.combined, features = c("Sox7","Sox17", "Gata2", "Kit", "Notch1","EPRC" , "Tek", "Bmp4"), min.cutoff = "q10", max.cutoff = "q90")

#APC 
FeaturePlot(rca.combined, features = c("H2-Aa","H2-Ab1", "H2-Eb1", "Cd74"), min.cutoff = "q10", max.cutoff = "q90")

#EndICLT - No results for Lyz2
FeaturePlot(rca.combined, features = c("C1qa","C1qb","C5ar1","Tnf", "Lyz2"), min.cutoff = "q10", max.cutoff = "q90")

##############################################


#rca.combined.Integrated <-  (object = rca.combined@assays$RNA@data, species = 'Mouse', tissue = c('Heart', 'Heart muscle'))
#gene_umap <- DimPlot(clu_gene_ann_TibiaIntegrated, reduction = "umap", label = TRUE, split.by = "orig.ident")
rca.combined.Integrated <- findcelltype(object = rca.combined.markers)

rca.combined.Integrated <- findcelltype (object = rca.combined.Integrated)
#DimPlot(rca.combined.Integrated)

#gene_umap <- DimPlot(rca.combined.Integrated, reduction = "umap", label = TRUE, split.by = "orig.ident")

#cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$tissue %in% c("Heart", "Heart muscle"), ]
#rca.combined.markers <- findmarkergene(object = rca.combined.markers, species = "Mouse",marker = cellmatch_new, tissue = "Heart", "Heart muscle")
#rca.combined.markers <- findcelltype(rca.combined.markers)

#gene_umap <- DimPlot(rca.combined.markers, reduction = "umap", label = TRUE, split.by = "orig.ident")



#############################
rca.combined1 <- rca.combined
rca.combined <- RenameIdents(object = rca.combined, '0' = 'E2', '1' = 'E8', '2' = 'E4', '3' = 'E1', '4' = 'SMCs', '5' = 'Mo1', '6' = 'E3', '7' = 'Fibro', '8' = 'E5', '9' = 'Mo2', '10' = 'Mo4', '11' = 'DC', '12' = 'T', '13' = 'E6', '14' = 'Mo3', '15' = 'E7')
rca.combined1 <- rca.combined
rca.combined <- rca.combined1

DimPlot(rca.combined, pt.size = 2.5, label = TRUE, label.size = 5, ncol=2) + NoLegend()
DimPlot(rca.combined, pt.size = 2.5, label = TRUE, label.size = 5, split.by = "orig.ident", ncol=2) + NoLegend()


#Generate a cluster tree
rca_combined_small<- BuildClusterTree(object = rca.combined)
PlotClusterTree(object = rca_combined_small)

#Organize clusters
levels(rca.combined@active.ident)
cluster_names <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "SMCs", "Fibro", "Mo1", "Mo2", "Mo3", "Mo4", "DC", "T")
rca.combined@active.ident <- factor(x = rca.combined@active.ident, levels = cluster_names)
cluster_names <- levels(rca.combined@active.ident) 
DimPlot(rca.combined, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

#ggsave("E1-E7.png", dpi = 1500)

#Crear factors to have levels in orig.ident
my_levels <- c("2D-R", "2D-L", "2W-R", "2W-L")
rca.combined@meta.data$orig.ident <- factor(x = rca.combined@meta.data$orig.ident, levels = my_levels)

#Color clusters ith custome colors
kolore <- c("#E6194B", "#04FA1C", "#481930", "#29FFFF", "#B76FFF", "#D7026A", "#FF9999", "#FFA500", "#000073", "#F032E6", "#FFE011", "#FFE011", "#FFE011", "#FFE011", "#CC7000", '#7F7F7F')
DimPlot(rca.combined,  pt.size = 2.5, label = FALSE) 
DimPlot(rca.combined, reduction = "umap", pt.size = 2.5, label = FALSE) + scale_color_manual(values = kolore) + NoLegend()
#ggsave("UMAP.tiff", units="in", width=7.02, height=4.62, dpi=1200)

DimPlot(rca.combined, reduction = "umap", split.by = "orig.ident", pt.size = 2.5, label = FALSE, ncol = 2) + NoLegend()+ scale_color_manual(values = kolore) + NoLegend()

#Generate a Dot-plot
mygenes <- c("Pecam1","Cdh5", "Icam2", "Tie1","Tagln","Cnn1", "Acta2","Snai1","Sox7","Sox17","Gata2","Kit","Notch1","Tek","Bmp4", "Sox7","Sox17","Gata2","Kit","Notch1","EPRC","Tek","Bmp4","H2-Aa","H2-Ab1","H2-Eb1","Cd74","C1qa","C1qb","C5ar1","Tnf","Lyz2")


DotPlot(rca.combined, 
        features = c("Pecam1","Cdh5", "Icam2", "Tie1","Tagln","Cnn1", "Acta2","Snai1","Sox7","Sox17","Gata2","Kit","Notch1","Tek","Bmp4", "Sox7","Sox17","Gata2","Kit","Notch1","EPRC","Tek","Bmp4","H2-Aa","H2-Ab1","H2-Eb1","Cd74","C1qa","C1qb","C5ar1","Tnf","Lyz2")
)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster.markers <- FindMarkers(rca.combined, ident.1 = "E2", ident.2 = "E8", min.pct = 0.25)
head(cluster.markers)

p1<- FeaturePlot(rca.combined, features = c("Pecam1"))+ NoLegend()

#Differential analysis all versus all
rca.markers <- FindAllMarkers(rca.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Make a HEATMAP of differential genes for each cluster
#rca.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rca.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- rca.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(rca.combined, features = top10$gene) + NoLegend()+ theme(axis.text=element_text(size=5))

#Combine all 4 Mo clustes in one
new.cluster.ids <- c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7", "EC8", "SMCs", 
                     "Fibro", "Mo", "Mo", "Mo", "Mo", "DC", "T")
names(new.cluster.ids) <- levels(rca.combined)
rca.combined <- RenameIdents(rca.combined, new.cluster.ids)

#Create a subset of Endothelial population
ec <- subset(rca.combined, idents = c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7", "EC8"))

DimPlot(rca.combined, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

#ggsave("E1-E7.png", dpi = 1500)

#saveRDS(rca.combined, file = "R_project.RDS")

#####################################################

FindAllMarkers(
  rca.combined,
  #ident.1 = "EC1", ident.2 = "EC2", ident.3 = "EC3", ident.4 = "EC4", ident.5 = "EC5", ident.6 = "EC6", ident.7 = "EC7", ident.8 = "EC8",
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = "data",
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  densify = FALSE,
)

Idents(rca.combined)
rca.combined$my_meta <- paste0(rca.combined$orig.ident,"-",rca.combined$seurat_clusters) 


VlnPlot(rca.combined, type = ifelse(test = split.plot, yes = "splitViolin"))
FeaturePlot(rca.combined, features = c("Map3k7cl","Eln"))

DoHeatmap(rca.combined, features = top10$gene) + NoLegend()
