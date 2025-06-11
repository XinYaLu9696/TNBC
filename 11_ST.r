dev.off()
rm(list=ls())
dir()
gc()

getwd()
setwd("ST_V1")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")
set.seed(1)
library(patchwork)
library(dplyr)
library(hdf5r)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(BayesSpace)
library(ArchR)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887") 

TNBC<-Load10X_Spatial(
  data.dir = '/data/V1/', 
  assay = "Spatial",
  slice = "slice1",  
  filter.matrix = TRUE, 
  to.upper = TRUE 
  )

TNBC
TNBC@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(TNBC@images[["slice1"]]@coordinates[["tissue"]])
TNBC@images[["slice1"]]@coordinates[["row"]] <- as.integer(TNBC@images[["slice1"]]@coordinates[["row"]])
TNBC@images[["slice1"]]@coordinates[["col"]] <- as.integer(TNBC@images[["slice1"]]@coordinates[["col"]])
TNBC@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(TNBC@images[["slice1"]]@coordinates[["imagerow"]])
TNBC@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(TNBC@images[["slice1"]]@coordinates[["imagecol"]])

plot1 <- VlnPlot(TNBC, features = "nCount_Spatial", pt.size = 0.5) + NoLegend()
plot2 <- SpatialFeaturePlot(TNBC, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

TNBC <- SCTransform(TNBC, assay = "Spatial", verbose = FALSE)
TNBC <- SCTransform(TNBC, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
TNBC <- NormalizeData(TNBC, verbose = FALSE, assay = "Spatial") 

TNBC <- GroupCorrelation(TNBC, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
TNBC <- GroupCorrelation(TNBC, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(TNBC, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(TNBC, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p3 <- plot_grid(p1, p2)
p3

TNBC <- RunPCA(TNBC, assay = "SCT", verbose = FALSE)  
TNBC <- FindNeighbors(TNBC, reduction = "pca", dims = 1:30)
TNBC <- FindClusters(TNBC, verbose = FALSE)
TNBC <- RunUMAP(TNBC, reduction = "pca", dims = 1:40)

p1 <- DimPlot(TNBC, reduction = "umap", label = TRUE,
              cols = paletteDiscrete(values = unique(TNBC@meta.data$seurat_clusters)))
p2 <- SpatialDimPlot(TNBC, label = TRUE, label.size = 3,
                     cols = paletteDiscrete(values = unique(TNBC@meta.data$seurat_clusters)))
p1
p2

SpatialFeaturePlot(TNBC, features = "FDCSP", pt.size.factor = 2)

markers.to.plot <- c( "KRT19","KRT18","CD24", "SOX4", 
                      "PTPRC", "PECAM1","VWF", "COL1A1","COL3A1","DCN","DES", "ACTA2","MYL9", "CD79A"
)

Idents(TNBC) = "seurat_clusters"
DotPlot(TNBC, 
        features =markers.to.plot, 
        dot.scale = 8) + RotatedAxis()

TNBC@meta.data$Celltype = "NA"
TNBC@meta.data$Celltype[which(TNBC@meta.data$seurat_clusters%in%c(1,2,5,6,8,9,10,11,13,14,15))] <- 'Epi'
TNBC@meta.data$Celltype[which(TNBC@meta.data$seurat_clusters%in%c(0,3,4,7, 12,16))] <- 'Immune_Stroma'

Idents(TNBC) = "Celltype"
p1 <- DimPlot(TNBC, reduction = "umap", label = TRUE,
              cols = paletteDiscrete(values = unique(TNBC@meta.data$Celltype)))
p2 <- SpatialDimPlot(TNBC, label = F, label.size = 3,
                     cols = paletteDiscrete(values = unique(TNBC@meta.data$Celltype)))

DotPlot(TNBC, 
        features =markers.to.plot, 
        dot.scale = 8) + RotatedAxis()
