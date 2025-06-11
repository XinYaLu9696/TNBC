getwd()
rm(list=ls())
dir()
gc()

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
library(Rmagic)
library(ggpubr)
set.seed(1)

analysis_parent_folder <- "results/epi"
setwd(analysis_parent_folder)
TNBC = readRDS("epithelial_all_samples_initial.rds")
sample_name <- "all_samples" 
execute_steps <- c(1,2,3)

vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(TNBC, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(TNBC@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(TNBC@meta.data$orig.ident))$Freq)))
	dev.off()
}

normalize_and_dim_reduce <- function (TNBC, sample_name){
	# stanfard seurat pipeline
	TNBC <- NormalizeData(TNBC, normalization.method = "LogNormalize", scale.factor = 10000)
	TNBC <- FindVariableFeatures(TNBC, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(TNBC)
	TNBC <- ScaleData(TNBC, features = all.genes)
	TNBC <- RunPCA(TNBC, features = VariableFeatures(object = TNBC))
	TNBC <- FindNeighbors(TNBC, dims = 1:30)
	TNBC <- FindClusters(TNBC, resolution = 1)
	TNBC <- RunUMAP(TNBC, dims = 1:30)
	return(TNBC)
}

plotUMAPandRunHarmony <- function(TNBC, run_harmony, version){
	if (run_harmony){
		library(harmony)
		TNBC <- RunHarmony(TNBC, "orig.ident") # will use PCA
		TNBC <- RunUMAP(TNBC, dims = 1:30, reduction = "harmony", reduction.name = "umapharmony")
		TNBC <- FindNeighbors(TNBC, reduction = "harmony", dims = 1:20)
		TNBC <- FindClusters(TNBC, resolution = 1.0)
	}
	return(TNBC)
}

# 1) Normalize and scale data
if (1 %in% execute_steps){
	TNBC <- normalize_and_dim_reduce(TNBC, sample_name)
}


# 2) Plot UMAPs and run harmony
if (2 %in% execute_steps){
	TNBC <- plotUMAPandRunHarmony(TNBC, TRUE, version = "initial")
}
DimPlot(TNBC, reduction = "umap",label = T,
        cols = paletteDiscrete(values = unique(TNBC@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
Myo=c("KRT17", "KRT14", "KRT5", "ACTA2", "MYL9", "MYLK", "MYH11","TAGLN") # myoepithelial cells
Lum=c("KRT19", "KRT18", "KRT8") # luminal cells
Hs=c("PRLR", "CITED1", "PGR", "PROM1", "ESR1")  
AV=c("MFGE8", "TF", "CSN3", "ELF5", "LTF")
Lp=c("KIT", "ALDH1A3", "CD14")
genes_to_check = list(
  Myo=Myo,
  Lum=Lum,
  Hs=Hs, 
  AV=AV,
  Lp=Lp ) 
genes_to_check = lapply(genes_to_check , str_to_upper)
p_all_markers=DotPlot(TNBC, 
                      features = genes_to_check,
                      scale = T,assay='RNA' )
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers

# 3) Remove bad clusters and redo analysis
if (3 %in% execute_steps){
	bad_clusters <- bad_clusters <- c(5)
	TNBC_new <- DietSeurat(subset(TNBC, subset = seurat_clusters %ni% bad_clusters))
	TNBC_new <- normalize_and_dim_reduce(TNBC_new, sample_name)
	TNBC_new <- FindClusters(TNBC_new, resolution = 2.0)
	TNBC <- TNBC_new
}

DimPlot(TNBC, reduction = "umap",label = T,
        cols = paletteDiscrete(values = unique(TNBC@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))


TNBC@meta.data$celltype_epi = "NA"
basal_clusters <- c(4,8,12,13,20)
basal_clusters_low <- c(28)
luminal_clusters <- c(0:3,5:7,14:16,21:24,26,29:36)
normal_clusters = c(9,10,11,17:19,25,27)
TNBC@meta.data[which(TNBC@meta.data$seurat_clusters %in% luminal_clusters),'celltype_epi'] <- "FDCSP_low_Luminal"
TNBC@meta.data[which(TNBC@meta.data$seurat_clusters %in% basal_clusters),'celltype_epi'] <- "FDCSP_high_basal"
TNBC@meta.data[which(TNBC@meta.data$seurat_clusters %in% basal_clusters_low),'celltype_epi'] <- "FDCSP_low_basal"
TNBC@meta.data[which(TNBC@meta.data$seurat_clusters %in% normal_clusters),'celltype_epi'] <- "FDCSP_medium_NEs"

p_FDCSP = FeaturePlot(TNBC,
            reduction = "umap",
            features = c("FDCSP"),
            sort.cell = TRUE,
            pt.size = 1,
            raster=FALSE,
            label = F)
p_FDCSP

pdf("./plot/p_FDCSP.pdf", width=7, height=6)
print(p_FDCSP)
graphics.off()

Idents(TNBC) = "celltype_epi"
p1 = DimPlot(TNBC, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(TNBC@meta.data$celltype_epi), set = "stallion", reverse = F))
p1

pdf("p1_epi.pdf", width=9, height=6)
print(p1)
graphics.off()

p11 = DimPlot(TNBC, reduction = "umap",label = F,split.by = "class",
        cols = paletteDiscrete(values = unique(TNBC@meta.data$celltype_epi), set = "stallion", reverse = F))
p11

pdf("p11_epi_split.pdf", width=25, height=6)
print(p11)
graphics.off()

p3 = DotPlot(TNBC, 
        features =c("FDCSP",
                    "KRT19", "KRT18", "KRT8",
                   "KRT17", "KRT14", "KRT5",  "MYL9", "MYLK"
                   ), 
        dot.scale = 8
        ) + RotatedAxis()
p3

pdf("dotplot_epi.pdf", width=8, height=4)
print(p3)
graphics.off()

p4 = FeaturePlot(TNBC,
            reduction = "umap",
            features = c("FDCSP"),
            sort.cell = TRUE,
            split.by = "class",
            pt.size = 1,
            label = F)
p4 

pdf("FeaturePlot_spp1_FDCSP.pdf", width=24, height=6)
print(p4)
graphics.off()

cell.prop<-as.data.frame(prop.table(table(TNBC@meta.data$celltype_epi,TNBC$class)))
colnames(cell.prop)<-c("cluster","sample","proportion")

p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))

cell_type_cols <-paletteDiscrete(values = unique(TNBC@meta.data$celltype_epi), set = "stallion", reverse = FALSE)
p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=1), 
        axis.ticks.length=unit(.1,"lines"),
        axis.ticks.x = element_line(size=0.05, colour = "black"),
        #axis.ticks.margin=unit(.4,"cm"),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"),
        #axis.ticks.x.bottom =  = element_line(size = 0.5),
        panel.grid = element_blank(),
        #scale_fill_distiller(palette = "Spectral"),
        #scale_fill_brewer(palette = 'Paired'),
        #scale_fill_manual(values=ccc),
        legend.position = 'right',
        legend.key.width = unit(0.2,'cm'),
        legend.key.height = unit(0.2,'cm'),
        legend.text = element_text(color = 'black',size = 16))
p

pdf("barplot_epi.pdf", width=7, height=6)
print(p)
graphics.off()

p_EGFR = FeaturePlot(TNBC,
            reduction = "umap",
            features = c("EGFR"),
            sort.cell = TRUE,
            pt.size = 1,
            raster=FALSE,
            label = F)
p_EGFR

pdf("./plot/p_EGFR_ok.pdf", width=7, height=6)
print(p_EGFR)
graphics.off()

p5 = FeaturePlot(TNBC, features = c("FDCSP", "EGFR"),             
            sort.cell = TRUE,
            pt.size = 0.8,blend = TRUE, blend.threshold=0)
p5 

pdf("./plot/FeaturePlot_FDCSP_EGFR.pdf", width=24, height=6)
print(p5)
graphics.off()

