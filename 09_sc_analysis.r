getwd()
rm(list=ls())
dir()
gc()

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
set.seed(1)
execute_steps <- c(1,2,3,4,5)
# Define variables
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- "./initial_clustering/doublet_analysis/"
analysis_parent_folder <- "./initial_clustering/"
setwd(analysis_parent_folder)
path_to_metadata <- "hubmap_htan_metadata_atac_and_rna_final.csv"
# Define sets and locations of files for initial processing
scRNA_data_path <- "/data/"
path_to_metadata <- "metadata.csv"
scRNA_set_names <- c("NORMAL","TN", "HER2", "ER")
scRNA_sets <- 
  list(c("GSM4909253","GSM4909254","GSM4909257","GSM4909261","GSM4909263","GSM4909265","GSM4909266","GSM4909268","GSM4909270","GSM4909271","GSM4909272","GSM4909274","GSM4909276"), 
      c("GSM4909281","GSM4909282","GSM4909283","GSM4909284"), c("GSM4909289", "GSM4909290","GSM4909291","GSM4909292","GSM4909293","GSM4909294"),
c("GSM4909296","GSM4909297","GSM4909298","GSM4909299","GSM4909300","GSM4909301","GSM4909302","GSM4909303","GSM4909304","GSM4909305","GSM4909306",
"GSM4909307","GSM4909309","GSM4909311","GSM4909313","GSM4909315","GSM4909317")) 

# Define functions
seurat_standard_normalize_and_scale <- function(TNBC, cluster, cluster_resolution){
	TNBC <- NormalizeData(TNBC, normalization.method = "LogNormalize", scale.factor = 10000)
	TNBC <- FindVariableFeatures(TNBC, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(TNBC)
	TNBC <- ScaleData(TNBC, features = all.genes)
	TNBC <- RunPCA(TNBC, features = VariableFeatures(object = TNBC))
	if (cluster){
		TNBC <- FindNeighbors(TNBC, dims = 1:20)
		TNBC <- FindClusters(TNBC, resolution = cluster_resolution)
	}
	TNBC <- RunUMAP(TNBC, dims = 1:20)
	return(TNBC)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
	TNBC.data <- Read10X(data.dir = data_directory)
	currentSample <- CreateSeuratObject(counts = TNBC.data, project = project_name, min.cells = 3, min.features = 40)
	currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
	pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
	dev.off()
	pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
	dev.off()
	currentSample <- subset(currentSample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000)
	currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
	nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  
	seu_TNBC <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	print(head(seu_TNBC@meta.data))
	seu_TNBC$doublet.class <- seu_TNBC[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
	seu_TNBC[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
	pann <- grep(pattern="^pANN", x=names(seu_TNBC@meta.data), value=TRUE)
	seu_TNBC$pANN <- seu_TNBC[[pann]]
	seu_TNBC[[pann]] <- NULL
	pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_TNBC, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
	dev.off()
	seu_TNBC <- subset(seu_TNBC, subset = doublet.class != "Doublet")
	pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_TNBC, reduction = "umap", cols = c("#D51F26")))
	dev.off()
	seu_TNBC <- DietSeurat(seu_TNBC, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
	return(seu_TNBC)
}

seurat_qc_plots <- function(TNBC, sample_name){
	pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(TNBC, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()
	pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(TNBC, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()
	pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(TNBC, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
	dev.off()
}


# 1) Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
if (1 %in% execute_steps){
	setwd(individual_qc_and_dublet_plot_location)
	for (j in 1:length(scRNA_set_names)){
		samples <- scRNA_sets[[j]]
		print(paste0(scRNA_data_path, samples[1], "/"))
		data_directory <- paste0(scRNA_data_path, samples[1], "/")
		sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1])
		seu_list <- c()
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path, samples[i], "/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i]))
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = scRNA_set_names[j])
		if (j==1){
			TNBC <- current_merge
		} else if (j>1){
			TNBC <- merge(TNBC, y = current_merge, project = "full_TNBC_project")
		}
	}
	setwd(analysis_parent_folder)
	TNBC[["percent.mt"]] <- PercentageFeatureSet(TNBC, pattern = "^MT-")
}

# 2) QC
if (2 %in% execute_steps){
	# create and set working directory to save qc plots
	if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
		dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	}
	setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	seurat_qc_plots(TNBC, sample_name)
	TNBC <- subset(TNBC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)
	setwd(analysis_parent_folder)
}

# 3) Add metadata
metadata1 <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
if (3 %in% execute_steps){
	metadata <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
	meta_data_types <- colnames(metadata)
	for (i in 2:length(meta_data_types)){
		identities <- TNBC[['orig.ident']]
		for (j in 1:length(metadata$Sample)){
			identities[identities==metadata$Sample[j]] <- metadata[j,meta_data_types[i]]
		}
		TNBC <- AddMetaData(TNBC, identities$orig.ident, col.name = meta_data_types[i])
	}
}

# 4) Normalize and scale data
if (4 %in% execute_steps){
	TNBC <- seurat_standard_normalize_and_scale(TNBC, TRUE, 1.0)
}

DimPlot(TNBC, reduction = "umap", cols = paletteDiscrete(values = unique(TNBC@meta.data$seurat_clusters), set = "stallion", reverse = FALSE),, label = T)
DimPlot(TNBC, reduction = "umap", group.by = "class", split.by = "class",cols = paletteDiscrete(values = unique(TNBC@meta.data$class), set = "stallion", reverse = FALSE))

DotPlot(TNBC, features =c( 'EPCAM','CD24',"SOX4","KRT18",  'COL1A1','COL3A1','MYL9',"DCN", 'PTPRC',"CD27","CD79A","CD3D","LYZ","CPA3", 'PECAM1','ENG',"PLVAP","CDH5" 
        ),  dot.scale = 8,) + RotatedAxis()

Idents(TNBC) = "orig.ident"
new.cluster.ids <- c("Stromal", 
                     "Stromal","Epithelial","Stromal","Stromal","Epithelial","Epithelial","Immune","Epithelial",
                     "Immune","Immune", "Immune","Immune", "Epithelial", "Stromal","Epithelial","Stromal",  "Epithelial","Epithelial","Epithelial", "Epithelial", "Epithelial",
			"Epithelial","Epithelial", "Stromal", "Epithelial","Stromal", "Stromal","Stromal","Stromal", "Immune","Epithelial","Immune","Immune","Epithelial","Immune","Immune"

)

identities <- as.character(TNBC@meta.data$seurat_clusters)
for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
TNBC <- AddMetaData(TNBC, identities, col.name = "CellTypeInitial")
Idents(TNBC) = "CellTypeInitial"

DimPlot(TNBC, reduction = "umap", cols = paletteDiscrete(values = unique(TNBC@meta.data$CellTypeInitial), set = "stallion", reverse = FALSE), label = T)
DotPlot(TNBC, 
        features =c('EPCAM','CD24',"SOX4","KRT18",
                    'COL1A1','MYL9', "DCN", "ACTA2", 
                    'PTPRC',"CD27","CD3D","CD79A","LYZ"
        ), 
        assay = 'RNA', 
        dot.scale = 8) + RotatedAxis()

Idents(TNBC) = "CellTypeInitial"
cell.prop<-as.data.frame(prop.table(table(Idents(TNBC),TNBC$class)))

colnames(cell.prop)<-c("cluster","sample","proportion")
library("ggplot2")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))

cell_type_cols <-paletteDiscrete(values = unique(TNBC@meta.data$CellTypeInitial), set = "stallion", reverse = FALSE)
p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=0), 
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


