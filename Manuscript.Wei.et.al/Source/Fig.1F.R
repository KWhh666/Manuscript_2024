################################# Fig.1F
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
install.packages("R.utils")
library(R.utils)

###Gunzip
setwd("where_you_save_the_filtered_data_files")
#setwd("/Users/karenwhh/Library/CloudStorage/Box-Box/HaohanWei/Thesis_research/VNG_sc_seq_w_Chang_lab/2023_KP3_KP_merged_analysis")

gunzip("2023-05-03_wt/filtered_feature_bc_matrix/barcodes.tsv.gz", remove=FALSE)
gunzip("2023-05-03_wt/filtered_feature_bc_matrix/features.tsv.gz", remove=FALSE)
gunzip("2023-05-03_wt/filtered_feature_bc_matrix/matrix.mtx.gz", remove=FALSE)

gunzip("2023-05-02_kp3/filtered_feature_bc_matrix/barcodes.tsv.gz", remove=FALSE)
gunzip("2023-05-02_kp3/filtered_feature_bc_matrix/features.tsv.gz", remove=FALSE)
gunzip("2023-05-02_kp3/filtered_feature_bc_matrix/matrix.mtx.gz", remove=FALSE)

###Read in matrix; create Seurat objects
setwd("2023-05-03_wt/filtered_feature_bc_matrix")

cts_wt <- ReadMtx(mtx = paste0('matrix.mtx'),
                  features = paste0('features.tsv'),
                  cells = paste0('barcodes.tsv'))
All_wt <- CreateSeuratObject(counts = cts_wt, features = features, barcodes = cells)

setwd("2023-05-02_kp3/filtered_feature_bc_matrix")

cts_kp3 <- ReadMtx(mtx = paste0('matrix.mtx'),
                   features = paste0('features.tsv'),
                   cells = paste0('barcodes.tsv'))
All_kp3 <- CreateSeuratObject(counts = cts_kp3, features = features, barcodes = cells)

###Sorting out neuronal cell clusters
setwd('where_you_want_to_save_the_analysis')

#Rename the cell id for each group 
All_wt@meta.data$orig.ident <- recode(All_wt@meta.data$orig.ident, 'SeuratProject' = "WT")
All_kp3@meta.data$orig.ident <- recode(All_kp3@meta.data$orig.ident, 'SeuratProject' = "KP3")

view(All_wt@meta.data)
view(All_kp3@meta.data)

#Filter mito features 
mito.features <- grep(pattern = "^mt-", x = rownames(x = All_wt), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = All_wt, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = All_wt, slot = 'counts'))
All_wt[['percent.mito']] <- percent.mito
VlnPlot(object = All_wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
All_wt <- subset(x = All_wt, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = All_wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
All_wt

mito.features <- grep(pattern = "^mt-", x = rownames(x = All_kp3), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = All_kp3, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = All_kp3, slot = 'counts'))
All_kp3[['percent.mito']] <- percent.mito
VlnPlot(object = All_kp3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
All_kp3 <- subset(x = All_kp3, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mito < 0.10)
VlnPlot(object = All_kp3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1)
All_kp3

#normalize
All_wt <- NormalizeData(All_wt)
All_wt <- FindVariableFeatures(All_wt, selection.method = "vst", nfeatures = 2000)
All_kp3 <- NormalizeData(All_kp3)
All_kp3 <- FindVariableFeatures(All_kp3, selection.method = "vst", nfeatures = 2000)

#Integrate
ref.list <- list(All_wt, All_kp3)
anchors <- FindIntegrationAnchors(object.list = ref.list,
                                  dims = 1:100)
seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:100)

seurat.integrated <- FindVariableFeatures(object = seurat.integrated)
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object= seurat.integrated)
seurat.integrated <- FindNeighbors(object = seurat.integrated, dims = 1:50)
seurat.integrated <- FindClusters(object = seurat.integrated)
KP3_merged_All <- RunUMAP(object = seurat.integrated, dims = 1:50)
DimPlot(KP3_merged_All, reduction = 'umap', label = TRUE)
DimPlot(KP3_merged_All, reduction = 'umap', split.by = 'orig.ident', label = TRUE)

#Select nodose neuron clusters
DefaultAssay(KP3_merged_All) <- "RNA"
VlnPlot(KP3_merged_All, features = c('Slc17a6','Snap25','Syn1','Phox2b','Cd74','Cldn5','Apoe','Prdm12','Trpv1'))

DefaultAssay(KP3_merged_All) <- "integrated"
KP3_merged_VSN <- subset(x = KP3_merged_All, (seurat_clusters == "5")|(seurat_clusters == "6")|(seurat_clusters == "7")|(seurat_clusters == "8")|(seurat_clusters == "9")|(seurat_clusters == "10")|(seurat_clusters == "13")|(seurat_clusters == "15")|(seurat_clusters == "16")|(seurat_clusters == "20")|(seurat_clusters == "21")|(seurat_clusters == "24"), invert = FALSE)
KP3_merged_VSN <- subset(x = KP3_merged_VSN, subset = Sprr1a < 0.1 & Ecel1 < 0.1)

#2nd round - Scale data, run PCA and UMAP
reduced_data <- FindVariableFeatures(object = KP3_merged_VSN)
reduced_data <- RunPCA(object= reduced_data)
reduced_data <- FindNeighbors(object = reduced_data, dims = 1:20)
reduced_data <- FindClusters(object = reduced_data)
reduced_data <- RunUMAP(object = reduced_data, dims = 1:20)
reduced_data <- ScaleData(reduced_data, features = rownames(reduced_data@assays$RNA@counts))
DimPlot(reduced_data, reduction = 'umap', label = TRUE)
DimPlot(reduced_data, reduction = 'umap', split.by = "orig.ident", label = TRUE)
