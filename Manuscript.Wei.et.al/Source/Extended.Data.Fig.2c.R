################################# Extended Data Fig.2c
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5

library(Seurat)
install.packages("devtools")
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
vagus_all <- merge(vagus010419_feature_bc_matrix, y = c(vagus051818_feature_bc_matrix, vagus070519V1_feature_bc_matrix, vagus070519V2_feature_bc_matrix, vagus071819V1_feature_bc_matrix, vagus071819V2_feature_bc_matrix, vagusRC3_feature_bc_matrix, vagusRC4_feature_bc_matrix),
                   add.cell.ids = ls()[5:12], #check ls() first
                   project = 'vagus')
view(vagus_all@meta.data)
vagus_all$sample <- rownames(vagus_all@meta.data)
vagus_all@meta.data <- separate(vagus_all@meta.data, col = 'sample', into = c('Batch', 'Barcode'),
                                sep = '_feature_bc_matrix_')
unique(vagus_all@meta.data$Batch) #sanity check

#optional, QC filtering 
vagus_all_filtered <- subset(vagus_all, subset = nCount_RNA > 800 &
                               nFeature_RNA > 500) 

#Regular Seurat object pipeline#
vagus_all <- NormalizeData(object = vagus_all)
vagus_all <- FindVariableFeatures(object = vagus_all)
vagus_all <- ScaleData(object = vagus_all)
vagus_all <- RunPCA(object = vagus_all)
ElbowPlot(vagus_all)
vagus_all <- FindNeighbors(object = vagus_all, dims = 1:20)
vagus_all <- FindClusters(object = vagus_all)
vagus_all <- RunUMAP(object = vagus_all, dims = 1:20)

#PLOT#
DimPlot(vagus_all, reduction = 'umap')

#Perform integration to correct for batch effects#
obj.list <- SplitObject(vagus_all, split.by = 'Batch')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

#Select integration features#
features <- SelectIntegrationFeatures(object.list = obj.list)

#Find integration anchors (CCA)#
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

#Integrate data#
seurat.integrated <- IntegrateData(anchorset = anchors)

#Scale data, run PCA and UMAP and visualize integrated data#
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object= seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

DimPlot(seurat.integrated, reduction = 'umap')

saveRDS(seurat.integrated, file = "All_VNG.rds")

#Generate ED Fig 2c
AllVNG <- readRDS("Rds/All_VNG.rds")
FeaturePlot(AllVNG, features = "QZ2", label = T)
Lung_VNG <- subset(x = AllVNG, idents = (17) , invert = FALSE)
Lung_VNG <- NormalizeData(object = Lung_VNG)
Lung_VNG <- FindVariableFeatures(object = Lung_VNG)
Lung_VNG <- ScaleData(object = Lung_VNG)
Lung_VNG <- RunPCA(object = Lung_VNG, verbose = FALSE)
ElbowPlot(Lung_VNG, ndims=50) 
Lung_VNG <- FindNeighbors(object = Lung_VNG, dims = 1:20)
Lung_VNG <- FindClusters(object = Lung_VNG, resolution = 1.0) 
Lung_VNG <- RunUMAP(object = Lung_VNG, dims = 1:20, verbose = FALSE)

DimPlot(Lung_VNG, reduction = 'umap')
FeaturePlot(Lung_VNG, features = c("Trpv1", "Npy2r", "P2ry1"), reduction = "umap", label = TRUE)
