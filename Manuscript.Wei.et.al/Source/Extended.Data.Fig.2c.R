################################# Extended Data Fig.2c
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5
###All_VNG.rds file has been uploaded here: https://upenn.box.com/s/oyidlmnw6z92mt0xudfhp4rnc1uv6k4h
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5

library(Seurat)
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
