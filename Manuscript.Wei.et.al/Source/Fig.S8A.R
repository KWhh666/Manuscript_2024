################################# Fig.S8A
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5
###Please refer to Extended.Data.Fig.2c.R for All_VNG.rds file generation

library(Seurat)

AllVNG <- readRDS("Rds/All_VNG.rds")
AllVNG@meta.data$NvsP <- "n"
Lung_VNG <- subset(AllVNG, idents = (17))
Lung_VNG <- NormalizeData(object = Lung_VNG)
Lung_VNG <- FindVariableFeatures(object = Lung_VNG)
Lung_VNG <- ScaleData(object = Lung_VNG)
Lung_VNG <- RunPCA(object = Lung_VNG, verbose = FALSE)
ElbowPlot(Lung_VNG, ndims=50) 
Lung_VNG <- FindNeighbors(object = Lung_VNG, dims = 1:20)
Lung_VNG <- FindClusters(object = Lung_VNG, resolution = 1.0) 
Lung_VNG <- RunUMAP(object = Lung_VNG, dims = 1:20, verbose = FALSE)

Lung_VNG@meta.data$LungNvsP <- "n"
Lung_Npy2r <- subset(Lung_VNG, Npy2r > 0.1)
Lung_P2ry1 <- subset(Lung_VNG, P2ry1 > 0.1)
Lung_VNG@meta.data[rownames(Lung_Npy2r@meta.data),]$LungNvsP <- 'Lung_Npy2r'
Lung_VNG@meta.data[rownames(Lung_P2ry1@meta.data),]$LungNvsP <- 'Lung_P2ry1'

genelist <- c("Calca", "Calcb", "Vip", "Tac1")

#heatmap
DoHeatmap(subset(Lung_VNG, subset = LungNvsP !="n"), group.by = "LungNvsP", features = genelist)