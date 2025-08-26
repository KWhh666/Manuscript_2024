################################# Fig.1i
library(Seurat)

reduced_data <- readRDS("Customized directory/Manuscript.Wei.et.al/Rds/reduced_all.rds")
DimPlot(reduced_data, reduction = 'umap', label = TRUE)
