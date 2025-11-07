################################# Fig.1i
reduced_all <- readRDS("Customized directory/Manuscript.Wei.et.al/Rds/reduced_all.rds")
DimPlot(reduced_all, reduction = 'umap', group.by = "orig.ident", shuffle = T)
