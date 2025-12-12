################################# ED.Fig.2c
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5
AllVNG <- readRDS("nodose.neuron.1.rds")
AllVNG <- UpdateSeuratObject(AllVNG)
LungVNG <- subset(AllVNG, subset = subpopulation.ident %in% c("K1", "K2", "K3", "L1", "A3"))
LungVNG <- FindVariableFeatures(object = LungVNG)
LungVNG <- ScaleData(object = LungVNG)
LungVNG <- RunPCA(object = LungVNG, verbose = FALSE)
ElbowPlot(LungVNG, ndims=50) 
LungVNG <- FindNeighbors(object = LungVNG, dims = 1:20)
LungVNG <- FindClusters(object = LungVNG, resolution = 1.0) 
LungVNG <- RunUMAP(object = LungVNG, dims = 1:20, verbose = FALSE)

p1 <- DimPlot(LungVNG, reduction = 'umap')
ggsave(
  filename = "changlung_UMAP.svg",
  plot = p1,              # your ggplot object name
  width = 7,
  height = 6,
)

DefaultAssay(LungVNG) <- "RNA"
p2 <- FeaturePlot(LungVNG, features = "Trpv1", reduction = "umap")
ggsave(
  filename = "changlung_Trpv1.svg",
  plot = p2,              # your ggplot object name
  width = 7,
  height = 6,
)
p3 <- FeaturePlot(LungVNG, features = "Npy2r", reduction = "umap")
ggsave(
  filename = "changlung_Npy2r.svg",
  plot = p3,              # your ggplot object name
  width = 7,
  height = 6,
)
p4 <- FeaturePlot(LungVNG, features = "P2ry1", reduction = "umap")
ggsave(
  filename = "changlung_P2ry1.svg",
  plot = p4,              # your ggplot object name
  width = 7,
  height = 6,
)
p5 <- FeaturePlot(LungVNG, features = "Kcng1", reduction = "umap")
ggsave(
  filename = "changlung_Kcng1.svg",
  plot = p5,              # your ggplot object name
  width = 7,
  height = 6,
)