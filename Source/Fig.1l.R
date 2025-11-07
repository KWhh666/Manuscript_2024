################################# Fig.1j
#Read in the reduced_all R data file first
DefaultAssay(reduced_all) <- "RNA"

Lung_VNG <- subset(x = reduced_all, seurat_clusters == 6)
nonLung_VNG <- subset(x = reduced_all, subset = Kcng1 == 0)

gene_list <- read.csv(file = "Customized directory/Manuscript.Wei.et.al/Raw_Txt/DEG_heatmap_key_All_Kcng1_padj0.05_FC1.csv", sep = "\t", header=TRUE)
my_genes <- gene_list$gene
my_genes <- my_genes[my_genes %in% rownames(Lung_VNG)]
my_genes <- my_genes[my_genes %in% rownames(nonLung_VNG)]
print(my_genes)

Lung_VNG <- ScaleData(Lung_VNG, features = my_genes)
nonLung_VNG <- ScaleData(nonLung_VNG, features = my_genes)

Lung_VNG@meta.data$group <- factor(
  Lung_VNG@meta.data$group,
  levels = c("WT", "KP3")  # ← specify your desired order here
)

nonLung_VNG@meta.data$group <- factor(
  nonLung_VNG@meta.data$group,
  levels = c("WT", "KP3")  # ← specify your desired order here
)

DoHeatmap(
  object = Lung_VNG,
  features = my_genes,
  group.by = "group",      # or "seurat_clusters", "condition", etc.
  label = TRUE
) +
  scale_fill_gradientn(colors = c("magenta", "black", "yellow"))

DoHeatmap(
  object = nonLung_VNG,
  features = my_genes,
  group.by = "group",      # or "seurat_clusters", "condition", etc.
  label = TRUE
) +
  scale_fill_gradientn(colors = c("magenta", "black", "yellow"))