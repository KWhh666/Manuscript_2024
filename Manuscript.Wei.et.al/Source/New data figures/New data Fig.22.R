#New data Fig.22b-d

###New data Fig. 21b:
merged_tumor_xe <- readRDS("Merged.All.Tumor_from_xe")
library(dplyr)

gene_of_interest <- "Nkx2-1"
batch_column <- 'batch'
merged_xe$gene_expression <- FetchData(merged_tumor_xe, vars = gene_of_interest)[,1]
merged_xe$gene_positive <- merged_tumor_xe$gene_expression > 0
table_summary <- merged_tumor_xe@meta.data %>%
  group_by(!!sym(batch_column)) %>%
  summarise(
    total_cells = n(),
    positive_cells = sum(gene_positive),
    percent_positive = 100 * positive_cells / total_cells,
    avg_expr_positive = ifelse(
      positive_cells > 0,
      mean(gene_expression[gene_positive], na.rm = TRUE),
      NA
    ),
    .groups = 'drop'
  )

gene_of_interest <- "Ascl1"
batch_column <- 'batch'
merged_xe$gene_expression <- FetchData(merged_tumor_xe, vars = gene_of_interest)[,1]
merged_xe$gene_positive <- merged_tumor_xe$gene_expression > 0
table_summary <- merged_tumor_xe@meta.data %>%
  group_by(!!sym(batch_column)) %>%
  summarise(
    total_cells = n(),
    positive_cells = sum(gene_positive),
    percent_positive = 100 * positive_cells / total_cells,
    avg_expr_positive = ifelse(
      positive_cells > 0,
      mean(gene_expression[gene_positive], na.rm = TRUE),
      NA
    ),
    .groups = 'drop'
  )


###New data Fig. 21c:
library(data.table)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(GSEABase)
library(ggpubr)
library(dplyr)

genelist <- list(AT2 = c('Abca3', 'Nkx2-1', 'Rab27a','Muc1'),
                 AT1 = c('Hopx', 'Ager', 'Aqp5', 'Pdpn'),
                 MAE = c('Foxj1', 'Rsph1', 'Cfap65', 'Dnah12'),
                 NE = c('Ascl1', 'Ncam1', 'Syp'))

DefaultAssay(merged_tumor_xe) <- "RNA"

for (i in names(genelist)){
  genes <- list(genelist[[i]])
  merged_tumor_xe <- AddModuleScore(
    object = merged_tumor_xe,
    features = genes,
    name = i,
    assay = "RNA",
    slot = "counts",
  )
  ids <- ncol(merged_tumor_xe@meta.data)
  colnames(merged_tumor_xe@meta.data)[ids] <- substring(colnames(merged_tumor_xe@meta.data)[ids],
                                                        1, nchar(colnames(merged_tumor_xe@meta.data)[ids]))
  print(i)
}

head(merged_tumor_xe@meta.data)

merged_tumor_xe@meta.data$AT2 <- NA
merged_tumor_xe@meta.data[rownames(merged_tumor_xe@meta.data),]$AT2 <- merged_tumor_xe@meta.data$AT21
head(merged_tumor_xe@meta.data)

merged_tumor_xe@meta.data$AT1 <- NA
merged_tumor_xe@meta.data[rownames(merged_tumor_xe@meta.data),]$AT1 <- merged_tumor_xe@meta.data$AT11
head(merged_tumor_xe@meta.data)

merged_tumor_xe@meta.data$MAE <- NA
merged_tumor_xe@meta.data[rownames(merged_tumor_xe@meta.data),]$MAE <- merged_tumor_xe@meta.data$MAE1
head(merged_tumor_xe@meta.data)

merged_tumor_xe@meta.data$NE <- NA
merged_tumor_xe@meta.data[rownames(merged_tumor_xe@meta.data),]$NE <- merged_tumor_xe@meta.data$NE1
head(merged_tumor_xe@meta.data)

meta <- merged_tumor_xe@meta.data
save(meta,file = 'Addmodulescore-tumor.csv')

# Make sure condition factor levels are correct
meta$condition <- factor(meta$condition, levels = c("PBS", "DT"))

# Plot aggregated data (each point = a batch_selection)
batch_data <- meta %>%
  group_by(batch_selection, condition) %>%
  summarise(NE1 = mean(NE1, na.rm = TRUE), .groups = 'drop')

p <- ggplot(batch_data, aes(x = condition, y = NE1, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1) +
  stat_compare_means(method = "wilcox.test") +
  theme_minimal() +
  ylab("Average Module Score per Batch") +
  ggtitle("Average NE Signature by Condition (per Batch)") +
  scale_fill_manual(values = c("PBS" = "blue", "DT" = "red"))
ggsave(filename = "vlnplot_NEsignature.png", plot = p, width = 6, height = 4, dpi = 300)


batch_data <- meta %>%
  group_by(batch_selection, condition) %>%
  summarise(AT21 = mean(AT21, na.rm = TRUE), .groups = 'drop')

p <- ggplot(batch_data, aes(x = condition, y = AT21, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1) +
  stat_compare_means(method = "wilcox.test") +
  theme_minimal() +
  ylab("Average Module Score per Batch") +
  ggtitle("Average AT2 Signature by Condition (per Batch)") +
  scale_fill_manual(values = c("PBS" = "blue", "DT" = "red"))
ggsave(filename = "vlnplot_AT2signature.png", plot = p, width = 6, height = 4, dpi = 300)


batch_data <- meta %>%
  group_by(batch_selection, condition) %>%
  summarise(MAE1 = mean(MAE1, na.rm = TRUE), .groups = 'drop')

p <- ggplot(batch_data, aes(x = condition, y = MAE1, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1) +
  stat_compare_means(method = "wilcox.test") +
  theme_minimal() +
  ylab("Average Module Score per Batch") +
  ggtitle("Average MAE Signature by Condition (per Batch)") +
  scale_fill_manual(values = c("PBS" = "blue", "DT" = "red"))
ggsave(filename = "vlnplot_MAEsignature.png", plot = p, width = 6, height = 4, dpi = 300)


batch_data <- meta %>%
  group_by(batch_selection, condition) %>%
  summarise(AT11 = mean(AT11, na.rm = TRUE), .groups = 'drop')

p <- ggplot(batch_data, aes(x = condition, y = AT11, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1) +
  stat_compare_means(method = "wilcox.test") +
  theme_minimal() +
  ylab("Average Module Score per Batch") +
  ggtitle("Average AT1 Signature by Condition (per Batch)") +
  scale_fill_manual(values = c("PBS" = "blue", "DT" = "red"))
ggsave(filename = "vlnplot_AT1signature.png", plot = p, width = 6, height = 4, dpi = 300)


###New data Fig.22d:
library(Seurat)
library(harmony)
library(dplyr)

merged_all_tumor <- subset(merged_xe, subset = cell_type == "Tumor")
merged_tumor_xe <- merged_all_tumor

merged_tumor_xe[["RNA"]]$data <- merged_tumor_xe[["RNA"]]$counts
Seurat <- FindVariableFeatures(object = merged_tumor_xe)
Seurat <- ScaleData(object = Seurat)
Seurat <- RunPCA(object= Seurat, npcs = 30)
Seurat <- RunHarmony(Seurat, group.by.vars = "batch")
Seurat <- RunUMAP(Seurat, reduction = "harmony", dims = 1:20)
Seurat <- FindNeighbors(Seurat, reduction = "harmony", dims = 1:20)
Seurat <- FindClusters(Seurat, resolution = 1)
merged_all_tumor_recluster <- Seurat

library(ggplot2)
P1 <- DimPlot(merged_all_tumor_recluster, reduction = 'umap', label = T, raster = F)
ggsave("All_clusters_res=1.pdf", plot = P1, width = 8, height = 6)

P2 <- DimPlot(merged_all_tumor_recluster, reduction = 'umap', group.by = "condition", raster = F)
ggsave("By_condition_res=1.pdf", plot = P2, width = 8, height = 6)

# Count cells per cluster per condition
cell_counts <- merged_all_tumor_recluster@meta.data %>%
  count(seurat_clusters, condition, name = "count")

# Normalize by total cells per condition
total_cells <- merged_all_tumor_recluster@meta.data %>%
  count(condition, name = "total")

cell_prop <- left_join(cell_counts, total_cells, by = "condition") %>%
  mutate(percentage = count / total)

cell_prop$condition <- factor(cell_prop$condition, levels = c('PBS','DT'))

# Plot
P_yw <- ggplot(cell_prop, aes(x = seurat_clusters, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  ylab("Percentage of Cells") +
  xlab("Cluster") +
  ggtitle("Cluster Composition by Condition") +
  theme_minimal()
ggsave("Cluster Composition by Condition_res=1.pdf", plot = P_yw, width = 8, height = 6)
