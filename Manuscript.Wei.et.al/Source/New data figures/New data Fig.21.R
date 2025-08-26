#New data Fig.21a, c-e

###load "merged_xe" from new data Fig. 15c first
merged_all_tumor <- subset(merged_xe, subset = cell_type == "Tumor")
merged_tumor_xe <- merged_all_tumor

###New data Fig. 21a:
library(ggrepel)

merged_tumor_xe@assays$RNA$data <- merged_tumor_xe@assays$RNA$counts

# Assume 'seurat_obj' is your Seurat object
set.seed(123)  # For reproducibility
Idents(merged_tumor_xe) <- "condition"
merged_tumor_xe_downsmp <- subset(x = merged_tumor_xe, downsample = 10000)
table(merged_tumor_xe_downsmp$condition)


markers <- FindMarkers(merged_tumor_xe_downsmp,
                       # test.use = "DESeq2",
                       # pseudocount.use = 0.1,
                       ident.1 = "DT",
                       ident.2 = "PBS",
                       logfc.threshold = 0,
                       min.pct = 0.01)

# filtered_markers <- markers %>% 
# filter(p_val_adj > 1e-300)

markers$gene <- rownames(markers)

# Significance labeling
markers$significance <- ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.5,
                               ifelse(markers$avg_log2FC > 0, "Up", "Down"), "NS")

# Pick top 20 most significant genes to label
top_genes <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::slice_head(n = 20)

# Volcano plot with labels
p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot - All tumor cells",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value") +
  theme(legend.position = "top") +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = Inf)

# Save
ggsave(filename = paste0("volcano_all_tumor_cells_PBSvsDT_downsampled", ".png"),
       plot = p, width = 6, height = 5)


###New data Fig. 21c:
load("~/tumorCell/env-tumor.integrated.RData")
DimPlot(tumor.integrated,label = T)
DefaultAssay(tumor.integrated) <- 'RNA'
marker.all <- FindAllMarkers(tumor.integrated,min.pct = .2)
table(marker.all$cluster)
#marker.all <- marker.all[-(grep('mt',marker.all$gene)),]
topG <- marker.all[marker.all$p_val_adj < 0.01& marker.all$avg_log2FC>1&
                     marker.all$pct.1>0.2,] %>%
  group_by(cluster) %>% top_n(n = 5, avg_log2FC) %>% as.data.frame()
DotPlot(tumor.integrated,features = unique(topG[order(topG$cluster),]$gene))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

topG <- marker.all[marker.all$p_val_adj < 0.01& marker.all$avg_log2FC>1&
                     marker.all$pct.1>0.2,] %>% as.data.frame()
topG <- topG[order(topG$avg_log2FC,decreasing = T),]
topG <- topG[order(topG$cluster,decreasing = F),]
table(topG$cluster)
topG$cellType <- paste0('cluster',topG$cluster)
write.table(topG[topG$cluster==1,],file = '~/cluster1_DEGs.tsv',quote = F,sep = '\t',row.names = F)
#addmodulescore####
DefaultAssay(xe_tumor)
for (i in c(unique(topG$cellType))){
  genes <- unique(topG[topG$cellType==i,]$gene)
  print(length(genes))
  genes <- list(genes)
  xe_tumor <- AddModuleScore(
    object = xe_tumor,
    features = genes,
    name = i)
  ids <- ncol(xe_tumor@meta.data)
  colnames(xe_tumor@meta.data)[ids] <- substring(colnames(xe_tumor@meta.data)[ids],
                                                 1, nchar(colnames(xe_tumor@meta.data)[ids]) - 1)
  print(i)
}
head(xe_tumor@meta.data)
VlnPlot(xe_tumor,features = paste0('cluster',0:6),pt.size = 0,ncol = 7,group.by = 'cell_type')
meta <- xe_tumor@meta.data
#save(meta,file = paste0('~/tumorCell/karen/env-addmodulescore-allDEGs.RData'))

table(meta$batch)
meta$condition <- 'PBS'
meta[meta$batch %in% c(4,5,6,7),]$condition <- 'DT'
table(meta$condition)
head(meta)
meta$condition <- factor(meta$condition,levels = c('PBS','DT'))
library(ggpubr)          # provides stat_compare_means()
cmp <- combn(levels(meta$condition), 2, simplify = FALSE)
dat <- meta[,c(paste0('cluster',0:6),'condition')]
dat <- reshape2::melt(dat)

pdf(file = '~/cluster1_AMS.pdf',width = 3,height = 3,useDingbats = F)
ggplot(dat[dat$variable=='cluster1',],aes(condition,value,fill=condition))+
  geom_boxplot(
    outlier.shape  = NA,            # hide outliers (already shown by violin)
    color          = "black") +
  facet_wrap(~variable,ncol = 4)+
  stat_compare_means(
    comparisons = cmp,              # list of pairwise comparisons
    method      = "wilcox.test",    # same as Seurat's default
    label       = "p.format"
  )
dev.off()


###New data Fig. 21d:
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Vector of genes
selected_genes <- c("E2f4", "Id2", "Xiap", "Hnf4a")

# Pull expression matrix and metadata
expr_matrix <- GetAssayData(merged_tumor_xe, slot = "counts")  # log-normalized data
meta <- merged_tumor_xe@meta.data
meta$cell <- rownames(meta)

# Melt data and merge with metadata
library(reshape2)
df_long <- as.data.frame(t(as.matrix(expr_matrix[selected_genes, ])))
df_long$cell <- rownames(df_long)
df_long <- merge(df_long, meta, by = "cell")

# Now reshape the data for ggplot
df_long_melted <- reshape2::melt(df_long,
                                 id.vars = c("cell", "condition", "batch_selection"), 
                                 variable.name = "gene",
                                 value.name = "expression")
df_long_melted$expression <- as.numeric(df_long_melted$expression)

df_summary <- df_long_melted %>%
  group_by(gene, batch_selection, condition) %>%
  summarise(avg_expr = mean(expression), .groups = "drop")

df_summary$condition <- factor(df_summary$condition, levels = c("PBS", "DT"))

# Define custom colors for the conditions
custom_colors <- c("PBS" = "blue", "DT" = "red")  # Replace with your actual condition names
jitter_colors <- c("PBS" = "black", "DT" = "black")  # Replace with your actual condition names

for (gene in selected_genes) {
  p <- ggplot(df_summary %>% filter(gene == !!gene), aes(x = condition, y = avg_expr, fill = condition)) +
    geom_violin(trim = FALSE, width = 0.9, alpha = 0.5) +
    geom_jitter(width = 0.1, size = 1.5, aes(color = condition), show.legend = FALSE) +  # optional: show each sample
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = jitter_colors) +
    theme_minimal() +
    ggtitle(gene) +
    ylab("Average expression per sample")
  
  ggsave(filename = paste0("vlnplot_", gene, ".png"), plot = p, width = 6, height = 4, dpi = 300)
}


###New data Fig. 21e:
library(dplyr)

gene_of_interest <- "Irf1"
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

gene_of_interest <- "Bak1"
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

gene_of_interest <- "Mlkl"
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