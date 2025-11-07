################################# ED. Fig.1l
p_clusters <- DimPlot(
  reduced_all, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, label.size = 3
) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(
  filename = file.path(out_dir, "UMAP_clusters.svg"),
  plot = p_clusters, width = 6.5, height = 5.5, device = svglite
)


library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

reduced_all <- readRDS("Customized directory/Manuscript.Wei.et.al/Rds/reduced_all.rds")

## ---- Parameters ----
GROUP_COL   <- "group"
CLUSTER_COL <- "seurat_clusters"
ASSAY_USED  <- "RNA"
PADJ_CUT    <- 0.05
LFC_CUT     <- 1
MIN_PCT     <- 0.10
LOGFC_TEST  <- 0.0
ID1         <- 'KP3'
ID2         <- 'WT'

obj <- reduced_all
DefaultAssay(obj) <- ASSAY_USED
Idents(obj) <- obj[[CLUSTER_COL]][,1]

## ---- Select two groups ----
groups_all <- unique(obj[[GROUP_COL]][,1])
if (length(groups_all) < 2) {
  stop("Fewer than two groups; cannot perform two-group DE.")
}
if (is.null(ID1) || is.null(ID2)) {
  groups_all <- sort(as.character(groups_all))
  ID1 <- groups_all[1]
  ID2 <- groups_all[2]
}
message(sprintf("Comparing groups: %s (ident.1) vs %s (ident.2)", ID1, ID2))

## ---- Helper: identify logFC column ----
get_lfc_col <- function(df) {
  if ("avg_log2FC" %in% colnames(df)) return("avg_log2FC")
  if ("avg_logFC"  %in% colnames(df)) return("avg_logFC")
  stop("Could not find logFC column (avg_log2FC / avg_logFC).")
}

## ---- DE analysis within each cluster ----
clusters <- levels(obj[[CLUSTER_COL]][,1])
res_list <- list()
summary_counts <- tibble(cluster = character(), up = integer(), down = integer(), total = integer())

for (cl in clusters) {
  cells_in_cl <- WhichCells(obj, idents = cl)
  sub <- subset(obj, cells = cells_in_cl)
  
  gs <- unique(sub[[GROUP_COL]][,1])
  if (!all(c(ID1, ID2) %in% gs)) {
    message(sprintf("Cluster %s: one of the groups (%s/%s) is missing; skip.", cl, ID1, ID2))
    summary_counts <- add_row(summary_counts, cluster = cl, up = 0L, down = 0L, total = 0L)
    next
  }
  
  Idents(sub) <- sub[[GROUP_COL]][,1]
  
  mk <- tryCatch({
    FindMarkers(
      sub, ident.1 = ID1, ident.2 = ID2, only.pos = FALSE,
      min.pct = MIN_PCT, logfc.threshold = LOGFC_TEST, test.use = "wilcox"
    ) %>% as.data.frame() %>% tibble::rownames_to_column("gene")
  }, error = function(e) {
    message(sprintf("Cluster %s: FindMarkers error; skip. Reason: %s", cl, e$message))
    NULL
  })
  
  if (is.null(mk) || nrow(mk) == 0) {
    summary_counts <- add_row(summary_counts, cluster = cl, up = 0L, down = 0L, total = 0L)
    next
  }
  
  lfc_col <- get_lfc_col(mk)
  
  mk <- mk %>%
    mutate(sig = ifelse(p_val_adj < PADJ_CUT & .data[[lfc_col]] >  LFC_CUT, "up",
                        ifelse(p_val_adj < PADJ_CUT & .data[[lfc_col]] < -LFC_CUT, "down", "ns")))
  
  up_n   <- sum(mk$sig == "up")
  down_n <- sum(mk$sig == "down")
  
  res_list[[cl]] <- mk
  summary_counts <- add_row(summary_counts,
                            cluster = cl,
                            up = as.integer(up_n),
                            down = as.integer(down_n),
                            total = as.integer(up_n + down_n))
}

## ---- UMAP centroids & labels ----
emb <- Embeddings(obj, reduction = "umap") %>% as.data.frame()
emb$cell <- rownames(emb)
emb[[CLUSTER_COL]] <- obj@meta.data[emb$cell, CLUSTER_COL, drop = TRUE]

centroids <- emb %>%
  group_by(.data[[CLUSTER_COL]]) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop") %>%
  rename(cluster = {{CLUSTER_COL}}) %>%
  left_join(summary_counts, by = "cluster") %>%
  mutate(label = sprintf("↑%d / ↓%d", up, down))

## ---- Plot with DEG counts ----
p <- DimPlot(obj, reduction = "umap", group.by = CLUSTER_COL, label = TRUE, label.size = 3) +
  ggtitle(sprintf("DEG counts per cluster (%s vs %s)", ID1, ID2)) +
  geom_label(
    data = centroids,
    aes(x = UMAP_1, y = UMAP_2, label = label),
    inherit.aes = FALSE,
    alpha = 0.85,
    label.size = 0.2,
    label.r = unit(0.15, "lines"),
    size = 3
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)

## ---- Outputs: summary table ----
summary_counts <- arrange(summary_counts, as.numeric(as.character(cluster)))
summary_counts

## ---- DEG count maps ----
library(ggplot2)
library(dplyr)
library(patchwork)

CLUSTER_COL <- "seurat_clusters"
emb <- Embeddings(reduced_all, "umap") %>% as.data.frame()
emb$cell <- rownames(emb)
emb$cluster <- reduced_all@meta.data[emb$cell, CLUSTER_COL, drop = TRUE]

summary_counts$cluster <- as.character(summary_counts$cluster)
emb$cluster <- as.character(emb$cluster)

emb2 <- left_join(emb, summary_counts, by = "cluster")

centroids2 <- emb2 %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2),
    up = mean(up), down = mean(down), .groups = "drop"
  )

## Upregulated map
p_up <- ggplot(emb2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = up), size = 0.1) +
  scale_color_gradientn(
    colours = c("#FEE0D2","#FC9272","#DE2D26","#A50F15"),
    name = "Upregulated genes"
  ) +
  geom_text(data = centroids2, aes(label = cluster), color = "black", size = 3) +
  labs(title = "Upregulated DEG counts per cluster") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## Downregulated map
p_down <- ggplot(emb2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = down), size = 0.1) +
  scale_color_gradientn(
    colours = c("#DEEBF7","#9ECAE1","#4292C6","#08519C"),
    name = "Downregulated genes"
  ) +
  geom_text(data = centroids2, aes(label = cluster), color = "black", size = 3) +
  labs(title = "Downregulated DEG counts per cluster") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf("~/VSN/p_degs_num.pdf", width = 8.2, height = 2.5)
p_up + p_down
dev.off()