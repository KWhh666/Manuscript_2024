#New data Fig. 15c

load("env-seurat_Bo_XR_matrix.RData")

xe$cell_type <- "n"
library(dplyr)

xe@meta.data <- xe@meta.data %>%
  mutate(cell_type = case_when(
    leiden_1 == 1 ~ "endothelial",
    leiden_1 == 2 ~ "T",
    leiden_1 == 3 ~ "epithelial",
    leiden_1 == 6 ~ "macrophage_6",
    leiden_1 == 8 ~ "AT2",
    leiden_1 == 9 ~ "AT1",
    leiden_1 == 11 ~ "macrophage_11",
    leiden_1 == 12 ~ "Tumor",
    leiden_1 == 13 ~ "neutrophil_13",
    leiden_1 == 16 ~ "B",
    leiden_1 == 17 ~ "neutrophil_17",
    leiden_1 == 18 ~ "neutrophil_18",
    leiden_1 == 19 ~ "neutrophil_19",
    leiden_1 == 20 ~ "neutrophil_20",
    leiden_1 == 21 ~ "neutrophil_21",
    leiden_1 == 4 ~ "unassigned",
    leiden_1 == 5 ~ "unassigned",
    leiden_1 == 7 ~ "unassigned",
    leiden_1 == 10 ~ "unassigned",
    leiden_1 == 14 ~ "unassigned",
    leiden_1 == 15 ~ "unassigned"
  ))

library(data.table)
lsTumor <- list()

for (n in 1:7){
  print(n)
  
  path <- paste0("Tumor_ROI_coordinates_KW/batch",n)
  
  files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
  
  lsTumor[[n]] <- data.frame()
  for (file in files) {
    
    data <- as.data.frame(fread(file))
    
    section <- gsub('_cells_stats.csv|.*/','',file)
    
    data$section <- section
    
    lsTumor[[n]] <- rbind(lsTumor[[n]],data)
    
    print(section)
    
  }
  rownames(lsTumor[[n]]) <- lsTumor[[n]]$`Cell ID`
  lsTumor[[n]]$batch_selection <- paste0(n,'-',lsTumor[[n]]$section)
}

batch1 <- lsTumor[[1]]
batch2 <- lsTumor[[2]]
batch3 <- lsTumor[[3]]
batch4 <- lsTumor[[4]]
batch5 <- lsTumor[[5]]
batch6 <- lsTumor[[6]]
batch7 <- lsTumor[[7]]

library(dplyr)
xe_batch1 <- subset(xe, subset = cell_id %in% batch1$`Cell ID`)
xe_batch1$section <- batch1[xe_batch1@meta.data$cell_id,]$section
xe_batch1$batch_selection <- batch1[xe_batch1@meta.data$cell_id,]$batch_selection
xe_batch2 <- subset(xe, subset = cell_id %in% batch2$`Cell ID`)
xe_batch2$section <- batch2[xe_batch2@meta.data$cell_id,]$section
xe_batch2$batch_selection <- batch2[xe_batch2@meta.data$cell_id,]$batch_selection
xe_batch3 <- subset(xe, subset = cell_id %in% batch3$`Cell ID`)
xe_batch3$section <- batch3[xe_batch3@meta.data$cell_id,]$section
xe_batch3$batch_selection <- batch3[xe_batch3@meta.data$cell_id,]$batch_selection
xe_batch4 <- subset(xe, subset = cell_id %in% batch4$`Cell ID`)
xe_batch4$section <- batch4[xe_batch4@meta.data$cell_id,]$section
xe_batch4$batch_selection <- batch4[xe_batch4@meta.data$cell_id,]$batch_selection
xe_batch5 <- subset(xe, subset = cell_id %in% batch5$`Cell ID`)
xe_batch5$section <- batch5[xe_batch5@meta.data$cell_id,]$section
xe_batch5$batch_selection <- batch5[xe_batch5@meta.data$cell_id,]$batch_selection
xe_batch6 <- subset(xe, subset = cell_id %in% batch6$`Cell ID`)
xe_batch6$section <- batch6[xe_batch6@meta.data$cell_id,]$section
xe_batch6$batch_selection <- batch6[xe_batch6@meta.data$cell_id,]$batch_selection
xe_batch7 <- subset(xe, subset = cell_id %in% batch7$`Cell ID`)
xe_batch7$section <- batch7[xe_batch7@meta.data$cell_id,]$section
xe_batch7$batch_selection <- batch7[xe_batch7@meta.data$cell_id,]$batch_selection

xe_batch1$condition <- "PBS"
xe_batch2$condition <- "PBS"
xe_batch3$condition <- "PBS"
xe_batch4$condition <- "DT"
xe_batch5$condition <- "DT"
xe_batch6$condition <- "DT"
xe_batch7$condition <- "DT"

#merge seurat datas
merged_all <- Reduce(function(x, y) merge(x, y), list(xe_batch1, xe_batch2, xe_batch3, xe_batch4, xe_batch5, xe_batch6, xe_batch7))
table(merged_all$condition)
merged_all@assays$RNA <- JoinLayers(merged_all@assays$RNA)
merged_xe <- merged_all

merged_all_T <- subset(merged_xe, subset = cell_type == "T")

gene_of_interest <- "Ifng"
batch_column <- 'batch_selection'
merged_all_T$gene_expression <- FetchData(merged_all_T, vars = gene_of_interest)[,1]
merged_all_T$gene_positive <- merged_all_T$gene_expression > 0
table_summary <- merged_all_T@meta.data %>%
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

gene_of_interest <- "Prf1"
batch_column <- 'batch_selection'
merged_all_T$gene_expression <- FetchData(merged_all_T, vars = gene_of_interest)[,1]
merged_all_T$gene_positive <- merged_all_T$gene_expression > 0
table_summary <- merged_all_T@meta.data %>%
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






