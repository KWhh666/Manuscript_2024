################################# ED.Fig.1d
library(pheatmap)

Combined_avg <- read.csv('Customized directory/Manuscript.Wei.et.al/Raw_Txt/Tumor_Bulk_Combined_gene_list_avg.csv', header = TRUE, row.names = 1)
pheatmap(log2(as.matrix(Combined_avg)+1), scale = "row", cluster_cols = F, cluster_rows = F, treeheight_row = 0)
