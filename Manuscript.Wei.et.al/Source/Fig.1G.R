################################# Fig.1G
#Read in the reduced_data R data file
DefaultAssay(reduced_data) <- "RNA"
genelist <- read.csv(file = "Raw_Txt/DEG_heatmap_key_All_Kcng1_padj0.05_FC0.5.csv", sep = "\t", header=TRUE)

DoHeatmap(subset(reduced_data, subset = seurat_clusters =="the_Kcng1+_cluster_selected_from_FigS1F"), group.by = "group", group.colors = c("blue", "red"), features = genelist$Genes)
DoHeatmap(subset(reduced_data, subset = seurat_clusters !="the_Kcng1+_cluster_selected_from_FigS1F"), group.by = "group", group.colors = c("blue", "red"), features = genelist$Genes)
