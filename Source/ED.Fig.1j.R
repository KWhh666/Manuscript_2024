################################# ED.Fig.1j
library(Seurat)
reduced_all <- readRDS("~/reduced_all.rds")
nodose.neuron.1 <- readRDS("~/nodose.neuron.1.rds")
nodose.neuron.1 <- UpdateSeuratObject(nodose.neuron.1)
DimPlot(nodose.neuron.1,label = T)
DefaultAssay(nodose.neuron.1)
DimPlot(reduced_all,label = T)
DefaultAssay(reduced_all)

anchors <- FindTransferAnchors(reference = nodose.neuron.1, query = reduced_all, 
                               dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = nodose.neuron.1$subpopulation.ident, 
                            dims = 1:30)
reduced_all <- AddMetaData(reduced_all, metadata = predictions)
DimPlot(reduced_all,group.by = 'predicted.id',label = T)
save(reduced_all,file="~/reduced_all_withLabelTransferSubpopulation.RData")

pdf("~/LT_UMAP_K123.pdf",width = 6.8,height = 5,useDingbats = F)
DimPlot(reduced_all,group.by = 'predicted.id',label = T)
dev.off()