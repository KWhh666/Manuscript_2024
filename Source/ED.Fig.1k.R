################################# ED.Fig.1k
#Read in the reduced_data R data from Fig.1i and All_VNG R data
DimPlot(All_VNG,label = T)
DefaultAssay(All_VNG)
DimPlot(reduced_all)
DefaultAssay(reduced_all)

anchors <- FindTransferAnchors(reference = All_VNG, query = reduced_all, 
                               dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = All_VNG$seurat_clusters, 
                            dims = 1:30)
reduced_all <- AddMetaData(reduced_all, metadata = predictions)
DimPlot(reduced_all,group.by = 'predicted.id',label = T)