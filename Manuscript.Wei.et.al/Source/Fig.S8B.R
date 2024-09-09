################################# Fig.S8B
###R data adapted from Zhao, Q., Yu, C.D., Wang, R. et al. A multidimensional coding architecture of the vagal interoceptive system. Nature 603, 878–884 (2022). https://doi.org/10.1038/s41586-022-04515-5
###Please refer to Extended.Data.Fig.2c.R for All_VNG.rds file generation

library(Seurat)

AllVNG <- readRDS("Rds/All_VNG.rds")
AllVNG@meta.data$NvsP <- "n"
Lung_VNG <- subset(AllVNG, idents = (17))
Lung_VNG <- NormalizeData(object = Lung_VNG)
Lung_VNG <- FindVariableFeatures(object = Lung_VNG)
Lung_VNG <- ScaleData(object = Lung_VNG)
Lung_VNG <- RunPCA(object = Lung_VNG, verbose = FALSE)
ElbowPlot(Lung_VNG, ndims=50) 
Lung_VNG <- FindNeighbors(object = Lung_VNG, dims = 1:20)
Lung_VNG <- FindClusters(object = Lung_VNG, resolution = 1.0) 
Lung_VNG <- RunUMAP(object = Lung_VNG, dims = 1:20, verbose = FALSE)

Lung_VNG@meta.data$LungNvsP <- "n"
Lung_Npy2r <- subset(Lung_VNG, Npy2r > 0.1)
Lung_P2ry1 <- subset(Lung_VNG, P2ry1 > 0.1)
Lung_VNG@meta.data[rownames(Lung_Npy2r@meta.data),]$LungNvsP <- 'Lung_Npy2r'
Lung_VNG@meta.data[rownames(Lung_P2ry1@meta.data),]$LungNvsP <- 'Lung_P2ry1'

genelist <- c("Calca", "Calcb", "Vip", "Tac1")

#violin plot
cell1 <- colnames(Lung_Npy2r)
length(cell1)
cell2 <- colnames(Lung_P2ry1)
length(cell2)
data <- Lung_VNG@assays$RNA@data
data[1:5,1:5]

datAll <- data.frame()
for (i in genelist){
  dat1 <- as.data.frame(data[i,cell1])
  head(dat1)
  colnames(dat1) <- 'Expression_level'
  dat1$type <- 'Npy2r'
  
  dat2 <- as.data.frame(data[i,cell2])
  head(dat2)
  colnames(dat2) <- 'Expression_level'
  dat2$type <- 'P2ry1'
  
  dat <- rbind(dat1,dat2)
  head(dat)
  
  dat$gene <- i
  datAll <- rbind(datAll,dat)
  datAll$gene <- factor(datAll$gene,levels = genelist)
}
head(datAll)
mytheme <- c()
library(ggpubr)
ggplot(datAll,aes(type,Expression_level,fill=type))+
  geom_violin()+
  stat_summary(fun=mean,geom = 'point',shape=95,size=3,fill='black',color='black')+
  ggtitle(label = 'pos')+guides(fill='none')+
  ylim(0,max(datAll$Expression_level))+
  theme_set(theme_bw() +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank()))+ 
  facet_wrap(~gene,nrow = 1)+
  stat_compare_means(method = "wilcox.test", label = "p.format")

