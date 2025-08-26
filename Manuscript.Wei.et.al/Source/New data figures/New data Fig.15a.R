#New data Fig. 15a

library(gridExtra)
library(ggplot2)
library(dplyr)

all <- read.csv("cluster_info_for_selected_cells_in_tumors.csv")

head(all)
datCellNum <- as.data.frame(table(all$leiden_1,all$section,all$batch,all$condition))

colnames(datCellNum) <- c('leiden', 'section','batch','condition','cellNum')
head(datCellNum)

datSecArea <- all %>% group_by(batch,condition,section) %>% summarise(AreaSum=sum(Area)) %>% as.data.frame()
head(datSecArea)

dim(datSecArea)
hist(log2(datSecArea$AreaSum))

datMerge <- merge(datCellNum,datSecArea,by=c('section','batch','condition'))
datMerge$CellNumPerArea <- datMerge$cellNum/datMerge$AreaSum
datMerge <- datMerge[order(datMerge$leiden),]
datMerge$cond <- factor(datMerge$cond, levels = c('PBS','DT'))

i='2'
mean_value <- datMerge[datMerge$leiden==i,] %>%
  group_by(cond) %>%
  summarise(mean_value = mean(CellNumPerArea))
t_test_result <- t.test(CellNumPerArea ~ cond, data = datMerge[datMerge$leiden==i,])
p_value <- t_test_result$p.value
x=datMerge[datMerge$leiden==i,]$CellNumPerArea
nmax <- as.numeric(quantile(x,probs=c(.99)))
temp <- data.frame(cluster=i, mean_DT=as.numeric(mean_value[mean_value$cond=='DT',2]),
                   mean_PBS=as.numeric(mean_value[mean_value$cond=='PBS',2]),
                   FC =as.numeric(mean_value[mean_value$cond=='DT',2])/as.numeric(mean_value[mean_value$cond=='PBS',2]),
                   p_value=p_value)

ggplot(datMerge[datMerge$leiden==i,],aes(cond,CellNumPerArea,fill=cond))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values = c("#0000FF", "#FF0000"))+
  ggtitle(label = paste0('leiden=1, cluster ',i))+
  scale_y_continuous(limits = c(0, nmax))+
  geom_text(data = mean_value, aes(x = cond, y = mean_value, label = paste("Mean =", round(mean_value, 4))), vjust = -1, color = "magenta") +
  annotate("text", x = 1.5, y = nmax, label = paste("p-value =", format(p_value, digits = 2)), size = 5, color = "red")
