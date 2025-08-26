################################# ED.Fig.13a and b

library(data.table)
library(survival)
library(survminer)
library(ggpubr)
library(TCGAbiolinks)
library(EDASeq)
library(tidyverse)
#survival plot####
#load clinical data which was downloaded from TCGA
load('Customized directory/Manuscript.Wei.et.al/TCGA/clin.luad.RData')#clin.luad <- GDCquery_clinic('TCGA-LUAD','clinical')
#res.ssgsea is output from GenePattern(https://cloud.genepattern.org/gp/pages/index.jsf)
res.ssgsea <- as.data.frame(fread('Customized directory/Manuscript.Wei.et.al/TCGA/ssGSEA.gct'))
tmp <- data.frame(t(res.ssgsea[1,-c(1,2)]))
colnames(tmp) <- 'res.ssgsea'
tmp$id <- rownames(tmp)
tmp <- tmp[order(tmp$res.ssgsea),]
#group pathents by ssGSEA score
num = 2
tmp$bin <- cut_number(tmp$res.ssgsea,n = num)
unique(tmp$bin)
tmp$type <- 'none'
tmp[tmp$bin==unique(tmp$bin)[1],]$type <- 'low'
tmp[tmp$bin==unique(tmp$bin)[num],]$type <- 'high'
table(tmp$type)

#ED. Fig. 9a: survival data
clin.data <- clin.luad
clin.data <- merge(clin.data, tmp,by.x='submitter_id',by.y='id')
clin.data <- clin.data[clin.data$type!='none',]
clin.data$days_to_event <- ifelse(clin.data$vital_status == "Dead",
                                  clin.data$days_to_death,
                                  clin.data$days_to_last_follow_up)
clin.data$event <- ifelse(clin.data$vital_status == "Dead", 1, 0)
surv_object <- Surv(time = clin.data$days_to_event, 
                    event = clin.data$event)
fit <- survfit(surv_object ~ clin.data$type,data = clin.data)
#plot
days <- 365*5
clin.data$n <- 1:nrow(clin.data)
ns <- na.omit(clin.data[clin.data$days_to_event<days,]$n)
plot(fit, col = c('red','blue'),  xlab = "Days", ylab = "Survival Probability",
     xlim = c(0,days), lwd=2.5)
legend("topright", legend = c(paste0(levels(factor(clin.data$type))[1],
                                     ', n = ',nrow(clin.data[clin.data$type==levels(factor(clin.data$type))[1],])),
                              paste0(levels(factor(clin.data$type))[2],
                                     ', n = ',nrow(clin.data[clin.data$type==levels(factor(clin.data$type))[2],]))),   
       col = c('red','blue'), lty = 1)
log_rank_test <- survdiff(surv_object[ns,] ~ clin.data[ns,]$type)
p_value <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
cat("Log-rank test p-value:", p_value)
text(x = 500, y = 0.2, labels = paste("p-value:", signif(p_value, 3)), col = "black")

#ED. Fig. 9b: correlation plot####
ls <- list(l1=c('VSN_T_N_TH','CD8'))
names <- ls[[1]]
tmp <- as.data.frame(t(res.ssgsea[res.ssgsea$Name %in% names,-c(1,2)]))
colnames(tmp) <- c('res.ssgsea1','res.ssgsea2')
ggplot(tmp,aes(res.ssgsea1,res.ssgsea2))+
  geom_point(color = "#6495ED",shape=1,size=1.5)+
  geom_smooth(method = "lm", se = T, color = "#1E90FF")+
  stat_cor(method = "pearson")+xlab(label = names[1])+ylab(label = names[2])+
  theme(panel.grid = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(),axis.title = element_text(),axis.text = element_text())

