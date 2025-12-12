################################# ED.Fig.8l and m
#The intermediate files 'VSN_sympathetic.gct' and 'CD8.gct' are available in the TCGA folder

#ED.Fig.8l: Survival probabilities of VSN+sympathetic###
library(ggpubr)
library(ggplot2)
library(survival)
#survival plot####
res.ssgsea <- as.data.frame(fread('VSN_sympathetic.gct'))
#res.ssgsea is output from GenePattern(https://cloud.genepattern.org/gp/pages/index.jsf)
rownames(res.ssgsea) <- res.ssgsea[,1]
tmp <- data.frame(t(res.ssgsea[,-c(1,2)]))
tmp$id <- rownames(tmp)
clin.data <- clin.luad
clin.data <- merge(clin.data, tmp,by.x='submitter_id',by.y='id')

clin.data$VSN_group <- ifelse(clin.data$VSN > median(clin.data$VSN, na.rm = TRUE), "VSNHigh", "VSNLow")
clin.data$sympathetic_group <- ifelse(clin.data$sympathetic > median(clin.data$sympathetic, na.rm = TRUE), "sympatheticHigh", "sympatheticLow")
# Optional: combine two groupings
clin.data$combined_group <- paste(clin.data$VSN_group, clin.data$sympathetic_group, sep = "_")
table(clin.data$combined_group )
#group patients by ssGSEA score
clin.data$type <- clin.data$combined_group
unique(clin.data$type)
for (i in c("VSNLow_sympatheticLow")){#"VSNHigh_sympatheticLow","VSNLow_sympatheticLow","VSNLow_sympatheticHigh"
  dat_sub <- clin.data[clin.data$type %in% c(i,'VSNHigh_sympatheticHigh'),]
  dat_sub$days_to_event <- ifelse(dat_sub$vital_status == "Dead",
                                  dat_sub$days_to_death,
                                  dat_sub$days_to_last_follow_up)
  dat_sub$event <- ifelse(dat_sub$vital_status == "Dead", 1, 0)
  surv_object <- Surv(time = dat_sub$days_to_event, 
                      event = dat_sub$event)
  fit <- survfit(surv_object ~ dat_sub$type,data = dat_sub)
  
  #plot
  days <- 365*5
  dat_sub$n <- 1:nrow(dat_sub)
  ns <- na.omit(dat_sub[dat_sub$days_to_event<days,]$n)
  plot(fit, col = c('red','blue'),  xlab = "Days", ylab = "Survival Probability",
       xlim = c(0,days), lwd=2.5)
  legend("bottom", legend = c(paste0(levels(factor(dat_sub$type))[1],
                                     ', n = ',nrow(dat_sub[dat_sub$type==levels(factor(dat_sub$type))[1],])),
                              paste0(levels(factor(dat_sub$type))[2],
                                     ', n = ',nrow(dat_sub[dat_sub$type==levels(factor(dat_sub$type))[2],]))),   
         col = c('red','blue'), lty = 1)
  
  log_rank_test <- survdiff(surv_object[ns,] ~ dat_sub[ns,]$type)
  p_value <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
  cat("Log-rank test p-value:", p_value)
  text(x = 500, y = 0.3, labels = paste("p-value:", signif(p_value, 3)), col = "black")
}
#ED Fig.8m: Box plot of CD8 signature####
res.ssgsea <- fread('CD8.gct')
tmp <- data.frame(t(res.ssgsea[2,-c(1,2)]))
colnames(tmp) <- 'CD8 signature'
dat <- merge(dat_sub,tmp,by.x='submitter_id',by.y='row.names')
table(dat$type)
dat$type <- factor(dat$type,levels = c('VSNLow_sympatheticLow','VSNHigh_sympatheticHigh'))
cmp <- combn(levels(dat$type), 2, simplify = FALSE)
ggplot(dat,aes(type,`CD8 signature`,fill=type))+
  geom_boxplot(outliers = F)+
  stat_compare_means(
    comparisons = cmp,              # list of pairwise comparisons
    method      = "wilcox.test",    # same as Seurat's default
    label       = "p.format"
  )+guides(fill='none')+scale_fill_manual(values = c('blue','red'))
