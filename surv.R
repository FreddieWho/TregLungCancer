setwd('~/tmp/lung/0714/')
library(dplyr)
library(survival)
library(survminer)

cli <- read.table('~/iDB/bulkRNA/xena_pan/Survival_SupplementalTable_S1_20171025_xena_sp',sep= '\t',header = T)
exp <- data.table::fread('~/iDB/bulkRNA/xena_pan/tcga_RSEM_gene_tpm',sep= '\t',header = T)
id <- read.table('~/iDB/bulkRNA/xena_pan/probeMap%2Fgencode.v23.annotation.gene.probemap',sep= '\t',header = T)
colnames(id)[1] <- 'sample'
id <- id[,1:2]

surv <- function(gene1 = c('FOXP3','CTLA4'),
                 gene2 = c('CD3E','CD3G','CD3D'),
                 lab1 = 'ADAM12-CTLA4+FOXP3',
                 lab2 = 'ADAM12+CTLA4+FOXP3'){
  meta <- subset(cli,cancer.type.abbreviation %in% c('LUAD'))
  rownames(meta) <- meta$sample
  expr1 <- subset(id,gene %in% gene1) %>% left_join(.,exp) %>% dplyr::select(-sample,-gene)
  expr1 <- expr1[,meta$sample[meta$sample %in% colnames(expr1)]] %>% colSums()
  expr2 <- subset(id,gene %in% gene2) %>% left_join(.,exp) %>% dplyr::select(-sample,-gene)
  expr2 <- expr2[,meta$sample[meta$sample %in% colnames(expr2)]] %>% colSums()
  
  expr <- expr1-expr2
  meta <- meta[names(expr),]
  meta$expr <- expr
  
  df_plot <- meta
  df_plot$geneStatus <- cut(df_plot$expr,c(-Inf,quantile(df_plot$expr,c(0.25,0.75)),Inf),c(lab1,'skip',lab2))
  df_plot <- subset(df_plot, geneStatus != 'skip')
  
  df_plot
  # survfit(Surv(time = df_plot$OS.time,event = df_plot$OS)~geneStatus, data=subset(df_plot, geneStatus != 'skip'))
}

dffit1 <- surv(c('FOXP3','CTLA4'),c('CD3E','CD3G','CD3D'),'CTLA4+FOXP3+ low','CTLA4+FOXP3 hi')
dffit2 <- surv(c('FOXP3','CTLA4','ADAM12'),c('CD3E','CD3G','CD3D'),'ADAM12+CTLA4+FOXP3 low','ADAM12+CTLA4+FOXP3 hi')

print(ggsurvplot(survfit(Surv(time = dffit1$OS.time,
                              event = dffit1$OS)~geneStatus, 
                         data=subset(dffit1, geneStatus != 'skip')),
                 pval = T,palette = 'aaas',conf.int = F))
print(ggsurvplot(survfit(Surv(time = dffit2$OS.time,
                              event = dffit2$OS)~geneStatus, 
                         data=subset(dffit2, geneStatus != 'skip')),
                 pval = T,palette = 'aaas',conf.int = F))


pdf('surv.pdf',height = 5,width = 7)
print(ggsurvplot(survfit(Surv(time = dffit1$OS.time,
                              event = dffit1$OS)~geneStatus, 
                         data=subset(dffit1, geneStatus != 'skip')),
                 pval = T,palette = 'aaas',conf.int = F))
print(ggsurvplot(survfit(Surv(time = dffit2$OS.time,
                              event = dffit2$OS)~geneStatus, 
                         data=subset(dffit2, geneStatus != 'skip')),
                 pval = T,palette = 'aaas',conf.int = F))
dev.off()


