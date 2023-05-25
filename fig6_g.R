setwd('~/tmp/lung/code/')
library(dplyr)
library(survival)
library(survminer)

cli <- read.table('~/iDB/bulkRNA/xena_pan/Survival_SupplementalTable_S1_20171025_xena_sp',sep= '\t',header = T)
exp <- data.table::fread('~/iDB/bulkRNA/xena_pan/tcga_RSEM_gene_tpm',sep= '\t',header = T)
id <- read.table('~/iDB/bulkRNA/xena_pan/probeMap%2Fgencode.v23.annotation.gene.probemap',sep= '\t',header = T)
colnames(id)[1] <- 'sample'
id <- id[,1:2]

meta <- subset(cli,cancer.type.abbreviation %in% c('LUAD'))
rownames(meta) <- meta$sample
expr1 <- subset(id,gene %in% c('FOXP3','CTLA4','ADAM12')) %>% left_join(.,exp) %>% select(-sample,-gene)
expr1 <- expr1[,meta$sample[meta$sample %in% colnames(expr1)]] %>% colSums()
expr2 <- subset(id,gene %in% c('CD3E','CD3G','CD3D')) %>% left_join(.,exp) %>% select(-sample,-gene)
expr2 <- expr2[,meta$sample[meta$sample %in% colnames(expr2)]] %>% colSums()

expr <- expr1-expr2
meta <- meta[names(expr),]
meta$expr <- expr

df_plot <- meta
df_plot$geneStatus <- cut(df_plot$expr,c(-Inf,quantile(df_plot$expr,c(0.25,0.75)),Inf),c('ADAM12-CTLA4+ Treg','skip','ADAM12+CTLA4+ Treg'))
df_plot <- subset(df_plot, geneStatus != 'skip')

fit <- survfit(Surv(time = df_plot$OS.time,event = df_plot$OS)~geneStatus, data=subset(df_plot, geneStatus != 'skip'))
print(ggsurvplot(fit,pval = T,palette = 'aaas',conf.int = F))

pdf('fig6_g.pdf',height = 5,width = 7)
ggsurvplot(fit,pval = T,palette = 'aaas')
dev.off()
