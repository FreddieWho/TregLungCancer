setwd('~/tmp/lung/0715/')

library(ggpubr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))


go <- read.csv('metascape.csv')
go <- go[grep('Summary',go$GroupID),]

go$GeneRatio <- sapply(strsplit(go$InTerm_InList,'/'),function(x) as.numeric(x[1]))/sapply(strsplit(go$InTerm_InList,'/'),function(x) as.numeric(x[2]))

pdf('f2_go.pdf',height = 8,width = 10)
ggplot(go %>% arrange(LogP) %>% mutate(Description = factor(Description,levels = rev(Description)))) +
  geom_bar(aes(x = -LogP,y = Description,fill = GeneRatio),stat = 'identity') +
  scale_fill_gradientn(colours = rev(getPalette(10)[1:5]),name = 'GeneRatio') +
  theme_bw() + ylab('') + xlab('-log10(p)') +
  theme(axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))
dev.off()
