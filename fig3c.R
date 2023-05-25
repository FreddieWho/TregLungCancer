setwd('~/tmp/lung/code/')

library(ggpubr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
library(stringr)


getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))

load('filtered_seurat.rda')
load('raw_seurat.rda')


f1_c <- function(obj,gene,gender,tissue){
  data <- data.frame(expr = subset(obj,Gender == gender & Tissue == tissue)@assays$RNA@data[gene,],
                     subset(obj,Gender == gender & Tissue == tissue)@meta.data)
  
  ggplot(data) +
    geom_violin(aes(x = majorCluster,y = expr,fill = Smoker,color = Smoker),
                alpha = 0.6) +
    stat_compare_means(aes(x = majorCluster,y = expr,group = Smoker),
                       label = 'p.signif',vjust = 1) + 
    geom_point(aes(x = majorCluster,y = expr,color = Smoker),
                 position=position_jitterdodge()) +
    ylab(paste0('Expression\nGender: ',gender,'; Tissue: ',tissue)) + xlab('') +  ggtitle(gene) +
    scale_color_manual(labels = c('Non-smoker','Smoker'),name = 'Smoking history',values = c('#00468B','#ED0000')) +
    guides(fill = F) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 15),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 15)) + coord_flip()
}

pdf('fig3c_vln.pdf',height = 6,width = 6)
f1_c(obj_filtered,'ADAM12','M','Tumor')
f1_c(obj_filtered,'FANK1','M','Tumor')
f1_c(obj_filtered,'ACTG2','M','Tumor')
dev.off()




