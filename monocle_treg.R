setwd('~/tmp/lung/0712/')

library(ggpubr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)

getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))


load('~/tmp/lung/code/filtered_seurat.rda')
load('~/tmp/lung/code/raw_seurat.rda')

treg <- subset(obj_filtered,majorCluster == 'CD4_C9-CTLA4')
treg$cutoff_0 <- 'CD4_C9-CTLA4_ADAM12-'
treg$cutoff_0[which(log(treg@assays$RNA@data['ADAM12',]) > -3)] <- 'CD4_C9-CTLA4_ADAM12+'
treg$subC <- treg$cutoff_0




library(monocle)
cds <- newCellDataSet(treg@assays$RNA@counts,
                      phenoData = treg@meta.data %>% new("AnnotatedDataFrame",data = .),
                      featureData = data.frame(gene_short_name = rownames(treg),
                                               row.names = rownames(treg)) %>% 
                        new("AnnotatedDataFrame",data = .),
                      expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~subC',cores = 30)

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]

cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

cds <- reduceDimension(cds, method = 'DDRTree')

cds <- orderCells(cds)

save(cds,file = 'monocle_treg.rda')
load('~/tmp/lung/0712/monocle_treg.rda')

pdf('f2_traj_treg.pdf',height = 6,width = 9)
plot_cell_trajectory(cds,color_by = 'subC',cell_size = 2,show_branch_points = F) +
  theme(legend.position = 'right',legend.text = element_text(size = 15)) +
  scale_color_d3(name = '',alpha = 0.5)
dev.off()

p1 <- plot_cell_trajectory(cds,markers = c('ADAM12'),cell_size = 2,use_color_gradient = T)+scale_color_viridis_c(name = 'Expression')
p2 <- plot_cell_trajectory(cds,markers = c('SDC4'),cell_size = 2,use_color_gradient = T)+scale_color_viridis_c(name = 'Expression')
p3 <- plot_cell_trajectory(cds,markers = c('TNFRSF9'),cell_size = 2,use_color_gradient = T) +scale_color_viridis_c(name = 'Expression')
p4 <- plot_cell_trajectory(cds,markers = c('ITGAE'),cell_size = 2,use_color_gradient = T) +scale_color_viridis_c(name = 'Expression')

ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="right",font.label = 20)
dev.off()

pdf('sf2_traj_fank1.pdf',height = 4,width = 4)
plot_cell_trajectory(cds,markers = c('FANK1'),
                     cell_size = 2,use_color_gradient = T)+
  scale_color_viridis_c(name = 'Expression')
dev.off()
