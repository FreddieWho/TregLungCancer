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


load('../code/filtered_seurat.rda')
load('../code/raw_seurat.rda')

treg <- subset(obj_filtered,majorCluster == 'CD4_C9-CTLA4')
treg$cutoff_0 <- 'CD4_C9-CTLA4_ADAM12-'
treg$cutoff_0[which(log(treg@assays$RNA@data['ADAM12',]) > -3)] <- 'CD4_C9-CTLA4_ADAM12+'
treg$subC <- treg$cutoff_0


monocle <- subset(obj_filtered,majorCluster %in% grep('^CD4_',unique(obj_filtered$majorCluster),value = T))
monocle$subC <- monocle$majorCluster
monocle@meta.data[rownames(treg@meta.data),'subC'] <- treg$subC


library(monocle)
cds <- newCellDataSet(monocle@assays$RNA@counts,
                       phenoData = monocle@meta.data %>% new("AnnotatedDataFrame",data = .),
                       featureData = data.frame(gene_short_name = rownames(monocle),
                                                row.names = rownames(monocle)) %>% 
                         new("AnnotatedDataFrame",data = .),
                       expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~subC',cores = 30)

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]

cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

cds <- reduceDimension(cds, method = 'DDRTree')

cds <- orderCells(cds)

cds <- orderCells(cds, root_state = GM_state(cds))

pdf('f2_traj_subC.pdf',height = 6,width = 11)
plot_cell_trajectory(cds,color_by = 'subC',cell_size = 2,show_branch_points = F) +
  theme(legend.position = 'top',legend.text = element_text(size = 15)) +
  scale_color_d3(name = '',alpha = 0.5) +
plot_cell_trajectory(cds,color_by = 'Pseudotime')
dev.off()

plot_cell_trajectory(cds,markers = c('TNFRSF9','ADAM12','SDC4'),cell_size = 2,use_color_gradient = T) 
save(cds,file = 'monocle_allCD4.rda')
