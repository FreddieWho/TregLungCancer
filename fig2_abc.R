setwd('~/tmp/lung/code/')

library(ggpubr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))


load('~/tmp/lung/code/filtered_seurat.rda')
load('~/tmp/lung/code/raw_seurat.rda')


treg <- subset(obj_filtered,majorCluster == 'CD4_C9-CTLA4')

pdf('fig2_a.pdf',width = 5,height = 5)
ggplot(data.frame(adam12 = treg@assays$RNA@data['ADAM12',])) +
  geom_histogram(aes(x = adam12)) +
  theme_bw() +
  xlab('ADAM12') +
  theme(axis.title = element_text(size = 20))
ggplot(data.frame(adam12 = treg@assays$RNA@data['ADAM12',] %>% log())) +
  geom_histogram(aes(x = adam12)) +
  theme_bw() +
  geom_vline(xintercept = -3,linetype = 'dashed') +
  xlab('log(ADAM12)') +
  theme(axis.title = element_text(size = 20))
dev.off()


treg$cutoff_0 <- 'CD4_C9-CTLA4_ADAM12-'
treg$cutoff_0[which(log(treg@assays$RNA@data['ADAM12',]) > -3)] <- 'CD4_C9-CTLA4_ADAM12+'
treg$subC <- treg$cutoff_0

pdf('vln_ADAM12+_-CTLA4_Treg_cutoff.pdf',height = 5,width = 5)
VlnPlot(treg,'ADAM12',group.by = 'cutoff_0')
dev.off()

Idents(treg) <- treg$cutoff_0
mk <- FindAllMarkers(treg,test.use = 'MAST',)
mk <- subset(mk, p_val_adj < 0.05)

write.csv(mk,file = 'mk_subCTLA4+_ADAM12Cut.csv')

obj_filtered$subC <- obj_filtered$majorCluster
obj_filtered@meta.data[rownames(treg@meta.data),'subC'] <- treg$cutoff_0


plotG <- c('ADAM12','SDC4','GCNT1','PIKFYVE','SLC16A1','IL1R1','CD177','ACSL1','ACSL4','HACD1','LPIN1','SPTLC2','SACM1L','BDH2','GM2A','INSIG2','IL17RB')
plotG <- plotG[plotG %in% mk$gene]

mk <- read.csv('mk_subCTLA4+_ADAM12Cut.csv')
adam12_sig <- subset(mk,p_val_adj < 0.05 & cluster == 'CD4_C9-CTLA4_ADAM12-')$gene
go <- enrichGO(adam12_sig,org.Hs.eg.db,keyType = 'SYMBOL')@result
go_sig <- subset(go,p.adjust < 0.05)
go_sig$Description <- factor(go_sig$Description, levels = rev(go_sig$Description))
pdf('f2_go.pdf',height = 6,width = 10)
ggplot(go_sig) +
  geom_point(aes(x = -log10(p.adjust),y = Description,size = Count/45,color = -log10(p.adjust))) +
  scale_colour_gradientn(colours = rev(getPalette(10)[1:4]),name = 'GeneRatio') +
  theme_bw() + ylab('') +
  theme(axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))
dev.off()

fig4_dot <- function(mca,features,color){
  mca$subC <- mca[['subC']]
  if(length(features) > 1){
    exp <- mca@assays$RNA@data[features,] %>% as.data.frame() %>% t() %>% as.data.frame()
  } else {
    exp <- mca@assays$RNA@data[features,] %>% as.data.frame()
    colnames(exp) <- features
  }
  
  exp$subC <- mca@meta.data[rownames(exp),'subC']
  
  trunc_z <- function(x) MinMax((x-mean(x))/sd(x),-3,3)
  iscale <- function(x) (x-min(x))/(max(x) - min(x))
  
  exp %>% melt() %>% group_by(subC,variable) %>% 
    summarise(mean = mean(value),frac = sum(value > 0)/length(value)) %>% 
    group_by(variable) %>% 
    mutate(mean_z = trunc_z(mean),
           mean_z_scale = iscale(trunc_z(mean)),
           frac_scale = iscale(frac)) -> plot_df
  
  print(unique(plot_df$subC))
  
  plot_df$subC <- factor(plot_df$subC,
                         levels = c('CD4_C9-CTLA4_ADAM12+',
                                    'CD4_C9-CTLA4_ADAM12-',
                                    grep('CTLA4',unique(plot_df$subC),invert = T,value = T) %>% sort()))
  ggplot(plot_df) +
    geom_point(aes(x = variable,y = subC, color = mean_z,size = frac)) +
    scale_colour_gradientn(colours = color,name = 'Expression') +
    scale_y_discrete(limits=rev(levels(plot_df$subC))) +
    scale_size(name = 'Proportion') +
    # facet_grid(~gType,space = 'free_x',scale = 'free_x') +
    theme_classic() + xlab('') + ylab('') +
    theme(axis.text.x = element_text(angle = 20,hjust = 1,size = 10),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_blank())
}
pdf('fig4_a.pdf',height = 5,width = 10)
fig4_dot(obj_filtered,c(plotG),color = rev(getPalette(10)))
dev.off()
residentSig <- c('ITGAE','STIM1')
# residentSig <- c('ITGAE','STIM1','ZNF683','STIM2')
pdf('dot_residentSig.pdf',height = 3,width = 6)
fig4_dot(treg,residentSig,color = rev(getPalette(10)))
dev.off()
