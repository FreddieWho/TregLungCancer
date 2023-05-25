tab <- read.table('~/tmp/lung/0709/GSEA_overlap_adam12-layn.tsv',sep = '\t',skip = 10,header = F,nrows = 10)
tab$V1 <- factor(tab$V1,levels = tab %>% arrange(desc(V7)) %>% pull(V1))

getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))

pdf('dot_gseaOverlap.pdf',height = 6,width = 10)
ggplot(tab) +
  geom_point(aes(x = -log10(V7),y = V1,color = V5),size = 3) +
  ylab('') + xlab('-log10(fdr)') +
  scale_colour_gradientn(colours = rev(getPalette(10)[1:4]),name = 'GeneRatio') +
  theme_bw() +
  theme(text= element_text(size = 15))
dev.off()
