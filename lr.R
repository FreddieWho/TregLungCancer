files <- list.files('../0713/','.tsv')


getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))


readlr <- function(file){
  tab <- read.table(paste0('../0713/',file),header = T,sep = '\t')
  ct <- strsplit(sub('.tsv','',file),'_') %>% unlist
  tab <- cbind(tab,str_split(tab$Ligand...Receptor,'---',simplify = T))
  colnames(tab)[3:4] <- c('lr1','lr2')
  
  tab$ct1 <- ct[1]
  tab$ct2 <- ct[2]
  
  tab
}

c(0.87684,5.9299e-06,0.87684,0.87684,0.87684,0.87684)

f1 <- readlr(files[3])
f1$ct1 <- 'CD8_C6-LAYN'
f1$ct2 <- 'CD4_C9-CTLA4_ADAM12-'
f1$lr <- paste(f1$lr1,f1$lr2,sep = '_')
f1$ct <- paste(f1$ct1,f1$ct2,sep = '___')

f2 <- readlr(files[4]);head(f2)
f2$ct1 <- 'CD8_C6-LAYN'
f2$ct2 <- 'CD4_C9-CTLA4_ADAM12+'
f2$lr <- paste(f2$lr1,f2$lr2,sep = '_')
f2$ct <- paste(f2$ct1,f2$ct2,sep = '___')

f3 <- readlr(files[5]);head(f3)
f3$ct1 <- 'CD4_C7-CXCL13'
f3$ct2 <- 'CD4_C9-CTLA4_ADAM12-'
f3$lr <- paste(f3$lr2,f3$lr1,sep = '_')
f3$ct <- paste(f3$ct1,f3$ct2,sep = '___')

f4 <- readlr(files[6]);head(f4)
f4$ct1 <- 'CD4_C7-CXCL13'
f4$ct2 <- 'CD4_C9-CTLA4_ADAM12+'
f4$lr <- paste(f4$lr2,f4$lr1,sep = '_')
f4$ct <- paste(f4$ct1,f4$ct2,sep = '___')

f5 <- readlr(files[1]);head(f5)
f5$ct1 <- 'CD8_C4-GZMK'
f5$ct2 <- 'CD4_C9-CTLA4_ADAM12-'
f5$lr <- paste(f5$lr1,f5$lr2,sep = '_')
f5$ct <- paste(f5$ct1,f5$ct2,sep = '___')

f6 <- readlr(files[2]);head(f6)
f6$ct1 <- 'CD8_C4-GZMK'
f6$ct2 <- 'CD4_C9-CTLA4_ADAM12+'
f6$lr <- paste(f6$lr1,f6$lr2,sep = '_')
f6$ct <- paste(f6$ct1,f6$ct2,sep = '___')

f7 <- readlr(files[7]);head(f7)
f7$ct1 <- 'CD8_C5-ZNF683'
f7$ct2 <- 'CD4_C9-CTLA4_ADAM12-'
f7$lr <- paste(f7$lr1,f7$lr2,sep = '_')
f7$ct <- paste(f7$ct1,f7$ct2,sep = '___')

f8 <- readlr(files[8]);head(f8)
f8$ct1 <- 'CD8_C5-ZNF683'
f8$ct2 <- 'CD4_C9-CTLA4_ADAM12+'
f8$lr <- paste(f8$lr1,f8$lr2,sep = '_')
f8$ct <- paste(f8$ct1,f8$ct2,sep = '___')

f <- dplyr::bind_rows(list(f1,f2,f3,f4,f5,f6,f7,f8))
f$z <- (f$contribution - mean(f$contribution))/sd(f$contribution)
f$q <- ifelse(f$ct == 'CD8_C6-LAYN___CD4_C9-CTLA4_ADAM12+',5.9299e-06,0.87684)
save(f,file = 'lr.rda')

pdf('lr.pdf',height = 5,width = 8)
f[grep('CCL5_SDC4|ADAM12$',f$lr),] %>% 
ggplot() +
  geom_point(aes(y = ct,x = lr,color = contribution,size = -log(q + 1e-3))) +
  scale_colour_gradientn(colors = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),name = 'Contribution',) +
  scale_size_continuous(range = c(4,7)) +
  theme_bw() + 
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 30,hjust = 1)) +
  xlab('') + ylab('')
dev.off()
