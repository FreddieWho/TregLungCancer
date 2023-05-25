library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggrepel)

# Loading raw data #############################################################
raw_meta <- read.table('../data/GSE99254.txt',header = T,stringsAsFactors = F)
patientInfo <- read.csv('../data/GSE99254_patientInfo.csv',stringsAsFactors = F)
raw_meta <- left_join(raw_meta,patientInfo)
rownames(raw_meta) <- raw_meta$UniqueCell_ID
raw_meta <- raw_meta[,-1]

raw_cts <- data.table::fread('../data/GSE99254_NSCLC.TCell.S12346.count.txt',
                             header = T,stringsAsFactors = F) %>% 
  as.data.frame()
## remove duplicated gene and na
raw_cts <- raw_cts[!duplicated(raw_cts$symbol),]
raw_cts <- raw_cts[!is.na(raw_cts$symbol),]
## set rownames
rownames(raw_cts) <- raw_cts$symbol

raw_cts <- raw_cts[,-c(1:2)]

# Creating Seurat object #######################################################
obj <- CreateSeuratObject(raw_cts,meta.data = raw_meta)
obj$Tissue <- case_when(substr(obj$sampleType,0,1) == 'N' ~ 'Adjacent',
                        substr(obj$sampleType,0,1) == 'T' ~ 'Tumor',
                        substr(obj$sampleType,0,1) == 'P' ~ 'Blood',)
obj$Gender <- sub('.*/','',obj$Age.Sex)
obj$Age <- as.numeric(sub('/[FM]','',obj$Age.Sex))
obj$RP_PCT <- PercentageFeatureSet(obj,'^RP[LS]')
save(obj,file = 'raw_seurat.rda')

# Repeat original results from the literature ##################################
obj_filtered <- subset(obj,majorCluster %in% grep('CD[48]_C',unique(obj$majorCluster),value = T))
obj_filtered <- NormalizeData(obj_filtered)
save(obj_filtered,file = 'filtered_seurat.rda')
