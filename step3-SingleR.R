######细胞注释-使用SingleR
#BiocManager::install('SingleR')
library(SingleR)
library(celldex)
library(dplyr)
library(stringr)
library(Seurat)
library(pheatmap)
library(ReactomeGSA)
library(ggplot2)
library(singleseqgset)
library(devtools)
#install_github("arc85/singleseqgset")
#BiocManager::install ('ReactomeGSA')
# singleR注释
#hpca.se <- HumanPrimaryCellAtlasData()
#save(hpca.se,file = 'hpca.RData')
load('hpca.RData')
load('bpe.RData')
unique(hpca.se$label.main)
unique(hpca.se$label.fine)
unique(bpe.se$label.main)
unique(bpe.se$label.fine)
#bpe.se <- BlueprintEncodeData()
#save(bpe.se,file = 'bpe.RData')
str(sce)
sce = readRDS('2-harmony/sce.all_int.rds')
anno <- SingleR(sce@assays$RNA$data,
                ref = list(BP=bpe.se,HPCA=hpca.se),
                labels = list(bpe.se$label.fine,hpca.se$label.main),
                clusters = sce@meta.data$seurat_clusters
)

plotScoreHeatmap(anno,clusters = anno@rownames,show_colnames = T)
table(anno$labels)

celltype = data.frame(ClusterID=rownames(anno), 
                      celltype=anno$labels, 
                      stringsAsFactors = F) 

sce@meta.data$singleR = celltype[match(sce@meta.data$seurat_clusters,celltype$ClusterID),'celltype']
table(sce$singleR)
source('scRNA_scripts/mycolors.R')
p1 = DimPlot(sce, reduction = "tsne",
             group.by = "singleR",label = T,cols = mycolors,label.box=T) 
p1
p+p1
ggsave('T_combine.pdf',width = 15,height = 7)


####自动注释也需要匹配人工验证

Idents(sce) = sce$singleR
VlnPlot(sce,
        features = c("CD3D","CD3E","CD4","CD8A",'MS4A1','IGHG1','MZB1', 'CD79A'),
        pt.size = 0,
        ncol = 4,
        cols=mycolors)

saveRDS(sce, "T_sce_celltype.rds")

setwd('../')
