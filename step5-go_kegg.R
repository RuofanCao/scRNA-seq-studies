## 
### ---------------
###
### Create: Jianming Zeng
### Date:  2023-01-16 20:53:22
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2023-01-16  First version 
###
### ---------------
rm(list=ls())
options(stringsAsFactors = F) 
getwd()
source('scRNA_scripts/lib.R') 
sce.all.int = readRDS('2-harmony_immune/sce.all_int.rds')
sce.all.int
load('phe_new.Rdata') 
phe = phe_new
colnames(phe)
sce.all.int=sce.all.int[,colnames(sce.all.int) %in% rownames(phe)]
identical(colnames(sce.all.int), rownames(phe))
sce.all.int@meta.data = phe
th =  theme(axis.text.x=element_text(angle=45,hjust = 1)) 
table(sce.all.int$celltype)
DimPlot(sce.all.int,group.by = 'celltype',
        repel = T,label = T)

#### 看各个亚群功能 ###### 
{
  sce = sce.all.int
  sel.clust = "celltype"
  sce <- SetIdent(sce, value = sel.clust)
  table(sce@active.ident) 
  markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                            min.pct = 0.25, 
                            thresh.use = 0.25)
  write.csv(markers,file='markers.csv')
  table(markers$cluster)
  library(dplyr)  
  library(clusterProfiler)  
  markers = markers[markers$p_val <0.01 & markers$avg_log2FC >1,]  
  table(markers$cluster)
  symbols_list = split(markers$gene,markers$cluster)
  source('scRNA_scripts/com_go_kegg_ReactomePA_human.R') 
  com_go_kegg_ReactomePA_human(symbols_list, pro='CD4T' )
}

###一直装不上org.Hs.eg.db这个包，修改一个镜像试试，几秒OK(这里选择的是Nanjing，678都行)
#ReactomePA四百多兆也是几秒的速度
#如果是install，可以使用chooseCRANmirror()
chooseBioCmirror()
#BiocManager::install('ReactomePA')
library('org.Hs.eg.db')
