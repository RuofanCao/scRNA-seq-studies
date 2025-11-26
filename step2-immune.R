#文中只取出了免疫细胞
sce = sce.all.int
set.seed(123456789)
phe = sce.all.int@meta.data
table(phe$RNA_snn_res.0.1) 
choose_cellIds=rownames(phe)[phe$RNA_snn_res.0.1 %in% c('0','1','2','4','5','6','7','8') ]

sce.all=sce.all.int[,choose_cellIds]
sce.all=CreateSeuratObject(
  counts = sce.all@assays$RNA$counts,
  meta.data = sce.all@meta.data
) 
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident) 

# 如果发现第一次降维聚类分群里面的干扰群
# 就首先去除干扰，然后把降维聚类分群重新跑一次
if(F){
  table(sce.all.int$RNA_snn_res.0.1)
  kp = !sce.all.int$RNA_snn_res.0.1 %in% c(6) 
  table(kp)
  trash_sce = sce.all.int[,!kp]
  save(trash_sce,file = 'trash_sce.Rdata')
  
  sce.all=sce.all.int[,kp]
  sce.all=CreateSeuratObject(
    counts = sce.all@assays$RNA$counts,
    meta.data = sce.all@meta.data
  )
}
sce.all
sp='human' 
###### step3: harmony整合多个单细胞样品 ######

if(T){
  dir.create("2-harmony_immune")
  getwd()
  setwd("2-harmony_immune")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all)
  setwd('../')
  
}

###### step4:  看标记基因库 ######
# 原则上分辨率是需要自己肉眼判断，取决于个人经验
# 为了省力，我们直接看  0.1和0.8即可

table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.3) 

###### step5: 人工手动命名 ###### 
sp='human'
tmp=sce.all.int@meta.data
colnames(sce.all.int@meta.data) 
table(sce.all.int$RNA_snn_res.0.3)
Idents(sce.all.int) = sce.all.int$RNA_snn_res.0.3

###查看默认Idents
Idents(sce.all.int)

# 付费环节 800 元人民币
# 如果是手动给各个单细胞亚群命名
###### 常见分群
# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# macrophages (CD68, CD163),
# 'CCL3L1' ,   #M2
# 'FABP4',  #M1
# Monocytes  (CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
# Photoreceptor cells (RCVRN), 
# Fibroblasts (FGF7, MME), 
# Neutrophil ('G0S2', 'S100A9','S100A8','CXCL8')
# Endothelial cells (PECAM1, VWF). 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
# immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 


#####注释亚群
#通常我们第一层次降维聚类分群：immune (CD45+,PTPRC),epithelial/cancer (EpCAM+,EPCAM),stromal (CD10+,MME,fibro or CD31+,PECAM1,endo)
#文章中一般第一次分群为常见的细胞群：上皮细胞（EPCAM、KRT19、CLDN4）、基质（PECAM1、CLO1A2、VWF）、增殖性（MKI67、STMN1、PCNA）、T（CD3D、CD3E、CD2）、B（CD79A，IGHG1，MS4A1），NK（KLRD1、GNLY、KLRF1）和髓系（CSF1R、CSF3R、CD68）细胞。

#绝大部分文章都是抓住免疫细胞亚群进行细分，包括淋巴系（T,B,NK细胞）和髓系（单核，树突，巨噬，粒细胞）的两大类作为第二次细分亚群。但是也有不少文章是抓住stromal 里面的 fibro 和endo进行细分
genes_to_check <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD8B",'CD4', ##T
                    "NKG7", "GNLY", "FCGR3", "ANCAM1", "CD160", #NK
                    "CD19", "MS4A1", "CD79B", "CD79A", #B
                    "MZB1", "XBP1", "JCHAIN", #Plasma
                    "TPSAB1", "TPSB2",'CPA3','CD38','KIT','NR4A2','NR4A1','TNFSF10', #Mast
                    "LYZ", "S100A8", "S100A9", "CD14", "CD68", #Myeloid
                    "PPBP", "MPIG6B", #Platelets
                    "CXCL8", "DCN", "SLPI", "CXADR",'FAP','CFD','KRT8', #Fibroblasts
                    "PECAM1", "VWF", "ENG", "MCAM", #Endothelial
                    "EPCAM", "SOX9", "MET", "EGFR", "MSLN", "YAP1", #Malignant
                    "MKI67", "BIRC5", "CDK1",'TOP2A') #Proliferative

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
rownames(sce.all.int)
length(rownames(sce.all.int))
p = DotPlot(sce.all.int, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p

###重要：注意人工肉眼对应，需要根据不同数据集修改分群所对应的细胞亚群
if(T){
  sce.all.int
  celltype=data.frame(ClusterID=0:15 ,
                      celltype=0:15 ) 
  #定义细胞亚群        
  celltype[celltype$ClusterID %in% c( 4,10,11,14 ),2]='B'
  celltype[celltype$ClusterID %in% c( 1,5,6,8 ),2]='CD4T'
  celltype[celltype$ClusterID %in% c( 0,2,13 ),2]='CD8T'
  celltype[celltype$ClusterID %in% c( 12 ),2]='Plasma'
  celltype[celltype$ClusterID %in% c( 7,9 ),2]='NK'
  celltype[celltype$ClusterID %in% c(3 ),2]='Myeloid'
  celltype[celltype$ClusterID %in% c( 15 ),2]='Proliferative'

  
  head(celltype)
  celltype
  table(celltype$celltype)
  sce.all.int@meta.data$celltype = "NA"
  
  ###上面选择的是RNA_snn_res.0.3，如果是别的resolution，是需要修改的
  for(i in 1:nrow(celltype)){
    sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  table(sce.all.int$celltype)
  Idents(sce.all.int)=sce.all.int$celltype
  
  
}

###查看亚群注释后marker表达是否对应的上，若不合理需重新检查上一步细胞亚群注释代码
p = DotPlot(sce.all.int, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p

###去除未被注释上的细胞亚群
if("celltype" %in% colnames(sce.all.int@meta.data ) ){
  
  sel.clust = "celltype"
  sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
  table(sce.all.int@active.ident) 
  
  phe=sce.all.int@meta.data
  save(phe,file = 'phe.Rdata')
  pdf('celltype-vs-orig.ident.pdf',width = 10)
  gplots::balloonplot(table(sce.all.int$celltype,sce.all.int$orig.ident))
  dev.off()
  dir.create('check-by-celltype')
  setwd('check-by-celltype')
  source('../scRNA_scripts/check-all-markers.R')
  setwd('../') 
  getwd()
}

library(paletteer)
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
umap1=DimPlot(sce.all.int, reduction = "umap",
              # split.by = "orig.ident",
              group.by = "RNA_snn_res.0.3",label = T,cols = color)
umap1


umap2=DimPlot(sce.all.int, reduction = "umap",
              # split.by = "orig.ident",
              group.by = "celltype",label = T,repel = T,cols = color)

umap2
umap1+umap2
