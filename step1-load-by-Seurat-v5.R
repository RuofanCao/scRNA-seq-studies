###
### Create: Jianming Zeng
### Date:  2023-12-31  
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2023-12-31   First version 
### 

rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
getwd()

####单个txt/csv/tsv数据读取###
library(data.table)
ct1=fread("GSE72056/GSE72056_melanoma_single_cell_revised_v2.txt.gz",data.table = F)
ct1[1:5,1:5]
ct1[,1]
ct1 = ct1[!duplicated(ct1$Cell),]
colnames(ct1)[1:5]
rownames(ct1) <- ct1$Cell
ct1 <- ct1[,-1]
head(ct1)[1:5,1:2]
anno1 = ct1[1:3,]
ct1 =ct1[-c(1:3),]
sce1=CreateSeuratObject(counts =  ct1 ,
                       min.cells = 5,
                       min.features = 300,)
dim(sce1)
sce1

##orig.ident(样品名),nCount_RNA(每个细胞中所有基因的测序读段总数。衡量每个细胞测序数据总量的指标)
#nFeature_RNA（每个细胞中检测到的基因数量）
sce1@meta.data$study = 'GSE72056'
Idents(sce1) = rownames(sce1)

sp='human'
# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
#目的是为了去除质量较差细胞，低质量细胞会形成独特的亚群，使分群结果变得复杂
dir.create("GSE72056/1-QC")
setwd("GSE72056/1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../../scRNA_scripts/qc.R')
sce.all = sce1
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()

###### step3: harmony整合多个单细胞样品 ######
#细胞筛选之后，需要进行harmony整合。
#前面在创建总的Seurat对象时，使用了merge函数对多个Seurat进行了简单的合并。merge只是按照行和列进行了合并，并未对数据进行其他处理。
#在拿到下游单细胞矩阵前，样本经历了多个实验环节的处理，时间、处理人员、试剂以及技术平台等变量都会成为混杂因素，可能会导致数据产生批次效应（batch effect）。
#目前单细胞测序常用的去批次算法有Harmony，CCA，RPCA,FastMNN,scVI等。这里采用Harmony进行演示。

if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../../') 
}



####第二个数据集
library(data.table)
ct2=fread("GSE115978/GSE115978_tpm.csv.gz",data.table = F)
ct2[1:5,1:5]
ct2[,1]
anno2 = fread("GSE115978/GSE115978_cell.annotations.csv.gz",data.table = F)
ct2 = ct2[!duplicated(ct2$V1),]
colnames(ct2)[1:5]
rownames(ct2) <- ct2$V1
ct2 <- ct2[,-1]
ct2[1:5,1:2]
sce2=CreateSeuratObject(counts =  ct2 ,
                        min.cells = 5,
                        min.features = 300,)
dim(sce2)
sce2

sce2@meta.data$study = 'GSE115978'

sp='human'
###### step2: QC质控 ######
dir.create("GSE115978/1-QC")
setwd("GSE115978/1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../../scRNA_scripts/qc.R')
sce.all = sce2
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()

###### step3: harmony整合多个单细胞样品 ######
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../../') 
}



####第三个数据集
anno3 = fread("GSE120575/GSE120575_patient_ID_single_cells.txt.gz",data.table = F)
anno3 = anno3[20:16311,]
colnames(anno3) <- anno3[1,]
anno3 <- anno3[-1,]
anno3 <- anno3[, !apply(is.na(anno3), 2, all)]
anno3 = na.omit(anno3)
ct3=read.csv("GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",  sep = "\t",quote = "")
ct3[1:5,1:5]
#可以明显看到第一列是错位的
tail(ct3)
ct3[1:5,16285:16292]
names = colnames(ct3)[-1]
names[988]
identical(names,anno3$title)
table(names == anno3$title)
# 发现不同的是符号
setdiff_names <- setdiff(names, anno3$title)
setdiff_anno3 <- setdiff(anno3$title, names)

names <- gsub("\\.", "-", names)
table(names == anno3$title)
ct3 <- ct3[-1,]
ct3[1:5,15700:15900]
colnames(ct3) <- names
ct3[1:5,1:2]
#删除有NA的行 ct3 = na.omit(ct3)
ct3 <- ct3[, !apply(is.na(ct3), 2, all)]
ct3 = na.omit(ct3)
identical(colnames(ct3),anno3$title)
#调整行名顺序与列名完全一致
p = identical(colnames(ct3),anno3$title);p
if(!p) anno3 = anno3[match(colnames(ct3),anno3$title),]


sce3=CreateSeuratObject(counts =  ct3 ,
                        min.cells = 5,
                        min.features = 300,)
dim(sce3)
names(sce3@assays$RNA@layers)
sce3[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce3, assay = "RNA", layer = "counts")
sce3

sce3@meta.data$study = 'GSE120575'
phe3 = sce3@meta.data
sp='human'
###### step2: QC质控 ######
dir.create("GSE120575/1-QC")
setwd("GSE120575/1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../../scRNA_scripts/qc.R')
sce.all = sce3
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()

###### step3: harmony整合多个单细胞样品 ######
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../../scRNA_scripts/harmony.R')
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../../') 
}




# 合并三个数据集：
GSE115978 = readRDS('GSE115978/2-harmony/sce.all_int.rds') 
table(GSE115978$orig.ident)
GSE120575 = readRDS('GSE120575/2-harmony/sce.all_int.rds')
table(GSE120575$orig.ident)
GSE72056 = readRDS('GSE72056/2-harmony/sce.all_int.rds')
table(GSE72056$orig.ident)
sceList = list(
  GSE115978 = CreateSeuratObject(
    counts = GSE115978@assays$RNA$counts
  ), 
  GSE120575 = CreateSeuratObject(
    counts = GSE120575@assays$RNA$counts
  ),
  GSE72056 = CreateSeuratObject(
    counts = GSE72056@assays$RNA$counts
  )
)
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = c('GSE115978','GSE120575','GSE72056')  ) 
names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )


#这三个数据集都有不同数量的样品，合并后就是一个大的对象了：
table(GSE115978$orig.ident)
table(GSE120575$orig.ident)
table(GSE72056$orig.ident)

table(sce.all$orig.ident) 


sp='human'
# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
#目的是为了去除质量较差细胞，低质量细胞会形成独特的亚群，使分群结果变得复杂
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()


# ###文章中的过滤条件，其他数据可以不运行
# sce.all.filt
# tmp=subset(sce.all.filt,subset = percent_mito < 20 &
#              # nCount_RNA <10000 &
#               nFeature_RNA <5000 &
#              nFeature_RNA>200)
# tmp
# sce.all=tmp
# sort(table(tmp$orig.ident))
# sce.all.filt =tmp


###### step3: harmony整合多个单细胞样品 ######
#细胞筛选之后，需要进行harmony整合。
#前面在创建总的Seurat对象时，使用了merge函数对多个Seurat进行了简单的合并。merge只是按照行和列进行了合并，并未对数据进行其他处理。
#在拿到下游单细胞矩阵前，样本经历了多个实验环节的处理，时间、处理人员、试剂以及技术平台等变量都会成为混杂因素，可能会导致数据产生批次效应（batch effect）。
#目前单细胞测序常用的去批次算法有Harmony，CCA，RPCA,FastMNN,scVI等。这里采用Harmony进行演示。

if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../') 
}
# 
# sce.all.int <- RunUMAP(sce.all.int,  dims = 1:15, 
#                      reduction = "harmony")
# seuratObj <- RunUMAP(sce.all.int,  dims = 1:15 )
# colnames(sce.all.int@meta.data)
# DimPlot(sce.all.int, reduction = "umap", group.by = "orig.ident" )+
#   DimPlot(seuratObj, reduction = "umap" , group.by = "orig.ident")

###### step4:  看标记基因库 ######
# 原则上分辨率是需要自己肉眼判断，取决于个人经验
# 为了省力，我们直接看  0.1和0.8即可

table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.8) 

getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.5')
setwd('check-by-0.5')
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.8')
setwd('check-by-0.8')
sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
load(file = 'qc-_marker_cosg.Rdata')
top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10))) 

library(dplyr)  
sink(file = 'for_act.txt')
cat(paste0('cluster',colnames(marker_cosg$names),':',
           unlist(apply(marker_cosg$names,2,function(x){
             paste(head(x),collapse=',')
           })),'\n'))
sink()
#a=data.table::fread('ACT_Annotation results_top 1.txt')
write.csv(marker_cosg$names,file = 'marker_cosg.csv') 
setwd('../') 
getwd()

last_markers_to_check


###### step5: 确定单细胞亚群生物学名字 ######
# 一般来说，为了节省工作量，我们选择0.1的分辨率进行命名
# 因为命名这个步骤是纯人工 操作
# 除非0.1确实分群太粗狂了，我们就选择0.8 


## 如果是视网膜
if(F){
  
  data=read.table("./paper.retina-cell_marker.txt",sep="\t",header=T)
  head(data)
  data$gene <- str_to_title(data$gene)
  selected_genes = split(data$gene,data$cell)
  selected_genes
  selected_genes = lapply(selected_genes, str_to_title)
  selected_genes
  p3 <- DotPlot(object = sce.all.int, features = selected_genes, 
                group.by = 'RNA_snn_res.0.1',
                assay = "RNA") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
    )+theme(
      strip.text= element_text(angle = 90, vjust = 0.5, hjust=0.5)
    )
  p3
  ggsave(p3,filename = "check_retina_markers_by_RNA_snn_res.0.1.pdf",
         units = "cm",width = 36,height = 20) 
  p3 <- DotPlot(object = sce.all.int, features = selected_genes, 
                group.by = 'RNA_snn_res.0.8',
                assay = "RNA") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
    )+theme(
      strip.text= element_text(angle = 90, vjust = 0.5, hjust=0.5)
    )
  p3
  ggsave(p3,filename = "check_retina_markers_by_RNA_snn_res.0.8.pdf",
         units = "cm",width = 36,height = 20) 
} 

## 如果是皮肤
if(F){
  
  # keratinocyte：SFN
  #melanocyte1 "MLANA" 
  selected_genes=c('MUC4', 'PI3', 'SIX3', # nose
                   'SCGB1A1', 'TFF3',   'KRT1',   'KRT10',    
                   'KRT5', 'TP63',   'DLK2',
                   'MKI67', 'TOP2A', 'CDC20' ,
                   'DMD','PUS7','PGAP1','IFI44L',"MLANA" ,'SFN',
                   'MUC5AC' , 'MUC5B',# secretory cells
                   'FOXJ1', 'TPPP3',  'SNTN',  # multiciliated cells 
                   'DEUP1', 'FOXN4',  'CDC20B' # deuterosomal cells
  )
  p1 <- DotPlot(sce.all.int, features = unique(str_to_upper(selected_genes)),
                assay='RNA' ,group.by = 'RNA_snn_res.0.1'  )  + coord_flip() # +th
  
  p1
  
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.8,sce.all.int$orig.ident))
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  
}

## 如果是大脑
if(F){
  
  astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3") 
  endothelial = c("CLDN5", "ABCB1", "EBF1") 
  excitatory = c("CAMK2A", "CBLN2", "LDB2") 
  inhibitory = c("GAD1", "LHFPL3", "PCDH15") 
  microglia = c("C3", "LRMDA", "DOCK8") 
  oligodendrocytes = c("MBP", "PLP1", "ST18")
  # 下面的 OPC是 上面的 oligodendrocytes 的前体细胞 
  OPC='Tnr,Igsf21,Neu4,Gpr17'
  Ependymal='Cfap126,Fam183b,Tmem212,pifo,Tekt1,Dnah12'
  pericyte=c(  'DCN', 'LUM',  'GSN' ,'FGF7','MME', 'ACTA2','RGS5')
  
  gene_list = list(
    Astro = astrocytes,
    Endo = endothelial,
    Excit = excitatory,
    Inhib = inhibitory,
    Mic = microglia,
    Oligo = oligodendrocytes,
    OPC= str_to_upper(trimws(strsplit(OPC,',')[[1]])),
    Ependymal= str_to_upper(trimws(strsplit(Ependymal,',')[[1]])) ,
    peri = pericyte
  )
  gene_list = lapply(gene_list , str_to_title)
  
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.1') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('brain-0.1.pdf',width=12)
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('brain-0.8.pdf',width=12)
  
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.8,sce.all.int$orig.ident))
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  
  Astrocytes = c("Slc1a2", "Slc1a3","Gpc5", "Prex2","Atp1a2") 
  Epithelial = c("Kitl","Thsd4","Gbp3","Rps8","Tshz2") 
  Oligodendrosytes = c("Mbp","St18","Rnf220","Plp1","Slc24a2") 
  Ependymal = c( "Dnah6","Dnah12","Cfap44","Spag17","Spag16") 
  Oligodendrosyte_progenitor = c("Lhfpl3","Tnr","Dscam","Xylt1","Pcdh15") 
  Contaminating_neurons = c("Meg3","Snhg11","Syt1","Cntnap2","Nrg3") 
  Immune_cells = c( "Ctsd","Hexb","Inpp5d","Ccr5","Cx3cr1")
  Choroid_plexus_epithelial_cells = c("Htr2c","Ttr","Otx2os1","Enpp2","Trpm3")
  Meningeal_cells=c("Fbxl7","Eya2","Bnc2","Cped1","Tmtc1")
  Endothelial =c("Flt1","Slco1a4","Mecom","Adgrl4","Cldn5")
  
  gene_list = list(Astrocytes,Epithelial,Oligodendrosytes,Ependymal,Oligodendrosyte_progenitor,
                   Contaminating_neurons,Immune_cells,Choroid_plexus_epithelial_cells,Meningeal_cells,Endothelial)
  
  names(gene_list)=trimws( strsplit('Astrocytes,Epithelial,Oligodendrosytes,Ependymal,Oligodendrosyte_progenitor,
                   Contaminating_neurons,Immune_cells,Choroid_plexus_epithelial_cells,Meningeal_cells,Endothelial',',')[[1]])
  p_all_markers=DotPlot(sce.all.int, 
                        group.by = 'RNA_snn_res.0.1',
                        features = gene_list,
                        scale = T,assay='RNA' )+  
    theme(axis.text.x=element_text(angle=45,hjust = 1))
  p_all_markers
  ggsave('check_paper_markers_RNA_snn_res.0.1.pdf',
         height = 8,width = 12)
  
}

## 如果是肾脏 
if(F){
  # proximal tubule (PT)
  # parietal epithelial cells (PEC)
  # thick ascending limb (TAL)
  # distal tubule (DCT1, DCT2)
  # connecting tubule (CNT)
  # collecting duct (PC, ICA, ICB)
  # endothelial cells (ENDO)
  # glomerular cell types (MES, PODO)
  # fibroblasts (FIB)
  # and a small population of leukocytes (LEUK) 
  
  celltype.markers <- c("CUBN","HAVCR1","SLC5A1","SLC5A2", # PT and PT-VCAM1+ markers
                        "CFH", # PEC
                        "SLC12A1", # TAL NKCC2
                        "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                        "SCNN1G","TRPV5", # DCT2/CNT ENaC
                        "CALB1", # CNT
                        "AQP2", # PC
                        "ATP6V0D2", # ICA and ICB
                        "SLC4A1","SLC26A7", # ICA
                        "SLC26A4", # ICB
                        "NPHS1","NPHS2", # PODO
                        "PECAM1","FLT1","EMCN", # ENDO
                        "CLDN5", # GEC
                        "ITGA8","PDGFRB", # MES
                        "ACTA2","CALD1", # FIB
                        "PTPRC") # WBC
  p2 <- DotPlot(sce.all.int, features = celltype.markers) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p1 <- DotPlot(sce.all.int, features = unique(str_to_upper(celltype.markers)),
                assay='RNA' ,group.by = 'RNA_snn_res.0.1'  )  + coord_flip() # +th
  
  p1
  
  p2 = DotPlot( sce.all.int,  features = celltype.markers, 
                group.by = 'RNA_snn_res.0.1') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('kidney-0.1.pdf',width=12)
  p2 = DotPlot( sce.all.int,  features = celltype.markers, 
                group.by = 'RNA_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('kidney-0.8.pdf',width=12)
  
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.8,sce.all.int$orig.ident))
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  
}

# 如果是胰腺 
if(F){
  # acinar ROIs (CEL, CPA1, CELA2B, and PRSS3)
  # ductal ROIs (MUC5B and MUC20).
  selected_genes=  c("MUC5AC", "TFF1",
                     'CEL', 'CPA1', 'CELA2B',   'PRSS3',
                     'MUC5B' , 'MUC20',
                     "MUC6", "TFF2",
                     "CHGB", "CHGA", # enteroendocrine cells (marked with CHGA;
                     "PGA3","PGA4",
                     "MKI67","TOP2A","BIRC5",
                     "MUC2","ITLN1","SPINK4",
                     "APOA1","ALPI","FABP1","APOA4",
                     "CEACAM5","CEACAM6",
                     "CD79A","CD19",
                     "CD2", "CD3D","CD3E","CD3G",
                     "DCN","PDPN",
                     "ACTA2","ACTN2","MYL2","MYH2",
                     "VWF","ENG",
                     "TPSAB1",
                     "CXCL3","IL8",
                     "OLFM4","PIGR","EPHB2","SOX9",
                     "CD14","CD163","CD68","CSF1R")
  p1 <- DotPlot(sce.all.int, features = unique(str_to_upper(selected_genes)),
                assay='RNA' ,group.by = 'RNA_snn_res.0.1'  )  + coord_flip() # +th
  
  p1
  ggsave('paper-markers-pancr.pdf',height =  8)
  
  gene_list <- list(
    "β cells" = c("INS", "PCSK1", "G6PC2"),
    "α cells" = c("GCG"),
    "δ cells" = c("SST"),
    "ε cells" = c("GHRL"),
    "pancreatic progenitor cells" = c("SOX9"),
    "proliferation cells" = c("MKI67"),
    "EC cells" = c("FEV")
    #, "polyhormonal endocrine cells" = c("GCG", "INS")
  )
  
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.1') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('pancreatic-0.1.pdf',width=12)
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.8') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('pancreatic-0.8.pdf',width=12)
  
  
  pdf('orig.ident-vs-RNA_snn_res.0.1.pdf',height =  12 )
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  dev.off()
}


source('scRNA_scripts/lib.R')
sce.all.int = readRDS('2-harmony/sce.all_int.rds')
sp='human'
tmp=sce.all.int@meta.data
colnames(sce.all.int@meta.data) 

###查看默认Idents
Idents(sce.all.int)
table(sce.all.int$RNA_snn_res.0.1) 
###设置为RNA_snn_res.0.8继续后续的单细胞亚群注释
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
table(sce.all.int$RNA_snn_res.0.1)
# pbmc_small <- BuildClusterTree(object = sce.all.int)
# plot(Tool(object = pbmc_small, slot = 'BuildClusterTree'))
# plot(pbmc_small@tools$BuildClusterTree)


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
                    "HBB", "HBA1", "HBA2", #Erythrocytes
                    "CXCL8", "DCN", "SLPI", "CXADR",'FAP','CFD','KRT8', #Fibroblasts
                    "PECAM1", "VWF", "ENG", "MCAM", #Endothelial
                    "MUC5AC", "TFF2", "TFF1", "GKN1", "GKN2", #Gastric
                    "ALB", "APOC1", "APOA2", "ADH1B", "FGB", "FABP1", "ALDOB", #Hepatocytes
                    "APOA1", "APOA4", "ALPI", "RBP2", "FABP2", "CDH17", "REG4", #Intestinal-like cells
                    "EPCAM", "SOX9", "MET", "EGFR", "MSLN", "YAP1", #Malignant
                    "MKI67", "BIRC5", "CDK1",'TOP2A') #Proliferative

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
rownames(sce.all.int)
length(rownames(sce.all.int))
'KRT8' %in% rownames(sce.all.int)
p = DotPlot(sce.all.int, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p

sce = sce.all.int
#文中只取出了免疫细胞
sce.all.int = sce

###重要：注意人工肉眼对应，需要根据不同数据集修改分群所对应的细胞亚群
if(T){
  sce.all.int
  celltype=data.frame(ClusterID=0:27 ,
                      celltype=0:27 ) 
  #定义细胞亚群        
  celltype[celltype$ClusterID %in% c( 0,1,2,3,4,9,16  ),2]='T/NK'
  celltype[celltype$ClusterID %in% c(   ),2]='NK'
  celltype[celltype$ClusterID %in% c(  6 ),2]='B'
  celltype[celltype$ClusterID %in% c( 11,15,26 ),2]='Plasma'
  celltype[celltype$ClusterID %in% c( 22 ),2]='mast'
  celltype[celltype$ClusterID %in% c( 5,8,17,19,20,21 ),2]='Myeloid'
  # 
  celltype[celltype$ClusterID %in% c( 18 ),2]='Platelets'
  celltype[celltype$ClusterID %in% c(  ),2]='Erythrocytes'
  ###其实27亚群根据目前的标志物无法看出具体是什么亚群，回到上一步看一下每个亚群的特异性基因表达
  #如果还不确定，可以走一下自动注释大致对应一下
  celltype[celltype$ClusterID %in% c( 27 ),2]='Fibroblasts'     
  celltype[celltype$ClusterID %in% c( 12 ),2]='Endothelial'  
  celltype[celltype$ClusterID %in% c( 10,25 ),2]='Gastric'
  celltype[celltype$ClusterID %in% c( 23 ),2]='Hepatocytes'
  celltype[celltype$ClusterID %in% c(  ),2]='Intestinal-like cells'     
  celltype[celltype$ClusterID %in% c( 7,13,14,24 ),2]='Malignant'  
  celltype[celltype$ClusterID %in% c(  ),2]='Proliferative' 
  
  head(celltype)
  celltype
  table(celltype$celltype)
  sce.all.int@meta.data$celltype = "NA"
  
  ###上面选择的是RNA_snn_res.0.1，如果是别的resolution，是需要修改的
  for(i in 1:nrow(celltype)){
    sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
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
              group.by = "RNA_snn_res.0.8",label = T,cols = color)
umap1


umap2=DimPlot(sce.all.int, reduction = "umap",
              # split.by = "orig.ident",
              group.by = "celltype",label = T,repel = T,cols = color)

umap2
umap1+umap2




