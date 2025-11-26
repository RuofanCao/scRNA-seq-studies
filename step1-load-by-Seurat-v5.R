rm(list=ls())
source('scRNA_scripts/lib.R')
getwd()

###### step1: 导入数据 ######   
dir='GSE184242_RAW/scRNA1/' 
samples=list.files('GSE184242_RAW/scRNA1/','^GSM')
samples
library(tidyverse)
library(readr)

sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  ct=fread( file.path(dir,pro),data.table = F)
  ct[1:4,1:4]
  ct[nrow(ct),1:4]
  length(unique(ct$Symbol))
  ct=ct[!duplicated(ct$GENE),]
  rownames(ct)=ct$GENE
  ct=ct[,-1]
  ct[1:4,1:4]
  
  sce =CreateSeuratObject(counts =  ct ,
                          project =   gsub('_tagged.dge.txt.gz','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  print(sce)
  return(sce)
})


sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  gsub('_tagged.dge.txt.gz','',gsub('^GSM[0-9]*_','',samples) )     )
table(sce.all$orig.ident)

names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 

# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )
 
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

# # 分组信息,通常是隐含在文件名,样品名字里面
# phe=str_split( colnames(sce.all),'[-_]',simplify = T)
# head(phe)
# table(phe[,2]) 
# sce.all$orig.ident=phe[,2]
table(sce.all$orig.ident)
sce.all
# sce.all = subset(sce.all,downsample=10000)
# sce.all
head(sce.all@meta.data, 10)
table(sce.all$orig.ident)

sp='human'
# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
getwd()

sce.all.filt
###根据文章中信息进行了严格的过滤，后续处理时上述代码就行，这里设置为F
if(F){
tmp=subset(sce.all.filt,subset = percent_mito < 20 &
              nCount_RNA <10000 &
             nFeature_RNA <7500&
             nFeature_RNA>500)
tmp
sce.all=tmp
sort(table(tmp$orig.ident))
sce.all.filt =tmp
}


###### step3: harmony整合多个单细胞样品 ######
 
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../')
  
}



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


# pbmc_small <- BuildClusterTree(object = sce.all.int)
# plot(Tool(object = pbmc_small, slot = 'BuildClusterTree'))
# plot(pbmc_small@tools$BuildClusterTree)

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

########文献临床信息
#500，532,XX 是lung；506是ovary；516是kidney；
#D是droped 类器官MOS，T是tissue

table(sce.all.int$RNA_snn_res.0.1)
table(sce.all.int@active.ident) 
sce.all.int <- SetIdent(sce.all.int, value = 'RNA_snn_res.0.1')
table(sce.all.int@active.ident) 
phe = sce.all.int@meta.data
table(sce.all.int$orig.ident)

##增加样本来源
phe$sample = phe$orig.ident
table(phe$sample)
phe$sample <- ifelse(str_detect(phe$sample,"500|532|XX"),"lung",
                ifelse(str_detect(phe$sample,"506"),"ovary","kidney"))
table(phe$sample)


##增加样本类型
phe$type = phe$orig.ident
table(phe$type)
phe$type <- ifelse(str_detect(phe$type,"D"),"MOS","Tissue")
table(phe$type)

sce.all.int@meta.data = phe


#####注释亚群
#通常我们第一层次降维聚类分群：immune (CD45+,PTPRC),epithelial/cancer (EpCAM+,EPCAM),stromal (CD10+,MME,fibro or CD31+,PECAM1,endo)
#文中提供的都是常见的细胞群：上皮细胞（EPCAM、KRT19、CLDN4）、基质（PECAM1、CLO1A2、VWF）、增殖性（MKI67、STMN1、PCNA）、T（CD3D、CD3E、CD2）、B（CD79A，IGHG1，MS4A1），NK（KLRD1、GNLY、KLRF1）和髓系（CSF1R、CSF3R、CD68）细胞。
dir.create("./3-Celltype")
setwd("./3-Celltype")
scRNA=sce.all.int
source('../scRNA_scripts/mycolors.R')

genes_to_check = c('CD45','CD10','CD31',
                   'EPCAM','KRT19','CLDN4','SCGB1A1',  #上皮
                   'KRT19', 'MKI67', #proliferating epithelial
                   'PECAM1' , 'CLO1A2', 'VWF',  #基质
                   'DCN', #mesenchyme 
                   'CDH5', 'PECAM1', 'VWF','CLDN5',  #内皮
                   'LUM' , 'FGF7', 'MME',  #成纤维
                   'CD3D', 'CD3E', 'CD8A', 'CD4','CD2', #T,proliferating lymphoid (CD3E, MKI67)
                   'AIF1', 'C1QC','C1QB','LYZ',  #巨噬
                   'MKI67', 'STMN1', 'PCNA',  #增殖
                   'CPA3' ,'CST3', 'KIT', 'TPSAB1','TPSB2',#肥大
                   'GOS2', 'S100A9','S100A8','CXCL8', #中性粒细胞
                   'KLRD1', 'GNLY', 'KLRF1','AREG', 'XCL2','HSPA6', #NK
                   'MS4A1','CD19', 'CD79A','IGHG1','MZB1', 'SDC1',  #B
                   'IGHD',  #MALT B
                   'CSF1R', 'CSF3R', 'CD68') #髓系

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(scRNA, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('check_markers.pdf',height = 11,width = 11)

####构建左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(scRNA)
axes <- result$axes
label <- result$label

umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
umap
ggsave('RNA_snn_res.0.1_umap.pdf',width = 9,height = 7)

##图中的标签和方框都是可以自定义的，例如下面这副删掉label和label.box
patient_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                      group.by = "sample") +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
patient_umap
ggsave('sample_umap.pdf',width = 9,height = 7)

umap+patient_umap


#####细胞生物学命名
celltype=data.frame(ClusterID=0:10,
                    celltype='Epithelial') 

# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c( 6),2]='Myeloid'
celltype[celltype$ClusterID %in% c( 5 ),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(0),2]='T'
celltype[celltype$ClusterID %in% c(9  ),2]='Endothelial'

table(scRNA@meta.data$RNA_snn_res.0.1)
table(celltype$celltype)

scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)
Idents(scRNA)=scRNA$celltype

table(scRNA$celltype)
# 过滤无法命名的 
scRNA = scRNA[,!scRNA$celltype %in% 0:100]
table(scRNA$celltype)

# 如果前面成功的给各个细胞亚群命名了
# 就可以运行下面的代码
scRNA
scRNA$celltype[scRNA$celltype=='']='unknown'
table(scRNA$celltype)
scRNA=scRNA[,scRNA$celltype !='unknown']
sce.all.int=scRNA
sp='human'
if("celltype" %in% colnames(scRNA@meta.data ) ){
  
  sel.clust = "celltype"
  sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
  table(sce.all.int@active.ident) 
  
  dir.create('check-by-celltype')
  setwd('check-by-celltype')
  source('../../scRNA_scripts/check-all-markers.R')
  setwd('../') 
  getwd()
  phe=sce.all.int@meta.data
  save(phe,file = 'phe.Rdata')
  pdf('celltype-vs-orig.ident.pdf',width = 10)
  gplots::balloonplot(table(sce.all.int$celltype,sce.all.int$orig.ident))
  dev.off()
}


# sce.all.int
# table(Idents(sce.all.int))
# sce.main.1000 = subset(sce.all.int,downsample=1000)
# sce.main.1000
# save(sce.main.1000,file = 'sub.sce.Rdata')

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)
celltype_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                       group.by = "celltype",label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
celltype_umap
ggsave('umap_by_celltype.pdf',width = 9,height = 7)


type_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                      group.by = "type") +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
type_umap
ggsave('type_umap.pdf',width = 9,height = 7)


type_umap + patient_umap+celltype_umap 
ggsave('combine_umap.pdf',width = 15,height = 7)


saveRDS(scRNA, "sce_celltype.rds")

table(scRNA$type)
scRNA$type = factor(scRNA$type,levels = c('Tissue','MOS'))
type_splitumap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                  split.by ="type",group.by = "celltype") 

type_splitumap


setwd('../')
