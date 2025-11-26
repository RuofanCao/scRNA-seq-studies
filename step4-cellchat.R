rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(clustree)
dir.create('4-cellchat/')
setwd('4-cellchat/')
Allsce=readRDS( "../2-harmony/sce.all_int.rds")
load('../phe.Rdata')
Allphe = phe
identical(colnames(Allsce),rownames(Allphe))
../2-harmony/sce.all_int.rds
source('../scRNA_scripts/lib.R')
source('../scRNA_scripts/mycolors.R')

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(CellChat)
library(igraph)
#devtools::install_github("jinworks/CellChat")

TNKsce=readRDS( "../step3-sub-cluster/sub-Tcells/delete_unknown_harmony/sce.all_int.rds")
load('../step3-sub-cluster/sub-Tcells/phe.Rdata')
TNKphe = phe
identical(colnames(TNKsce),rownames(TNKphe))
TNKsce@meta.data = TNKphe
table(TNKsce$celltype)

NKsce=readRDS( "../step3-sub-cluster/sub-Tcells/gene_score/NKsce.rds")
NKphe = NKsce@meta.data
table(NKsce$celltype)

DCsce = readRDS( "../step3-sub-cluster/sub-myeloids/DC/delete_unknown_harmony/sce.all_int.rds")
load('../step3-sub-cluster/sub-myeloids/DC/phe.Rdata')
DCphe = phe
identical(colnames(DCsce),rownames(DCphe))
DCsce@meta.data = DCphe
table(DCsce$celltype)


Myeloidsce = readRDS( "../step3-sub-cluster/sub-myeloids/delete_unknown_harmony/sce.all_int.rds")
load('../step3-sub-cluster/sub-myeloids/phe.Rdata')
Myeloidphe = phe
identical(colnames(Myeloidsce),rownames(Myeloidphe))
Myeloidsce@meta.data = Myeloidphe
table(Myeloidsce$celltype)

Treg_phe = sce@meta.data
TNK_phe = TNKsce@meta.data
All_phe = Allsce@meta.data
colnames(All_phe)
All_phe$celltype3 = All_phe$celltype
table(All_phe$celltype3)

####把Treg中的注释转换到TNK中
table(NKphe$celltype)
table(TNKphe$celltype)
table(DCphe$celltype)
table(Myeloidphe$celltype)
table(Allphe$celltype)
for (i in 1:nrow(DCphe)) {
  Myeloidphe$celltype[rownames(Myeloidphe) == rownames(DCphe)[i]] <- DCphe$celltype[i]
}
##检查一下是否对应
a = Myeloidphe[ Myeloidphe$celltype == 'cDC2',]
identical(rownames(a),rownames(DCphe[DCphe$celltype == 'cDC2',]))

for (i in 1:nrow(TNKphe)) {
  Allphe$celltype[rownames(Allphe) == rownames(TNKphe)[i]] <- TNKphe$celltype[i]
}

for (i in 1:nrow(Myeloidphe)) {
  Allphe$celltype[rownames(Allphe) == rownames(Myeloidphe)[i]] <- Myeloidphe$celltype[i]
}

table(Allphe$celltype)

##这样合并之后有一些在后续分群中被delete的亚群无法合适的再分群


Allsce@meta.data = Allphe

sce = Allsce
colnames(Allphe)
table(sce$celltype)
Idents(sce) = sce$celltype
colnames(sce@meta.data)

names(sce@assays$RNA@layers)
#sce[["RNA"]]$counts 
# Alternate accessor function with the same result
#LayerData(sce, assay = "RNA", layer = "counts")
sce
gene = c("XCR1",  "XCL1", "XCL2")
##vlnplot
p1 <- VlnPlot(sce,features =  gene,
              group.by  = "celltype",
              flip = T,stack = T )
p1 + NoLegend()

Allphe = sce@meta.data

save(Allphe,file = 'celltype_input_phe.Rdata')


# CellChat需要两个用户输入:一个是细胞的基因表达数据，另一个是细胞标签。
data.input <- GetAssayData(sce, layer = 'data')
meta <- sce@meta.data[,c("orig.ident","celltype")]
colnames(meta) <-  c("group","celltypes")
table(meta$celltypes)

identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "celltypes")
levels(cellchat@idents)

####设置配-受体相互作用数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# 使用所有CellChatDB数据进行细胞-细胞通信分析。
CellChatDB.use <- CellChatDB 

# 在构建的cellchat中设定需要使用的数据库
cellchat@DB <- CellChatDB.use


####预处理细胞-细胞通讯分析的表达数据
cellchat <- subsetData(cellchat) 
# future::plan("multisession", workers = 1) # do parallel
#devtools::install_github('immunogenomics/presto')
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# 默认情况下,cellchat使用object@data.signaling进行网络推断
# 此外提供了projectData函数,将基因投射到PPI，开发者说PPI并不会或导致很少的伪通讯
cellchat <- smoothData(cellchat, adj = PPI.human)


####计算通信概率并推断通信网络
#Compute the communication probability and iC3D1er cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#过滤，有些细胞亚群细胞数量过少被过滤掉之后会造成数据不一致而无法进行后续分析
cellchat <- filterCommunication(cellchat, min.cells = 10)

####提取信号通路水平的细胞通讯表：
cellchat <- computeCommunProbPathway(cellchat) #计算信号通路水平上的通讯概率
df.netp <- subsetCommunication(cellchat, slot.name = 'netP') #得到信号通路水平细胞通讯表
head(df.netp)

####细胞互作关系展示：
#计算细胞对间通讯的数量和概率强度
cellchat <- aggregateNet(cellchat)
#数据提取，subsetCommunication函数
df.net <- subsetCommunication(cellchat)

save(df.net,file = "df.net.Rdata")

head(df.net)


####可视化细胞-细胞间通信网络
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("Number_weights.pdf", width = 15, height = 9)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#互作的数量和互作的权重。其中每个颜色代表了不同的细胞，箭头代表了顺序，线的粗细代表了数量/权重。


####细分亚组的circle图
mat <- cellchat@net$weight
par(mfrow = c(2,5), xpd=TRUE,mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i]
  )
}

##单独查看和文章中所挑选的trNK对应的NK2
par(mfrow = c(1,1), xpd=TRUE,mar = c(1, 1, 1, 1))
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[26, ] <- mat[26, ]

pdf("NK_circle_plot.pdf", width = 13, height = 9)
netVisual_circle(mat2, vertex.weight = groupSize, 
                 weight.scale = T, arrow.size=0.05,
                 arrow.width=1, edge.weight.max = max(mat), 
                 title.name = rownames(mat)[26])
dev.off()


#### 计算配-受体对信号通路的贡献并可视化
cellchat@netP$pathways
levels(cellchat@idents) 
pathways.show <- "XCR"
netAnalysis_contribution(cellchat, signaling = pathways.show)

# 可视化由单个配体-受体对介导的细胞间通讯
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show,
                            geneLR.return = FALSE)
LR.show <- pairLR[1,] # show one ligand-receptor pair

# plot
#vertex.receiver = seq(1,3) # a numeric vector
pdf("XCR_XCL1_plot.pdf", width = 9, height = 7)
netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = pairLR[1,], 
                     layout = "circle")
dev.off()

pdf("XCR_XCL2_plot.pdf", width = 9, height = 7)
netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = pairLR[2,], 
                     layout = "circle")
dev.off()


####基于配-受体结果进一步可视化

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:3), 
                 targets.use = c(4:10), remove.isolate = FALSE)
ggsave("bubbleplot.pdf",width = 8,height = 16)

# 还可以增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = seq(23:29), 
                 targets.use = c(10:12), 
                 signaling = c("XCR"),
                 remove.isolate = FALSE)
ggsave("XCR_dotplot.pdf",width = 13,height = 8)


####识别对某些细胞群的输出或输入信号贡献最大的信号
# 特定的signaling
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
pathway = cellchat@netP$pathways
###选取了前30条，发现我们关注的XCR并不在此列，所以特地加上了一条“XCR”
pathways.show =  c(pathway[1:30],'XCR')

htout <- netAnalysis_signalingRole_heatmap(cellchat, 
                                           pattern = "outgoing",
                                           signaling = pathways.show)
htout

htcome <- netAnalysis_signalingRole_heatmap(cellchat, 
                                            pattern = "incoming",
                                            signaling = pathways.show)
htcome
pdf("outgoing_incoming_heatmap.pdf", width = 15, height = 8)
htout + htcome
dev.off()


setwd('../')
