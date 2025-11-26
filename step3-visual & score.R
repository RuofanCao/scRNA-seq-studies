rm(list=ls())
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(ggsci) 
library(patchwork) 
library(ggsci)
library(ggpubr)
library(RColorBrewer) 

getwd()
dir.create('step3-visual & score/')
setwd('step3-visual & score/') 

## 1. 载入 ---- 
sce.all = readRDS('../3-Celltype/sce_celltype.rds')
sp='human' 
load('../step2-paper-figures/phe.Rdata') 
identical(rownames(phe) , colnames(sce.all))
sce.all@meta.data = phe
# sel.clust = "celltype"
# sce.all <- SetIdent(sce.all, value = sel.clust)
Idents(sce.all) = sce.all$celltype
table(sce.all@active.ident) 
DimPlot(sce.all) 

table(Idents(sce.all))  
unique(sce.all$celltype)
length(unique(sce.all$celltype))
# sce.all=sce.all[,sce.all$celltype != '13']

###可以设置因子使每个图形的细胞亚群顺序一致
ord = c("Epithelial", "Fibroblast" , "T","Myeloid","Endothelial")
all(ord %in% unique(sce.all$celltype) )
sce.all$celltype = factor(sce.all$celltype ,levels = ord)
table(sce.all$celltype)
DimPlot(sce.all)

marker <- c('CD3E','KRT19','PECAM1','CD68','MKI67','CD79A','DCN','KIT')
FeaturePlot(sce.all,features = marker,
            cols = c("lightgrey" ,"#DE1F1F"),
            ncol=4,order = T,raster = T )
# 注意 raster参数 ,超过100,000个细胞时，设置为TRUE
ggsave('FeaturePlot_marker.pdf',width = 12,height = 8)


###使每一个基因的颜色阈值范围调整一致[0-6]
p1 <- FeaturePlot(sce.all, features = marker, combine = FALSE,
                  ncol=4,order = T,raster = T  )
#colours = c('lightgrey', "#DE1F1F")
fix.sc <- scale_color_gradientn( colours = c('lightgrey', "#DE1F1F"),  limits = c(0, 6))
#+NoLegend()+NoAxes()
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
ggsave('FeaturePlot_marker1.pdf',width = 12,height = 8)

##批量画基因,注意图例的范围不同(这里遇到过一个大bug -- 谨慎)
FeaturePlot(sce.all, features =marker, 
            cols = c("lightgrey", 'red'),
            order = T,raster = T ,
            ncol = 3 ) & NoLegend() & NoAxes() & theme(
              panel.border = element_rect(color = "black", size = 1)
            )
ggsave('FeaturePlot_marker2.pdf',width = 12,height = 8)

#Featureplot还可以把两个基因画在同一个图中,看右上角可以发现黄色越深的地方两个基因叠加越多
FeaturePlot(sce.all, features = c('S100A9','S100A8'),
            cols = c("lightgrey", "green", "orange"),
            order = T,raster = T ,
            blend=T,blend.threshold=0)
ggsave('FeaturePlot_marker3.pdf',width = 12,height = 7)


#Featureplot还可以把三个基因画在同一个图中
# 提取umap坐标
umap_df <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
umap_df$cluster <- as.factor(sce.all$celltype)
head(umap_df)
# 提取基因表达数据并与umap坐标合并
# gene_df <- as.data.frame(GetAssayData(object = sce.all, slot = "data")[c('S100A9','S100A8','CXCL8'), ])
# 下面的 sce.all@assays$RNA$counts 如果是在Seurat的v4里面需要修改为
# sce.all@assays$RNA@counts， 千万有留意版本问题
gene_df <-  sce.all@assays$RNA$counts[c('S100A9','S100A8','CXCL8'),]
gene_df <-  t(as.data.frame(gene_df))
library(ggnewscale)
##如果有运行错误，可能是电脑和服务器有什么差异（同样代码遇到过bug，暂未查原因
merged_df <- cbind(gene_df, umap_df)
head(merged_df)
colnames(merged_df)
source('../scRNA_scripts/mycolors.R')
plot = ggplot(merged_df, vars = c("umap_1", "umap_2", 
                           'S100A9','S100A8','CXCL8'), 
       aes(x = umap_1, y = umap_2, colour = S100A9)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = S100A8), size=0.3, alpha=0.7) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CXCL8), size=0.3, alpha=0.1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish)+
  theme_classic()

plot2 <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T,label.box = T, cols = color,repel = T) +
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))

plot + plot2
ggsave('FeaturePlot_marker4.pdf',width = 14,height = 7)



####密度图可视化umap####
# BiocManager::install("Nebulosa")
# install.packages("ggrastr")
library(Nebulosa)
library(ggrastr)
##免疫检查点
checkpoint <- c("CD274", "PDCD1", "TGFB1")
##常见marker
marker <- c('EPCAM','CDH1','FCER1A','LYZ','CD3E','IL7R','FAP','PDGFRA')

MOSsce = sce.all[, sce.all$type %in% c( 'MOS' )]
Tissuesce = sce.all[, sce.all$type %in% c( 'Tissue' )]

plots_MOSsce <- list()
plots_Tissuesce <- list()

for (target in checkpoint) {
  plots_MOSsce[[target]] <- plot_density(MOSsce, target, pal = "magma", reduction = 'umap')
  
  plots_Tissuesce[[target]] <- plot_density(Tissuesce, target, reduction = 'umap')
}

combined_plots <- c(plots_MOSsce, plots_Tissuesce)
names(combined_plots) <- targets

pdf('plot1.pdf',height =  10,,width=14 )
combined_all <- ggarrange(plotlist = combined_plots, ncol = 3, nrow = 2)
combined_all
dev.off()



for (target in marker) {
  plots_MOSsce[[target]] <- plot_density(MOSsce, target, pal = "magma", reduction = 'umap')
  
  plots_Tissuesce[[target]] <- plot_density(Tissuesce, target, reduction = 'umap')
}

pdf('plot2.pdf',height =  13,,width=10 )
combined_MOSsce <- ggarrange(plotlist = plots_MOSsce, ncol = 2, nrow = 4)
combined_MOSsce
dev.off()

pdf('plot3.pdf',height =  13,,width=10 )
combined_Tissuesce <- ggarrange(plotlist = plots_Tissuesce, ncol = 2, nrow = 4)
combined_Tissuesce
dev.off()





#####基因集评分 score ----
rm(list=ls())
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(clusterProfiler)
source('../scRNA_scripts/lib.R')
source('../scRNA_scripts/mycolors.R')
dir.create("score")
setwd('score/')
sce = readRDS('../../3-Celltype/sce_celltype.rds')
sce
Idents(sce) = sce$celltype
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

##这里可以选择需要哪一类型
# H：包含了由多个已知的基因集构成的超基因集，每个H 类别的基因集都对应多个基础的其他类别的基因集
# C1：包含人类每条染色体上的不同cytoband区域对应的基因集合。根据不染色体编码进行二级分类
# C2：包含了在线通路数据库、PubMed 出版物和领域专家知识的精选基因集
# C3：基于对 miRNA 种子序列和预测的转录因子结合位点的基因靶点预测的调控靶基因组
# C4：包含了计算机软件预测出来的基因集合，主要是和癌症相关的基因
# C5：包含了由相同GO术语注释的基因集
# C6：包含了致癌特征基因集：直接从来自癌症基因扰动的微阵列基因表达数据中定义
# C7：包含了免疫特征基因集：代表免疫系统内的细胞状态和扰动
# C8：包含了细胞类型特征基因组：从人体组织单细胞测序研究中确定的簇标记中收集
#genesets <- msigdbr(species = "Homo sapiens", category = "C2") 
genesets <- msigdbr(species = "Homo sapiens") 

a = genesets[str_detect(genesets$gs_name,'HYPOXIA'),]

unique(a$gs_name)

m_hypoxia = a[a$gs_name == 'HALLMARK_HYPOXIA',]

#评分计算
#BiocManager::install("UCell")
library(UCell)
m_hypoxia <- list(m_hypoxia$gene_symbol)#将基因整成list
names(m_hypoxia)[1] <- 'hypoxia'
?AddModuleScore_UCell
hypoxia_score <- AddModuleScore_UCell(sce,
                                      features=m_hypoxia,
                                      name="hypoxia")
colnames(hypoxia_score@meta.data)
###改不改名字都没关系，和后面对应就行
colnames(hypoxia_score@meta.data)[19] <- 'hypoxia' 

library(ggpubr)
###最好最开始就存好临床信息
phe = hypoxia_score@meta.data
df<- FetchData(hypoxia_score,vars = c("type","hypoxia"))
table(df$type)

colnames(df)
colnames(df)[2] = "Hallmark_hypoxia_score"

# kruskal_result <- kruskal.test(Hallmark_hypoxia_score ~ Procedure, data = df)
# p_value <- kruskal_result$p.value
#小提琴图
p1 = ggplot(df,aes(x=type,y=Hallmark_hypoxia_score,fill=type))+
  geom_violin(color='black',size=1)+#小提琴
  theme_classic() + 
  theme(text = element_text(size=10, colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+ 
  theme(legend.position="none") +  
  xlab(NULL)+
  scale_fill_manual(values = mycolors) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +#箱线图
  stat_compare_means(method = "kruskal.test", #统计方法
                     aes(label = "p.format"), #显示方式
                     label.x = 1, label.y = 0.4,#位置
                     size = 5) #大小
p1
ggsave('score_violinplot.pdf',width = 9,height = 7)

setwd('../')


