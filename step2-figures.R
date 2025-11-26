rm(list=ls())
options(stringsAsFactors = F)

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
dir.create('figures/')
setwd('figures/') 

## 1. 载入，并且理解数据 ---- 
sce.all = readRDS('../2-harmony/sce.all_int.rds')
sp='human' 
load('../phe.Rdata')
identical(rownames(phe) , colnames(sce.all))
sce.all@meta.data = phe
table(sce.all$celltype)
sel.clust = "celltype"
sce.all <- SetIdent(sce.all, value = sel.clust)
table(sce.all@active.ident) 
DimPlot(sce.all) 

table(Idents(sce.all))  
unique(sce.all$celltype)
length(unique(sce.all$celltype))
# sce.all=sce.all[,sce.all$celltype != '13']

###可以设置因子使每个图形的细胞亚群顺序一致
ord = c("Epithelial"  ,  "Macrophages"  , "CD4T"   ,      
        "Neutrophils"  , "Acinar"   ,     "Monocytes" ,   
         "Proliferative", "NK"     ,       "B"   ,         
         "DC"     ,       "Mast"  ,        "Plasma"  ,     
         "CD8T"   ,       "Fibroblasts" ,  "Ductal"  ,     
         "Endothelial" ,  "Endocrine")
#ord = unique(sce.all$celltype) 
all(ord %in% unique(sce.all$celltype) )
sce.all$celltype = factor(sce.all$celltype ,levels = ord)
table(sce.all$celltype)
DimPlot(sce.all) 

## 2. 检查内置分群 ----
sce = sce.all
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
##不加label
p1_umap=DimPlot(sce, reduction = "umap",
                # split.by = "orig.ident",
                group.by = "celltype",label = T,label.size = 0,cols = color)
p1_umap
ggsave(p1_umap,filename = "umap_by_celltype_without_label.pdf",
       width = 9,height = 7)

##加label，加label.box
p2_umap=DimPlot(sce, reduction = "umap",
                # split.by = "orig.ident",
                group.by = "celltype",label = T,label.box = T, repel = T,cols = color)

p2_umap
ggsave(p2_umap,filename = "umap_by_celltype_with_label.pdf",
       width = 9,height = 7)

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  

##设置左下角坐标轴，这里使用的是umap，如要使用tsne，要去Bottom_left_axis.R中修改umap为tsne
source('../scRNA_scripts/Bottom_left_axis.R')
axe = left_axes(sce)
axes = axe$axes
label = axe$label
p3_umap =DimPlot(sce, reduction = "umap",cols = color,pt.size = 1,
                     group.by = "celltype",label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')
p3_umap
ggsave(p3_umap,filename = "umap_by_celltype_with_left_axes.pdf",
       width = 9,height = 7)

# ord = c( "cycle" ,  'endo' 
#          , 'Microglial',"Mac" ,"Mono","NK" ,'Tcells','Bcells')
selected_genes <-  c(    'MKI67' , 'TOP2A',  
                         'EPCAM' , 'KRT19','KRT7', # epi 
                         
                         'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                      'PECAM1', 'VWF',  ## endo 
                     'CD68', 'CD163', 'CD14',  
                     'TPSAB1' , 'TPSB2',
                     'C1QA',  'C1QB', # mac  # myeloids  
                     'S100A9', 'S100A8', 'MMP19',# monocyte
                     "GNLY","GZMA","NKG7","GZMK", # NK
                     'CD3D', 'CD3E', 'CD4','CD8A',  # Tcells 
                     'CD19', 'CD79A', 'MS4A1' , # Bcells 
                     'IGHG1', 'MZB1', 'SDC1'
                   
) 
selected_genes <-  unique(str_to_upper(selected_genes))
p <- DotPlot(sce, features = selected_genes,
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="check_marker_by_celltype.pdf",width = 12 ,height = 10)



## 3. 各式各样的umap  ---- 
colnames(sce@meta.data)
table(sce$orig.ident)
sce$tissue = ifelse(grepl('AdjNorm',sce$orig.ident),'AdjNorm','PDAC')
table(sce$tissue)
####需要根据每篇文章中的自定义临床信息
###重要：这里的group需根据自己想要的临床信息列设置
# library(stringr)
# phe = sce@meta.data
# table(phe$orig.ident) 
#sce$group = ifelse(grepl('Control',phe$orig.ident),'control','treat')
sce$group = sce$tissue
table(sce$group) 
phe = sce.all@meta.data
#注意颜色的修改，看个人喜欢什么色系
source('../scRNA_scripts/mycolors.R')
plot1 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.3",label = T,cols = mycolors,repel = T)+  
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
plot2 <- DimPlot(sce, reduction = "umap", group.by = "celltype",label = F,cols = mycolors,repel = T) +
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
#cols = sample(color,10)
plot3 <- DimPlot(sce, reduction = "umap", group.by = "group",label = T,cols = mycolors,repel = T) +
  theme(plot.title = element_text(size=12),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10))

colnames(phe)
##如果多个临床信息一起，可以循环
if(F){
  
  p_list = lapply(c('tissue','patient','tissue_p'), function(x){
    DimPlot(sce, reduction = "umap", group.by = x ,
            label = T,cols = color,repel = T) +
      theme(plot.title = element_text(size=12),
            axis.title = element_text(size=10),
            axis.text = element_text(size=10),
            legend.text = element_text(size=10))
  })
  p_list[[1]]
  
  fig1 = (p_list[[1]]+p_list[[2]])/(p_list[[3]])
}
fig1
ggsave(fig1,filename = "umap_0.3_celltype_group.pdf",units = "cm",width = 22,height = 12)
fig2 <- plot1|plot2|plot3
fig2
ggsave(fig2,filename = "umap_cli_group.pdf",width = 15,height = 7)
 

sce$group = phe$tissue_NT
plot4 <- DimPlot(sce, reduction = "umap", group.by = "celltype",split.by = "group",
                 label = F,cols = color) + theme(plot.title = element_text(size=12),
                                                 axis.title = element_text(size=10),
                                                 axis.text = element_text(size=10))

plot4
ggsave(plot4,filename = "plot4.pdf",units = "cm",width = 22,height = 10)


#UMAP更容易将相似的亚群聚集在一块，而tSNE则是注重将不同亚群区分开。两种方法在选择上没有绝对的好坏，只要能体现目标亚群间的差异性即可。
plot5 <- DimPlot(sce, reduction = "tsne", group.by = "celltype",split.by = "group",
                 ncol=3,label = F,cols = color) + 
  facet_wrap(~group, nrow = 2)+
  theme(plot.title = element_text(size=12),
                                                 axis.title = element_text(size=10),
                                                 axis.text = element_text(size=10))

plot5


fig3 <- (plot1|plot2|plot3)/ plot4 + plot_layout(heights = c(1:2))
fig3
ggsave(fig2,filename = "fig1-2.pdf",units = "cm",width = 22,height = 15)


plot5 <- DimPlot(sce, reduction = "umap", group.by = "celltype",
                 split.by = "group",ncol = 5,label = F,cols = color) 
plot5
ggsave(plot5,filename = "fig2-2.pdf",
       width = 14,height = 7)


## 4. 已知的基因  ----
pro = 'fig3'
library(dplyr)  
library(paletteer)

cg = c( 'PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
        'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
        'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
        'IFNG', 'CCL4', 'CCL3' ,
        'PRF1' , 'NKG7',  'KLRB1','NCR1', # NK 
        'CD19', 'CD79A', 'MS4A1' , 
        'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
        'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
        'IGHG1', 'MZB1', 'SDC1',
        'CD68', 'CD163', 'CD14', 
        'MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7', 
        'MAF','APOE','FOLR2','RELB','BST2','BATF3',
        'TPSAB1' , 'TPSB2',  # mast cells,
        'RCVRN','FPR1' , 'ITGAM' ,
        'C1QA',  'C1QB',  # mac
        'S100A9', 'S100A8', 'MMP19',# monocyte
        'FCGR3A','JCHAIN',
        'LAMP3', 'IDO1','IDO2',## DC3 
        'CD1E','CD1C','Cd209a','Clec9a', # DC2
        'FGF7','MME', 'ACTA2', ## fibo 
        'DCN', 'LUM',  'GSN' , ##  fibo 
        'PDGFRB', 'PDGFRA',
        'CSPG4','GJB2', 'RGS5','ITGA7',
        'MKI67' , 'TOP2A', 
        'PECAM1', 'VWF',  ## endo 
        'EPCAM' , 'KRT19',  
        'AGER','SFTPA1','SCGB1A1','KRT17','TPPP3',
        'KRT4','KRT14','KRT8','KRT18',
        'PROM1', 'ALDH1A1' )
cg
library(stringr)
genes_to_check = unique(str_to_upper(cg))

# 基本上代替热图的小提琴图
table(Idents(sce.all))
p_all_markers=DotPlot(sce.all, 
                      group.by  = "celltype",
                      features = genes_to_check,
                      scale = T,assay='RNA' ) + coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggplot2::ggsave(paste0(pro,'_all_genes_to_check-DotPlot.pdf'),
                height = 15,width = 6)

source('../scRNA_scripts/mycolors.R')
p1 <- VlnPlot(sce.all,features =  selected_genes,
              group.by  = "celltype",
              flip = T,stack = T,cols = my36colors)
p1 + NoLegend()
ggplot2::ggsave(paste0(pro,'_genes_to_check-VlnPlot-heatmap.pdf'),height = 12,width = 6)



## 4. COSG的基因  ---- 
table(Idents(sce))
table(Idents(sce.all))

###COSG比findallmarker更加快速
if(F){
  Idents(sce)
  table(sce$celltype)
  Idents(sce) = sce$celltype
  
  if (!file.exists('sce.markers.csv')) {
    sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)
    write.csv(sce.markers,file='sce.markers.csv')
  } else {
    
    sce.markers = read.csv('sce.markers.csv',row.names = 1)
  }
}


if(T){
  sce= sce.all
  pro = 'cosg_celltype_'
  library(COSG)
  marker_cosg <- cosg(
    sce,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  # width <-0.006*dim(sce)[2];width
  # height <- 0.25*length(top_10)+4.5;height
  
  width <- 15+0.5*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_10);height
  
  # levels(  Idents(sce) )  = 0:length(levels(  Idents(sce) ))
  DoHeatmap( subset(sce,downsample=100), top_10 , 
             size=3)
  
  ggsave(filename=paste0(pro,'DoHeatmap_check_top10_markers_by_clusters.pdf') ,
         # limitsize = FALSE,
         units = "cm",width=width,height=height)
  width <- 8+0.6*length(unique(Idents(sce)));width
  height <- 8+0.2*length(top_10);height
  DotPlot(sce, features = top_10 ,
          assay='RNA'  )  + coord_flip() +FontSize(y.text = 4)
  ggsave(paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'),
         units = "cm",width=width,height=height)
  
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
  width <- 15+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  
  DoHeatmap( subset(sce,downsample=100), top_3 ,
             size=3)
  ggsave(filename=paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf') ,
         units = "cm",width=width,height=height)
  
  width <- 8+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  DotPlot(sce, features = top_3 ,
          assay='RNA'  )  + coord_flip()
  ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),width=width,height=height)
  
  
}
 

getwd()
pro = 'cosg_celltype_'
load(file = paste0(pro,'_marker_cosg.Rdata'))

top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3))) 
sce.Scale <- ScaleData(subset(sce.all,downsample=100),features =  top_3  )  
table(sce.Scale$celltype)
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
 
table(Idents(sce.Scale)) 

ll =  as.list(as.data.frame(apply(marker_cosg$names,2,head,3)))
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
ll
###如果raster=T热图就会有点模糊
DoHeatmap(sce.Scale,
          features = unlist(ll),
          group.by = "celltype",
          raster = F,
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(filename = "fig4-top5-marker_pheatmap.pdf",
       width = 15,height = 8)

head(marker_cosg)
## Top10 genes
library(dplyr)  
sink(file = 'for_act.txt')
cat(paste0('cluster',colnames(marker_cosg$names),':',
           unlist(apply(marker_cosg$names,2,function(x){
             paste(head(x),collapse=',')
           })),'\n'))
sink()
#a=data.table::fread('ACT_Annotation results_top 1.txt')
write.csv(marker_cosg$names,file = 'marker_cosg.csv') 

### 4.1 . 针对top基因绘制样本分类的小提琴图 ---- 
ll
# sce$group = sce$tissue_p
table(sce$group)
pdf('fig4-compare_sub_top10_by_VlnPlot.pdf')
lapply(names(ll), function(x){
  cg = ll[[x]]
  VlnPlot(sce,features = cg,group.by = "celltype",split.by = "group",flip = T,stack = T,
          split.plot =T )+scale_fill_d3()
  
})
dev.off()


### 4.2 . 针对top基因绘制FeaturePlot ---- 
pdf('fig4-FeaturePlot-top10.pdf',width = 8,height = 8)
lapply(names(ll), function(x){
  # x=1
  cg = ll[[x]]
  FeaturePlot(sce,features = cg ,ncol = 3,order = T,raster = T ) & 
    theme(plot.title = element_text(size=12),
          axis.title = element_text(size=10),
          axis.text = element_text(size=10)) & 
    scale_colour_continuous(low = 'lightgrey',high = 'red') 
  
})
dev.off()


### 4.3 .. 针对top基因绘制小提琴图 ----
#基因可以自定义例如gene = c('EPCAM','CD3D','CD79A','CD79B','MNDA','LYZ','PECAM1','DCN')
p1 <- VlnPlot(sce ,features =  unique(unlist(ll)),flip = T,
              stack = T,
              group.by = "celltype",
              split.by = "celltype")+
  scale_color_npg()

p1
ggsave(p1,filename = "fig4_Vlnplot_top5.pdf",width = 8,height = 25)




## 5 看比例  ---- 
#load('../phe.Rdata')
head(phe) 
phe = sce@meta.data
colnames(phe)
gplots::balloonplot(
  table(phe$group,phe$orig.ident)
)

###

###首先看看比例饼图
head(sce@meta.data)
table(sce$celltype)
mynames <-   table(sce$celltype) %>% names()
myratio <-  table(sce$celltype) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "celltype",col = mycolors)


### 5.1 批量组合不同变量（分组，样品）看细胞类型比例  ---- 
cal_table = function(x,y,prefix ){
  # x = phe$orig.ident
  # y = phe$celltype
  library(sur)
  library(reshape2)
  tbl =  table(x,y)
  pdf(paste0(prefix,'-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( tbl )
  dev.off() 
  df = dcast(as.data.frame(tbl),x~y)
  head(df)
  write.csv(  df ,paste0(prefix,'-table.csv'))
  
  # ptb = round(sur::percent.table(x,y),2)
  ptb = round(100*tbl/rowSums(df[,-1]),2)
  
  pdf(paste0(prefix,'-percent-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( ptb )
  dev.off()
  write.csv(  dcast(as.data.frame(ptb),x~y) ,paste0(prefix,'-percent-table.csv')) 
  
}
cal_table(phe$orig.ident,phe$celltype,prefix = 'celltype-vs-orig.ident')
cal_table(phe$group,phe$celltype,prefix = 'celltype-vs-group')

cal_table(phe$group,phe$seurat_clusters,prefix = 'seurat_clusters-vs-group')

cal_table(phe$group,phe$RNA_snn_res.0.8,prefix = 'RNA_snn_res.0.8-vs-group')


### 5.2 不同分组的细胞类型比例折线图  ---- 
x='celltype';y='group' 
plot_data <- data.frame(table(phe[, y ],
                              phe[, x ]))
head(plot_data)
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
colnames(plot_data) <- c("group","cluster","Freq","Total","Percentage") 
head(plot_data)
ggplot(plot_data,aes(group,Percentage,group= cluster,color =cluster))+geom_point(size=4)+
  geom_line(position = position_dodge(0.1),cex=2)+theme_bw()+theme_test(base_size = 30) 

ggsave("celltype-vs-group-percent-plot.pdf",width = 9,height = 7)

### 5.3 不同分组的细胞类型比例箱线图  ---- 
x = phe$orig.ident
y = phe$celltype
tbl =  table(x,y)
df = dcast(as.data.frame(tbl),x~y)
ptb = round(100*tbl/rowSums(df[,-1]),2)
df = as.data.frame(ptb)
colnames(df)=c('orig.ident','celltype','Percentage')
head(df)
gpinfo=unique(phe[,c('orig.ident','group')]);gpinfo
df=merge(df,gpinfo,by='orig.ident')
head(df)

p <- ggplot(df,aes(x=group,y=Percentage))+
  geom_boxplot(outlier.alpha=0, aes(fill=celltype))+facet_grid(~celltype,scales = "free")+
  geom_jitter(aes(x = group, y = Percentage,color=celltype))+ 
  scale_color_npg() +
  scale_fill_npg() +
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold",angle = 45, hjust=1, vjust=0.5,size = 14), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill = "steelblue3", linetype = "solid"),
        strip.text = element_text(face = "bold", color = "black",hjust = 0.5, size = 12,),
        plot.title=element_text(face="bold.italic",size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(face = "bold",size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p   
pro='Percentage-of-each--celltype-in-orig.ident-by-group'
ggsave(paste0(pro,"_boxplot.pdf"),width = 12,height = 6)



### 5.4 不同分组的细胞类型比例堆积柱状图  ---- 
library(tidyr)
library(reshape2)
source('../scRNA_scripts/mycolors.R')
colnames(phe)
table(phe$group)

tb=table(phe$group, phe$celltype)
head(tb)
library (gplots) 
balloonplot(tb)
bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
#write.csv(bar_per,file = "celltype_by_group_percent.csv")
colnames(bar_per)

library(ggthemes)
p1 = ggplot(bar_per, aes(x = percent, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=mycolors)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))

p1

###当然也可以把坐标换为细胞数量
p2 = ggplot(bar_per, aes(x = Freq, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Total cell number')+
  scale_fill_manual(values=mycolors)+
  theme_classic()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p2

p1+p2

ggsave(filename="prop.pdf",width = 14,height = 7)



#Nebulosa是一个基于核密度估计的R软件包，用于可视化单个细胞的数据。做密度图可视化的一个重要的作用是，可以筛选鉴定某些基因双阳的细胞群。可以看到，我们能够看到CD4阳性的细胞、CD8阳性的细胞，还能显示双阳的细胞在哪，可以说是对seurat的一个拓展，还是比较实用的。
#BiocManager::install("Nebulosa")
library(Nebulosa)

####CD4发挥免疫作用的一系列基因list，可以添加
plot_density(sce, "CD4")
plot_density(sce, "CD8A")
###不加reduction会默认tsne
plot_density(sce, "CD8A",reduction = 'umap')
####CD8发挥免疫作用的一系列基因list
#也可以用Featureplot，右上角可以黄色越深的地方两个基因叠加越多
FeaturePlot(sce, features = c('CD3D','CD8A'),
            cols = c("lightgrey", "green", "orange"),
            order = T,raster = T ,
            blend=T,blend.threshold=0)
ggsave('FeaturePlot_1.pdf',width = 12,height = 7)

#Featureplot还可以把三个基因画在同一个图中
# 提取umap坐标
umap_df <- as.data.frame(sce@reductions$umap@cell.embeddings)
umap_df$cluster <- as.factor(sce$celltype)
head(umap_df)
# 提取基因表达数据并与umap坐标合并
# gene_df <- as.data.frame(GetAssayData(object = sce, slot = "data")[c('S100A9','S100A8','CXCL8'), ])
# 下面的 sce@assays$RNA$counts 如果是在Seurat的v4里面需要修改为
# sce@assays$RNA@counts， 千万有留意版本问题
gene_df <-  sce@assays$RNA$counts[c('CD4','CD8A','CD3E'),]
gene_df <-  t(as.data.frame(gene_df))
library(ggnewscale)
##如果有运行错误，可能是电脑和服务器有什么差异（同样代码遇到过bug，暂未查原因
merged_df <- cbind(gene_df, umap_df)
head(merged_df)
colnames(merged_df)
source('../scRNA_scripts/mycolors.R')
plot = ggplot(merged_df, vars = c("umap_1", "umap_2", 
                                  'CD4','CD8A','CD3E'), 
              aes(x = umap_1, y = umap_2, colour = CD4)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CD8A), size=0.3, alpha=0.7) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CD3E), size=0.3, alpha=0.1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish)+
  theme_classic()

plot2 <- DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,label.box = T, cols = mycolors,repel = T) +
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))

plot + plot2
ggsave('FeaturePlot_2.pdf',width = 14,height = 7)


getwd()
setwd('../')


