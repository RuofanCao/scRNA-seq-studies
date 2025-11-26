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
dir.create('step2-paper-figures/')
setwd('step2-paper-figures/') 

## 1. 载入 ---- 
sce.all = readRDS('../3-Celltype/sce_celltype.rds')
sp='human' 
#load('../phe.Rdata') 
phe = sce.all@meta.data
identical(rownames(phe) , colnames(sce.all))
sce.all@meta.data = phe
colnames(phe)
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

## 2. 检查内置分群 ----
source('../scRNA_scripts/mycolors.R')
sce = sce.all
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),paletteer_d("awtools::a_palette"),paletteer_d("awtools::mpalette"))
p2_umap=DimPlot(sce, reduction = "umap",
                # split.by = "orig.ident",
                group.by = "celltype",label = T,label.size = 5,cols = color)
p2_umap
ggsave('umap_by_celltype_without_label.pdf',width = 9,height = 7)


p2_umap=DimPlot(sce, reduction = "umap",
                # split.by = "orig.ident",
                group.by = "celltype",repel = T,cols = mycolors,
                label = T,label.box = T)
p2_umap
ggsave('umap_by_celltype_with_label.pdf',width = 9,height = 7)


####构建左下角坐标轴umap图
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)
axes <- result$axes
label <- result$label
umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = "celltype",label = T,label.box = T) +
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
ggsave('umap_by_celltype_with_label_Bottom_left_axis.pdf',width = 9,height = 7)



th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  
# ord = c( "cycle" ,  'endo' 
#          , 'Microglial',"Mac" ,"Mono","NK" ,'Tcells','Bcells')
selected_genes <-  c(    'MKI67' , 'TOP2A',  
                      'PECAM1', 'VWF',  ## endo 
                     'CD68', 'CD163', 'CD14',  'C1QA',  'C1QB', # mac  # myeloids  
                     'S100A9', 'S100A8', 'MMP19',# monocyte
                     "GNLY","GZMA","NKG7","GZMK", # NK
                     'CD3D', 'CD3E', 'CD4','CD8A',  # Tcells 
                     'CD19', 'CD79A', 'MS4A1' , # Bcells 
                     'IGHG1', 'MZB1', 'SDC1',
                     'LUM' , 'FGF7', 'MME'
                   
) 
selected_genes <-  unique(str_to_upper(selected_genes))
p <- DotPlot(sce, features = selected_genes,
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="check_marker_by_celltype.pdf",width = 9,height = 7)


## 3. 各式各样的umap  ---- 
colnames(sce@meta.data)
table(sce$orig.ident)
library(stringr)
phe = sce@meta.data
table(phe$orig.ident)
colnames(phe)
####需要根据每篇文章中的自定义临床信息

#sce@meta.data = phe 
###重要：这里的group需根据自己想要的临床信息列设置
sce$group = phe$type
table(sce$group)
#注意颜色的修改，看个人喜欢什么色系
plot1 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,repel = T)+  
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
plot2 <- DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,label.box = T, cols = color,repel = T) +
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
plot3 <- DimPlot(sce, reduction = "umap", group.by = "group",label = T,cols = mycolors,repel = T) +
  theme(plot.title = element_text(size=12),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10))

##如果多个临床信息一起，可以循环
colnames(phe)
#c('tissue','sex','age','treatment')
if(F){
  
  p_list = lapply(c('sample','type'), function(x){
    DimPlot(sce, reduction = "umap", group.by = x ,
            label = T,cols = sample(color,6),repel = T) +
      theme(plot.title = element_text(size=12),
            axis.title = element_text(size=10),
            axis.text = element_text(size=10),
            legend.text = element_text(size=10))
  })
  p_list[[1]]
  
  fig1 = (p_list[[1]]+p_list[[2]])
  #fig1 = (p_list[[1]]+p_list[[2]])/(p_list[[3]]+p_list[[4]])
}
fig1
fig1 <- plot1|plot2|plot3
fig1
ggsave(fig1,filename = "fig1.pdf",width = 22,height = 7)
 
plot4 <- DimPlot(sce, reduction = "umap", group.by = "celltype",split.by = "group",
                 label = F,cols = mycolors) + theme(plot.title = element_text(size=12),
                                                 axis.title = element_text(size=10),
                                                 axis.text = element_text(size=10))

plot4
ggsave(plot4,filename = "fig2.pdf",width = 10,height = 7)


fig2 <- (plot1|plot2|plot3)/ plot4 + plot_layout(heights = c(1:2))
fig2
ggsave(fig2,filename = "fig1-2.pdf",width = 22,height = 11)


plot5 <- DimPlot(sce, reduction = "umap", group.by = "celltype",
                 split.by = "orig.ident",ncol = 5,label = F,cols = color) 
plot5
ggsave(plot5,filename = "fig2-2.pdf",units = "cm",width = 22,height = 40)



## 4. 已知的基因  ----
pro = 'fig3'
library(dplyr)  
library(paletteer)
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))

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
table(Idents(sce))
p_all_markers=DotPlot(sce, 
                      group.by  = "celltype",
                      features = genes_to_check,
                      scale = T,assay='RNA' ) + coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggplot2::ggsave(paste0(pro,'_all_genes_to_check-DotPlot.pdf'),
                height = 15,width = 8)

 
p1 <- VlnPlot(sce,features =  selected_genes,
              group.by  = "celltype",
              flip = T,stack = T,cols = color)
p1 + NoLegend()
ggplot2::ggsave(paste0(pro,'_genes_to_check-VlnPlot-heatmap.pdf'),height = 8,width = 7)

sce.Scale <- ScaleData(subset(sce,downsample=100),features =  selected_genes  )  
DoHeatmap(sce.Scale,
          features =  selected_genes,
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(paste0(pro,filename = "_genes_to_check-pheatmap.pdf") ,width = 12,height = 7)


## 4. COSG的基因  ---- 
table(Idents(sce))
table(Idents(sce.all))

###比findallmarker更加快速
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
  
  
  ## Top5 genes dotplot
  library(dplyr)  
  top5 <- unique(as.character(apply(marker_cosg$names,2,head,5)))
  top5
  # width <-0.006*dim(sce)[2];width
  # height <- 0.25*length(top_10)+4.5;height 
  top5_dotplot <- DotPlot(sce, features = top5 )+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  top5_dotplot
  ggsave('markers_top5_dotplot.pdf',width = 11,height = 7)
  
  
}
 


getwd()
pro = 'cosg_celltype_'
load(file = paste0(pro,'_marker_cosg.Rdata'))

top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10))) 
sce.Scale <- ScaleData(subset(sce,downsample=100),features =  top_10  )  
table(sce.Scale$celltype)
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
 
table(Idents(sce.Scale)) 

ll =  as.list(as.data.frame(apply(marker_cosg$names,2,head,10)))
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
ll
DoHeatmap(sce.Scale,
          features = unlist(ll),
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(filename = "fig4-top10-marker_pheatmap.pdf",width = 13,height = 8)



### 4.1 . 针对top基因绘制样本分类的小提琴图 ---- 
ll
##记得改group
sce$group = phe$type
table(sce$group )
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
p1 <- VlnPlot(sce ,features =  unique(unlist(ll)),flip = T,
              stack = T,
              group.by = "celltype",
              split.by = "celltype")+
  scale_color_npg()

p1
ggsave(p1,filename = "fig4_Vlnplot_top10.pdf",width = 8,height = 25)




## 5 看比例  ---- 
phe=sce@meta.data
head(phe) 
table(phe$group,phe$orig.ident)


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

cal_table(phe$group,phe$RNA_snn_res.0.1,prefix = 'RNA_snn_res.0.1-vs-group')


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
ggsave(paste0(pro,"_boxplot.pdf"),width = 12,height = 8)



### 5.4 不同分组的细胞类型比例堆积柱状图  ---- 
library(tidyr)
library(reshape2)
colnames(phe)
table(phe$patient)

##增加样本来源
phe$type1 = phe$orig.ident
table(phe$type1)
phe$type1[str_detect(phe$type1,"500")] <- "lung1"
phe$type1[str_detect(phe$type1,"532")] <- "lung2"
phe$type1[str_detect(phe$type1,"XX")] <- "lung3"
phe$type1[str_detect(phe$type1,"506")] <- "ovary"
phe$type1[str_detect(phe$type1,"516")] <- "kidney"
table(phe$type1)
table(phe$type1,phe$group)

##后面画分组柱状图会用到
phe$Procedure1 <- paste(phe$type1,phe$group,sep = '_')
table(phe$Procedure1)
sce@meta.data = phe

tb=table(sce$group, sce$celltype)
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
col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#BBBE00","#41BED1")
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
  scale_fill_manual(values=col)+
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
  scale_fill_manual(values=col)+
  theme_classic()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p2

p1+p2

ggsave(filename="prop.pdf",width = 14,height = 7)



###分组箱线图
###根据需要换列Procedure1
tb=table(sce$Procedure1, sce$celltype)
head(tb)
library (gplots) 
balloonplot(tb)
bar_data <- as.data.frame(tb)
library(stringr)
bar_data$group =  str_split( bar_data$Var1,'_',simplify = T)[,1]
#bar_data$group = factor(bar_data$group,levels = c("Resection",'EUS-Biopsy',"Liver-Biopsy"))

bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 

mycolors
#"#E64A35" "#4DBBD4" "#01A187" "#6BD66B" "#3C5588" "#F29F80"
library(ggpubr)
ggplot(bar_per, aes(x = group, y = percent))+ 
  labs(y="Cell Proportion", x =  NULL)+  
  geom_boxplot(aes(fill = Var2), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0)+ 
  scale_fill_manual(values = c( "#E64A35", "#4DBBD4", "#6BD66B", "#3C5588" ,"#F29F80")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) + 
  guides(fill=guide_legend(title = "Tissue"))+ 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T) 
ggsave(filename="celltype_group_prop.pdf",width = 13,height = 7)


###分组条形图
library(ggplot2)
library(ggpubr)
library(dplyr)
unique_types <- unique(phe$type1) 
dir.create("plots")
for (type in unique_types) {
  choose_cellIds <- rownames(phe)[phe$type1 == type]
  sce1 <- sce[, choose_cellIds]
  
  tb <- table(sce1$type, sce1$celltype)
  bar_data <- as.data.frame(tb)
  
  bar_per <- bar_data %>%
    group_by(Var1) %>%
    mutate(sum = sum(Freq)) %>%
    mutate(percent = Freq / sum)
  
  p <- ggplot(bar_per, aes(x = Var2, y = percent, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
    theme_bw() +
    labs(title = paste("", type), y = "Percentage", x = "Celltype") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.y = element_text(size = 12, color = "black")) +
    scale_fill_manual(values = c("#3176B7", "#F78000"))
  
  ggsave(filename = paste0("plots/", type, "_plot.pdf"), plot = p, width = 10, height = 7)
}



#saveRDS(sce, "sce_celltype.rds")
save(phe,file = 'phe.Rdata')

setwd('../')

