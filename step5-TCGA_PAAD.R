rm(list=ls())
###根据TCGA-PAAD数据区分正常上皮和恶性上皮
dir.create("5-TCGA_PAAD")
setwd('5-TCGA_PAAD/')
source('../scRNA_scripts/mycolors.R')
library(tidyverse)
library(data.table) 
proj = "TCGA-PAAD"
dir.create("input")

##由于网络限制，这里我下载好了数据存放在input文件夹
##如果需要自己下载，把F改为T
chooseBioCmirror()
chooseCRANmirror()
if(T){
  download.file(url = paste0("https://gdc-hub.s3.us-east-1.amazonaws.com/download/",proj, ".star_counts.tsv.gz"),destfile = paste0("input/",proj,".htseq_counts.tsv.gz"))  ##表达数据
  download.file(url = paste0("https://gdc-hub.s3.us-east-1.amazonaws.com/download/",proj, ".clinical.tsv.gz"),destfile = paste0("input/",proj,".GDC_phenotype.tsv.gz")) ##临床数据
  download.file(url = paste0("https://gdc-hub.s3.us-east-1.amazonaws.com/download/",proj, ".survival.tsv.gz"),destfile = paste0("input/",proj,".survival.tsv")) ##生存数据
}

clinical = read.delim(paste0("input/",proj,".GDC_phenotype.tsv.gz"),fill = T,header = T,sep = "\t")
surv = read.delim(paste0("input/",proj,".survival.tsv"),header = T) 
head(surv) #生存数据os和os.time

### 1.处理表达矩阵和分组信息
#### 1.1 表达矩阵
#dat = read.table(paste0("input/",proj,".htseq_counts.tsv.gz"),check.names = F,row.names = 1,header = T)
dat <- data.table::fread(paste0("input/",proj,".htseq_counts.tsv.gz"),
                         data.table = F) #全部是symbol
dat[1:5,1:5]
dat[60650:60660,1:5]
a = dat[60550:60660,]
rownames(dat) = dat$Ensembl_ID
a = dat
a = a[,-1]
exp = a
#仅去除在所有样本里表达量都为零的基因
exp = exp[rowSums(exp)>0,]
nrow(exp)

#仅保留在一半以上样本里表达的基因
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
nrow(exp)


#表达矩阵行名ID转换
library(stringr)
head(rownames(exp))
library(AnnoProbe)
rownames(exp) = substr(rownames(exp), 1, 15)
##rownames(exp) = str_split(rownames(exp),"\\.",simplify = T)[,1];head(rownames(exp))
#annoGene(rownames(exp),ID_type = "ENSEMBL") 
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
#if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray",upgrade = FALSE,dependencies = TRUE)
library(tinyarray)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
exp[1:4,1:4]


#### 1.2 分组信息
#根据样本ID的第14-15位，给样本分组（tumor和normal）
table(str_sub(colnames(exp),14,15))  
table(str_sub(colnames(exp),14,16)) 
#01A：癌症组织
#01B：福尔马林浸泡样本
#02A：复发组织
#06A：转移组织
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')  
Group = factor(Group,levels = c("normal","tumor"))
table(Group)

# group = make_tcga_group(exp)
# table(group)
#保存数据
save(exp,Group,proj,clinical,surv,file = paste0(proj,".Rdata"))


#表达矩阵只需要tumor数据，不要normal，将其去掉，新表达矩阵数据命名为exprSet；
#临床信息需要进一步整理，成为生存分析需要的格式，新临床信息数据命名为meta。
#不同癌症的临床信息表格列名需要根据实际情况修改。
load('../5-TCGA_PAAD/TCGA-PAAD.Rdata')
library(stringr)

### 1.整理表达矩阵
#### 1.1去除normal样本
table(Group)
exprSet = exp[,Group=='tumor']
ncol(exp)
ncol(exprSet)

#### 1.2 基因过滤
#再次进行基因过滤.
#(1)标准1：至少要在50%的样本里表达量大于0（最低标准）。
k = apply(exprSet,1, function(x){sum(x>0)>0.5*ncol(exprSet)});table(k)
exprSet = exprSet[k,]
nrow(exprSet)

#(2)标准2：至少在一半以上样本里表达量＞10(其他数字也可，酌情调整)
k = apply(exprSet,1, function(x){sum(x>10)>0.5*ncol(exprSet)});table(k)
exprSet = exprSet[k,]
nrow(exprSet)

#### 1.3 使用logCPM或logTPM数据
exprSet[1:4,1:4]
exprSet=log2(edgeR::cpm(exprSet)+1)
exprSet[1:4,1:4]

#### 1.4 样本与病人
#有的病人会有两个或两个以上的肿瘤样本，就有重复。两种可行的办法：
#（1）以病人为中心，对表达矩阵的列按照病人ID去重复，每个病人只保留一个样本。
exprSet = exprSet[,sort(colnames(exprSet))]
k = !duplicated(str_sub(colnames(exprSet),1,12));table(k)
exprSet = exprSet[,k] 
ncol(exprSet)
#（2）以样本为中心，如果每个病人有多个样本则全部保留。(删掉上面这一段代码即可)

### 2.整理生存信息和临床信息
#### 2.1数据合并
library(dplyr)
head(surv) #生存数据os和os.time
colnames(clinical)
colnames(surv)
identical(surv$sample,clinical$sample)
meta = left_join(surv,clinical,by = c("sample"= "sample"))
nrow(meta)
length(unique(meta$sample))
meta = distinct(meta,sample,.keep_all = T)

#### 2.2 样本过滤
#去掉生存信息不全或者生存时间小于30天的样本，样本纳排标准不唯一，且差别很大
k1 = meta$OS.time >= 30;table(k1)
k2 = !(is.na(meta$OS.time)|is.na(meta$OS));table(k2)
meta = meta[k1&k2,]

#### 2.3 选列、简化列名
tmp = data.frame(colnames(meta))
meta = meta[,c(
  'sample',
  'OS',
  'OS.time',
  'treatment_id.treatments.diagnoses',
  'submitter_id.treatments.diagnoses' ,
  'treatment_type.treatments.diagnoses',
  'treatment_or_therapy.treatments.diagnoses' 
)]
colnames(meta)[1:3]=c('ID','event','time')
str(meta)

#### 2.4 简化、规范内容
#(1) 结局事件
#生存分析的输入数据里，要求结局事件必须用0和1表示，1表示阳性结局。
#xena的数据是整理好的，其他来源的需要自行检查和整理。
table(meta$event)

#(2) 生存时间
range(meta$time)
meta$time = meta$time/30
range(meta$time)


### 3.实现表达矩阵与临床信息的匹配 
rownames(meta) <- meta$ID
head(rownames(meta))
head(colnames(exprSet))
colnames(exprSet) = str_sub(colnames(exprSet),1,16)
s = intersect(rownames(meta),colnames(exprSet));length(s)
exprSet = exprSet[,s]
meta = meta[s,]
dim(exprSet)
dim(meta)
identical(rownames(meta),colnames(exprSet))
save(meta,exprSet,proj,file = paste0(proj,"_sur_model.Rdata"))


### 生存分析
proj = "TCGA-PAAD"
load(paste0(proj,"_sur_model.Rdata"))
exprSet[1:4,1:4]
str(meta)
###这里有一种方式是直接使用单细胞亚群里的差异分析得到每个细胞亚群的差异基因集，统一进行GSVA分析后，循环得到每个细胞亚群在TCGA里的生存情况
#根据gsva结果高低分组后批量生存分析----
sce=readRDS( "../2-harmony/sce.all_int.rds")
load('../4-cellchat/celltype_input_phe.Rdata')
sce@meta.data = Allphe
pro = 'cosg_celltype_'
library(COSG)
library(survival)
library(survminer)
library(ggstatsplot)
library(gplots)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
marker_cosg <- cosg(
  sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)
save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))

head(marker_cosg$names)
deg_list=as.list(marker_cosg$names)
names(deg_list)
deg_list
gs=lapply(deg_list,toupper)
geneset<-GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
},gs,names(gs)))
geneset

#3.run gsva----
X=as.matrix(exprSet)
gsva.res <- gsvaParam(X, geneset) 
es.max<-gsva(gsva.res)
es.max[1:4,1:4]
pheatmap(es.max)
es.max[1:4,1:4]
splots<-list()
g=1
phe = meta
for(i in names(deg_list)){
  #i=names(deg_list)[1]
  subset=paste0('cluster_',i)
  print(subset)
  v=as.numeric(es.max[i,])#每一个亚群表达量。
  sub_group<-ifelse(v<0,"low","high")#如果表达量小于0的话，就定义为low。gsva处理过表达量。0.几左右
  table(sub_group)
  phe$sub_group=sub_group
  #Fitsurvivalcurves
  require("survival")
  fit<-survfit(Surv(time,event)~sub_group,data=phe)
  library("survminer")
  survp<-ggsurvplot(fit,data=phe,
                    surv.median.line="hv",
                    pval=TRUE,
                    conf.int=TRUE,
                    risk.table=TRUE,
                    tables.height=0.2,
                    tables.theme=theme_cleantable(),
                    palette="jco",
                    ggtheme=theme_bw(),
                    title=subset)
  print(survp)
  splots[[g]]<-survp
  g=g+1
}

length(splots)
x1=ceiling(sqrt(length(splots)))
y1=x1

all_plot<-arrange_ggsurvplots(splots,
                              print=F,
                              ncol=x1,nrow=y1,
                              risk.table.height=0.3,
                              surv.plot.height=0.7)
#all_plot
x2=5*x1
y2=5*y1
prefix=''
pro=''
ggsave(all_plot,#path=prefix,
       filename=paste0(pro,'all_survival_plot.pdf'),
       width=x2,height=y2)


### 2.KM-plot
library(survival)
library(survminer)
library(readxl)
###上面是更加直观的取差异基因计算，当然如果你有更好的标记亚群的基因，可以直接单独评分，只看感兴趣的基因集和亚群
##NK_C1群的基因集
NK_sheet1 <- read_excel("../elife-92672-fig6-data1.xlsx")
NK_sheet1 = as.data.frame(NK_sheet1)
NK_sheet1 = as.data.frame(NK_sheet1[-1,]) 
colnames(NK_sheet1) = 'NK_gene'
library(tidyverse)
library(GSVA)
geneset = list(NK_sheet1$NK_gene)
X=as.matrix(exprSet)
gsva.res <- gsvaParam(X, geneset) 
expr_geneset <- gsva(gsva.res)
dim(expr_geneset)
identical(colnames(expr_geneset) ,rownames(phe))
phe$NK_score=expr_geneset[1,]
#这里用中位数为示例，实际根据个人需要
phe$NK_score<-ifelse(phe$NK_score<median(phe$NK_score),"low","high")
table(phe$NK_score)
fit<-survfit(Surv(time,event)~NK_score,data=phe)
p=ggsurvplot(fit, data = phe,
                 surv.median.line = "hv",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 palette = "jco",
                 ggtheme = theme_bw(),
                 title = 'TCGA-PAAD NK signature', newpage = FALSE)
p
pdf("NK_sur.pdf", width = 9, height = 7)
print(p, newpage = FALSE)
dev.off()

##不显著，看原文中可能也是使用的最佳cutoff
phe$NK_score=expr_geneset[1,]
sur.cut<-surv_cutpoint(phe,
                       time='time',
                       event='event',
                       variables='NK_score')
summary(sur.cut)
plot(sur.cut, "NK_score", palette = "npg")
sur.cat<-surv_categorize(sur.cut)
head(sur.cat)
sfit<-survfit(Surv(time,event)~NK_score,data=sur.cat)
survp2<-ggsurvplot(sfit,data=phe,
                   surv.median.line="hv",
                   pval=TRUE,
                   conf.int=TRUE,
                   risk.table=TRUE,
                   tables.height=0.2,
                   tables.theme=theme_cleantable(),
                   palette="jco",
                   ggtheme=theme_bw(),
                   title='TCGA-PAAD NK signature')
survp2
pdf("NK_cutoff_sur.pdf", width = 9, height = 7)
print(survp2, newpage = FALSE)
dev.off()

###是否治疗条件下的生存分析
colnames(phe)
colnames(clinical)
colnames(clinical)[grepl("treat|treatment", colnames(clinical), ignore.case = TRUE)]
table(clinical$prior_treatment.diagnoses)
table(phe$treatment_type.treatments.diagnoses)
table(phe$submitter_id.treatments.diagnoses)
table(phe$treatment_or_therapy.treatments.diagnoses)

data_str <- phe$treatment_or_therapy.treatments.diagnoses
# 去除中括号
data_str <- gsub("\\[|\\]", "", data_str)
# 使用strsplit分割每个字符串元素，并去除引号
split_data <- unlist(lapply(data_str, function(x) {
  x <- strsplit(x, ", ")[[1]]
  x <- gsub("'", "", x)
  return(x)
}))

result_df <- data.frame(matrix(split_data, ncol = 2, byrow = TRUE))
result_df
result_df$X1 <- gsub("not reported", "no", result_df$X1, fixed = TRUE)
result_df$X2 <- gsub("not reported", "no", result_df$X2, fixed = TRUE)


data <- phe$treatment_type.treatments.diagnoses
data <- gsub("\\[|\\]", "", data)
data <- gsub("'", "", data)
split_data1 <- strsplit(data, ", ")
result_df1 <- do.call(rbind, split_data1) %>% as.data.frame()
data = cbind(result_df,result_df1)
data$treat1 = paste(data$X1,data$V1,sep = '_')
data$treat2 = paste(data$X2,data$V3,sep = '_')
table(data$treat1)
table(data$treat2)
data$group <- ifelse(data$treat1 == 'yes_Radiation Therapy' | data$treat2 == 'yes_Radiation Therapy',"IR_yes","IR_no")
table(data$group)
phe$group = data$group


##基因生存分析
#将肿瘤样本分为两组，所有患者的靶基因表达根据中位数分别为高低
gs=c('NCAM1','CCL5')
dim(exprSet)
colnames(exprSet)
rownames(exprSet)
exprSet = as.data.frame(exprSet)
exprSet['NCAM1',]
exp['NCAM1',]
###果然前面过滤的太多了，只剩下七千多个基因竟然没有CD56
# splots <- lapply(gs, function(g){
#   phe$gene=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
#   # quantiles <- quantile(exprSet[g, ])
#   # phe$gene <- ifelse(exprSet[g, ] > quantiles[4], 'high', 'low')
#   sfit1=survfit(Surv(time, event)~gene, data=phe)
#   ggsurvplot(sfit1,pval =TRUE, data = phe,
#              title = g)
# }) 
# arrange_ggsurvplots(splots, print = TRUE,  
#                     ncol = 2, nrow = 2)

#只能单独看一个基因了
identical(rownames(phe),colnames(exprSet))
phe$CCL5=ifelse(as.numeric(exprSet['CCL5',]) >median(as.numeric(exprSet['CCL5',])),'high','low')
table(phe$CCL5)
sfit1=survfit(Surv(time, event)~CCL5, data=phe)
ggsurvplot(sfit1,pval =TRUE, data = phe,title = 'CCL5')
#单独看是没有差异的，虽然也没有根据最佳cutoff来选定
#下面直接根据两个条件来分组
phe$CCL5=ifelse(as.numeric(exprSet['CCL5',]) >median(as.numeric(exprSet['CCL5',])),'high','low')

colnames(phe)
colnames(phe)[12] = 'IR_treatment'
table(phe$IR_treatment)
###要注意group顺序，排列先后的问题
phe$NK_IR <- with(phe, ifelse(NK_score >= 0.7886529 & IR_treatment == "IR_yes", "NK high + RTx",
                  ifelse(NK_score >= 0.7886529, "NK high",
                         ifelse(NK_score < 0.7886529 & IR_treatment == "IR_yes", "NK low + RTx",'NK low'))))
table(phe$NK_IR)
sfit1=survfit(Surv(time, event)~NK_IR, data=phe)
p = survminer::ggsurvplot(sfit1,pval =TRUE,
                          data = phe,
                          legend = c(0.8,0.8),
                          title = 'TCGA-PAAD',
                          legend.title='Expression',
                          xlab="Time_month",
                          surv.median.line = 'hv',
                          risk.table = T,
                          legend.labs = c("NK high", "NK high + RTx","NK low", "NK low + RTx")
)
p
p2 = p$plot+
  theme(plot.title = element_text(hjust = 0.5))
p2
ggsave(p2,
       filename='NK_IR_sur.pdf',
       width=9,height=7)


setwd('../')



