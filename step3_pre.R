phe = sce.all.int@meta.data
GSE72056 = phe[grepl("GSE72056", rownames(phe)), ]
GSE72056$sample = sub("^GSE72056_", "", rownames(GSE72056))
anno1_new = as.data.frame(t(anno1))
anno1_new$anno = anno1_new$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`
anno1_new$anno <- ifelse(anno1_new$anno == 1, "T",
                         ifelse(anno1_new$anno == 2, "B",
                                ifelse(anno1_new$anno == 3, "Macro",
                                       ifelse(anno1_new$anno == 4, "Endo",
                                              ifelse(anno1_new$anno == 5, "CAF",
                                                     ifelse(anno1_new$anno == 6, "NK", anno1_new$anno))))))


head(anno1_new)
anno1_new$an = anno1_new$`malignant(1=no,2=yes,0=unresolved)`
table(anno1_new$an)
anno1_new$an <- ifelse(anno1_new$an == 2, "malignant",
                         ifelse(anno1_new$an == 1, "non-malignant",'unresolved'))
table(anno1_new$an)
anno1_new = anno1_new[rownames(anno1_new) %in% GSE72056$sample,]
identical(rownames(anno1_new),GSE72056$sample)
anno1_new$celltype = GSE72056$celltype
GSE72056$orig.ident = GSE72056$sample
GSE72056$anno = anno1_new$anno
GSE72056$an = anno1_new$an
GSE72056 = GSE72056[,-17]
GSE72056$treatment = 'unknown'
GSE72056$therapy = 'unknown'
GSE72056$response = 'unknown'



GSE115978 = phe[grepl("GSE115978", rownames(phe)), ]
GSE115978$sample = sub("^GSE115978_", "", rownames(GSE115978))
anno3_new = anno2[anno2$cells %in% GSE115978$sample,]
identical(anno2$cells,GSE115978$sample)
p = identical(anno2$cells,GSE115978$sample);p
if(!p) anno2 = anno2[match(GSE115978$sample,anno2$cells),]
#GSE115978$orig.ident = GSE115978$sample
GSE115978$anno = anno3_new$cell.types
GSE115978$an = 'unknown'
GSE115978$treatment = anno3_new$treatment.group
GSE115978 = GSE115978[,-17]
GSE115978$therapy = 'unknown'
GSE115978$response = 'unknown'



GSE120575 = phe[grepl("GSE120575", rownames(phe)), ]
GSE120575$sample = sub("^GSE120575_", "", rownames(GSE120575))
anno3_new = anno3[anno3$title %in% GSE120575$sample,]
identical(anno3$title,GSE120575$sample)
p = identical(anno3$title,GSE120575$sample);p
if(!p) anno3 = anno3[match(GSE120575$sample,anno3$title),]
#GSE120575$orig.ident = GSE120575$sample
GSE120575$anno = 'unknown'
GSE120575$an = 'unknown'
GSE120575$treatment = anno3_new$`characteristics: patinet ID (Pre=baseline; Post= on treatment)`
GSE120575$therapy = anno3_new$`characteristics: therapy`
GSE120575$response = anno3_new$`characteristics: response`
GSE120575 = GSE120575[,-17]

phe_new = rbind(GSE115978,GSE120575,GSE72056)

save(phe_new,file = 'phe_new.Rdata')

