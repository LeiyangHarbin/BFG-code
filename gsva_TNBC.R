library(GEOquery)
library(GSVA)

load("D:/SDQ/TNBC/NEW_data/TNBC_ov1.rda")
load("D:/SDQ/TNBC/7-计算risk-score/celltypegene.rda")
#对TNBC计算每个数据集的ssGSEA
gsva_TNBC=list()
for (j in 1:length(TNBC_ov1)) {
  a1=exprs(TNBC_ov1[[j]])
  ce=rownames(a1)
  ce1=gsub(" ","_", ce)#将基因名字中的空格和-和.替换成"_"
  ce2=gsub("-","_",ce1)
  ce3=gsub("\\.","_",ce2)
  rownames(a1)=ce3
  a2=t(scale(t(a1)))
  gsva_es1<-gsva(a2,cellgene_type,method="ssgsea",abs.ranking=F)#计算免疫细胞的富集得分
  gsva_TNBC[j]=list(gsva_es1)
  names(gsva_TNBC)[j]=names(TNBC_ov1)[j]
}
save(gsva_TNBC,file = "D:/SDQ/TNBC/NEW_data/gsva_TNBC.rda")
