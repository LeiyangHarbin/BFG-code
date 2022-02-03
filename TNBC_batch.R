
load("D:/SDQ/TNBC/NEW_data/gsva_TNBC.rda")
load("D:/SDQ/TNBC/NEW_data/TNBC_ov1.rda")
cname=list()
batch=c()
sample_name=c()
for (i in 1:length(gsva_TNBC)) {
  nam=row.names(gsva_TNBC[[i]])#得到数据集的细胞类型
  cname[i]=list(nam)
  batch=c(batch,rep(names(gsva_TNBC[i]),times=length(gsva_TNBC[[i]][1,])))#得到每个样本所属的数据集
  sample_name=c(sample_name,colnames(gsva_TNBC[[i]]))
}
cellname=Reduce(intersect,cname)#取交集得到共同的细胞类型
#得到样本共同细胞类型的得分
gsva_mat1=gsva_TNBC[[1]][cellname,]
for (i in 2:length(gsva_TNBC)) {#对每个数据集提取共同细胞类型的得分，然后合并这些数据为一个数据集
  mat=gsva_TNBC[[i]][cellname,]
  gsva_mat1=cbind(gsva_mat1,mat)
}
gsva_mat=gsva_mat1[,sample_name]#列按样本名排序

#批次效应
library(sva)
combat_edata=ComBat(dat=gsva_mat, batch=batch, mod=NULL)#使用combat函数消除批次，batch为批次变量

save(combat_edata,file = "D:/SDQ/TNBC/NEW_data/combat_edata.rda")
