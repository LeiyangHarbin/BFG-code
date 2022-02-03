#                                       #选后面六个数据集为训练集
library(GEOquery)
library(survminer)
library(survival)
library(survcomp)
library(ggplotify)
load("D:/SDQ/TNBC/NEW_data/NEW_data/combat_edata_NES.rda")
load("D:/SDQ/TNBC/NEW_data/NEW_data/TNBC_ov1.rda")
status=c()
times=c()
samples=c()
for (k in 7:length(TNBC_ov1)) {  #得到每个样本的生存时间和状态，训练集为7-12个数据集
  phe=pData(TNBC_ov1[[k]])
  tim=phe$Overall.survival.months
  sta=as.character(phe$Overall.survival)
  samp=row.names(phe)
  times=c(times,tim)
  status=c(status,sta)
  samples=c(samples,samp)
}
status[status == "alive"] <- 0 #将alive和dead转化成0,1
status[status == "dead"] <- 1
status=as.integer(status)
Sta_Time_train=data.frame(samples = samples,status=status,times=times)#将训练集样本的生存时间和状态写入数据框
Sta_Time_train=na.omit(Sta_Time_train)
Sta_Time_train=Sta_Time_train[-which(Sta_Time_train[,3]==0),]

save(Sta_Time_train,file = "D:/SDQ/TNBC/NEW_data/NEW_data/生存分析/edata训练集生存时间和状态.rda")

train_data=as.data.frame(t(combat_edata_NES[,as.character(Sta_Time_train$samples)]))
train_data1=cbind(train_data,samples=rownames(train_data))
train_data2=merge(train_data1,Sta_Time_train, by  = "samples")
rownames(train_data2)=train_data2[,1]
train_data=train_data2[-1]

ce=colnames(train_data)
ce1=gsub(" ","_", ce)#将细胞名字中的空格和-和.替换成"_"
ce2=gsub("-","_",ce1)
ce3=gsub("\\.","_",ce2)
colnames(train_data)=ce3

#对训练集进行单因素cox
train_cox_NES<-matrix(nrow=ncol(train_data)-2,ncol=6)
for (j in 1:(ncol(train_data)-2)) {#对去除批次效应的训练集ssGSEA得分计算cox
  Bcox<-coxph(Surv(times, status)~as.numeric(as.character(train_data[,j]))>median(as.numeric(as.character(train_data[,j]))),data=train_data)
  summcph<-summary(Bcox)
  train_cox_NES[j,1]<-summcph$conf.int[1]
  train_cox_NES[j,2]<-summcph$conf.int[3]
  train_cox_NES[j,3]<-summcph$conf.int[4]
  train_cox_NES[j,4]<-as.matrix(summcph$logtest)[3]
  train_cox_NES[j,5]<-as.matrix(summcph$sctest)[3]
  train_cox_NES[j,6]<-summcph$coefficients[5]
}
rownames(train_cox_NES)=colnames(train_data)[1:175]#矩阵行名为细胞类型
colnames(train_cox_NES)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")#列名

train_log=train_cox_NES[train_cox_NES[,5]<0.01,]#筛选训练集p值小于0.01的数据集
save(train_log,file = "D:/SDQ/TNBC/NEW_data/NEW_data/edata_train_0.01unicox.rda")



####glmnet包说明文档https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#cox
###https://cosx.org/2016/10/data-mining-1-lasso/
###P值小于0.01计算risk得分---------------------------------------------------------------------------------------------------
library(glmnet)
#LASSO cox分析缩减变量
opar=par(no.readonly = T)
par(family="serif",ps="15")

set.seed(1000)#设置随机种子
sig_train_data=as.matrix(train_data[,rownames(train_log)])#得到变量矩阵，必须是矩阵，否则下面计算报错
TS=cbind(time=train_data$times,status=train_data$status)#设置响应变量，cox分析用时间和状态

train_glm=glmnet(sig_train_data,TS,family = "cox")#列 %Dev 代表了由模型解释的残差的比例,它在 0 和 1 之间，越接近 1 说明模型的表现越好
train_glm

coe_glm=coef(train_glm, s=c(train_glm$lambda[13]))#由%Dev结果可以知道当lambda=0.035880时,%Dev 最接近1
act_index <- which(coe_glm != 0)#筛选基因
act_coe <- coe_glm[act_index]
sig_glmname_glm=row.names(coe_glm)[act_index]#得到缩减后变量
sig_glmname_glm #只做了一次拟合的结果，可能过拟合

pdf("D:/SDQ/TNBC/NEW_data/NEW_data/lasso_coef.pdf",width=6.0,height=6.0)
par(family="serif",ps="15")
plot(train_glm, xvar="lambda", label=F)
#abline(v=log(train_cvglm$lambda.min),lty=2)
dev.off()



#采用交叉验证（cross validation）拟合进而选取模型，对模型的性能有一个更准确的估计
train_cvglm=cv.glmnet(sig_train_data,TS,family = "cox")#交叉验证，计算最佳的lambda值
#基因筛选，采用coef函数即可，有相应参数的gene则被保留，采用λ使用的是lambda.min
coe.min=coef(train_cvglm,s=train_cvglm$lambda.min)
act_index <- which(coe.min != 0)
act_coe <- coe.min[act_index]
sig_glmname=row.names(coe.min)[act_index]#得到缩减后变量
sig_glmname

pdf("D:/SDQ/TNBC/NEW_data/NEW_data/lasso_dev.pdf",width=6.0,height=6.0)
par(family="serif",ps="15")
plot(train_cvglm,xlab="Log Lambda")
dev.off()

par(opar)

lasso1=as.ggplot(~plot(train_glm, xvar="lambda", label=F))
lasso2=as.ggplot(~plot(train_cvglm,xlab="Log Lambda"))
save(lasso1,lasso2,file = "D:/SDQ/TNBC/NEW_data/NEW_data/LASSO_plot.rda")

save(sig_glmname,file = "D:/SDQ/TNBC/NEW_data/NEW_data/edata_train_0.01sig_glmname.rda")

#计算risk值
HR=as.vector(train_cox_NES[sig_glmname,1])#得到HR
sig_train_NES=train_data[,sig_glmname]#得到显著细胞的ssGSEA得分

risk_score_train=c()
for (i in 1:nrow(sig_train_NES)) {
  rscore=0
  for (j in 1:ncol(sig_train_NES)) {
    rscore=rscore+log(HR[j])*sig_train_NES[i,j]
  }
  risk_score_train=c(risk_score_train,rscore)
}

risk_score_train_glm=as.matrix(risk_score_train)
rownames(risk_score_train_glm)=rownames(sig_train_NES)
colnames(risk_score_train_glm)="risk_score"
save(risk_score_train_glm,file = "D:/SDQ/TNBC/NEW_data/NEW_data/edata_0.01trainglm_riskScore.rda")




