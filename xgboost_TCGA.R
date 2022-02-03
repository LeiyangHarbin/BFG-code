####         https://github.com/tqchen/xgboost/blob/master/R-package/demo/tweedie_regression.R
###     http://codewithzhangyi.com/2018/06/01/XGBOOST%E4%BD%BF%E7%94%A8%E8%AF%B4%E6%98%8E/

library(xgboost)
library(survival)
library(survminer)
library(Matrix)
library(Ckmeans.1d.dp)
library(GEOquery)
library(pROC)

load("D:/SDQ/TNBC/NEW_data/NEW_data/TNBC_ov1.rda")
load("D:/SDQ/TNBC/NEW_data/NEW_data/生存分析/edata训练集生存时间和状态.rda")
load("D:/SDQ/TNBC/NEW_data/NEW_data/edata_0.01trainglm_riskScore.rda")
load("D:/SDQ/TNBC/NEW_data/NEW_data/combat_edata_NES.rda")
load("D:/SDQ/TNBC/NEW_data/NEW_data/edata_train_0.01sig_glmname.rda")

sample1=rownames(pData(TNBC_ov1[[12]]))
na1=intersect(sample1,rownames(risk_score_train_glm))
risk_score_TCGA=cbind(risk=risk_score_train_glm[na1,1],samples=na1)
risk_ST_TCGA=merge(risk_score_TCGA,Sta_Time_train,by="samples")#合并训练集的生存数据和risk数据
rownames(risk_ST_TCGA)=risk_ST_TCGA[,1]
risk_ST_TCGA=risk_ST_TCGA[,-1]
risk_ST_TCGA[,1]=as.numeric(as.character(risk_ST_TCGA[,1]))
risk_ST_TCGA=na.omit(risk_ST_TCGA)
cutoff=surv_cutpoint(risk_ST_TCGA, time = "times", event = "status", variables= "risk",minprop = 0.3, progressbar = TRUE)
cut=cutoff$cutpoint$cutpoint
high=rownames(risk_ST_TCGA[which(risk_ST_TCGA[,1]>= cut),])
low=rownames(risk_ST_TCGA[which(risk_ST_TCGA[,1]< cut),])
NES=t(combat_edata_NES[sig_glmname,c(high,low)])
NES=cbind(NES,label=rep(c(1,0),times=c(length(high),length(low))))
set.seed(123)
samp=sample(1:nrow(NES),nrow(NES)*0.7,replace=F)#抽取70%为训练集
train_NES=NES[samp,]
train_NES1 <- Matrix(train_NES[,-11],sparse=T) # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
test_NES=NES[-samp,]
test_NES1 <- Matrix(test_NES[,-11],sparse=T) # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵

#训练模型
dtrain <- xgb.DMatrix(data =train_NES1, label = train_NES[,11]) # 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
param <- list(max_depth=6, eta=0.3,  objective='binary:logistic') # 定义模型参数,基本上默认的参数就行，根据需要改变objective
#bst_train = xgb.train( data = dtrain, params = param,nrounds = 100, nthread = 2) # 构造xgboost模型,nround是迭代次数
# tree_train<- xgb.dump(bst_train,with_stats = T) # 显示计算过程，查看树结构
# tree_train
# xgb.plot.tree(model=bst_train,trees = 0)#查看第一个树
# names <- dimnames(data.matrix(train_NES[,1:175]))[[2]] # 获取特征的真实名称
# #计算特征重要性，对特征进行重要性排序
# importance_train<- xgb.importance(feature_names = names, model = bst_train)  
# head(importance_train)
# xgb.ggplot.importance(importance_train,top_n = 5)
# xgb.plot.importance(importance_matrix =importance_train,top_n = 5)
# xgb.plot.shap(data=train_NES,model=bst_train, top_n = 5,n_col = 2, pch = 16, pch_NA = 17)


#交叉验证,找到最佳的迭代次数
set.seed(123)
cv_res <- xgb.cv(params = param, data = dtrain,nrounds = 100,early_stopping_rounds=10,nfold=10, 
                 metrics=list("rmse","auc"),prediction = TRUE)
cv_nround=cv_res$best_iteration#最佳的迭代次数
bst = xgboost(data =train_NES1, label = train_NES[,11], max_depth = 6,
              eta = 0.3, nthread = 2, nrounds = cv_nround, objective='binary:logistic')#构造模型
names <- colnames(train_NES1) # 获取特征的真实名称
#计算特征重要性，对特征进行重要性排序
importance_train<- xgb.importance(feature_names = names, model = bst)   # importance 就是 信息增益
head(importance_train)
p1=xgb.ggplot.importance(importance_train,top_n = 10,n_clusters = 1)
p2<-p1#+ggplot2::ylab("Frequency")
p3<-p2+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=11,family="serif",color="black",face="bold"),#y轴刻度内容调整
             axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
             legend.title=element_text(size=11,family="serif",color="black",face="bold"),
             legend.text=element_text(size=9,family="serif",color="black",face="bold"),
             legend.position = "none"
)
p4<-p3+theme(panel.background = element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
             panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
             panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p5<-p4+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p5
ggsave("D:/SDQ/TNBC/NEW_data/NEW_data/建立预测模型/xgboost_TCGA-import.pdf",width = 7,height = 7)
#xgb.plot.importance(importance_matrix =importance_train,top_n = 5)
pdf(file = "D:/SDQ/TNBC/NEW_data/NEW_data/建立预测模型/xgboost_TCGA-shap.pdf",width = 9,height = 5)
par(family="serif",ps="15")
xgb.plot.shap(data=train_NES1,model=bst, top_n = 5, pch = 17,n_col = 5,pch_NA = 19,col="palevioletred1",col_loess="cyan2")
dev.off()

#计算混淆矩阵
dtest <- xgb.DMatrix(data =test_NES1, label = test_NES[,11])
pred=predict(bst,newdata = dtest)

#绘制ROC曲线
xgboost_roc <- roc(test_NES[,11],round(pred))
plot.roc(xgboost_roc, print.auc=T, print.auc.y=0.5,col="#1c61b6", percent=TRUE,main='ROC curve')
pdf(file = "D:/SDQ/TNBC/NEW_data/NEW_data/建立预测模型/xgboost_TCGAroc.pdf",width = 7,height = 7)
par(family="serif",ps="15")
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("white", "white"), 
     max.auc.polygon=F,auc.polygon.col="lightskyblue", 
     print.thres=TRUE,main='TCGA cohort')
dev.off()


#计算测试集精确度、召回率
xgb.result <- table(test_NES[,11],round(pred))
(xgb.result[1,1]+xgb.result[2,2])/length(test_NES[,11])#准确率
xgb.result[2,2]/(xgb.result[1,2]+xgb.result[2,2]) # 精确度
xgb.result[2,2]/(xgb.result[2,1]+xgb.result[2,2]) # 召回率


