if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA",force = TRUE)

library(ggplot2)
library(GEOquery)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

load("D:/SDQ/TNBC/NEW_data/NEW_data/TNBC_ov1.rda")
meta=exprs(TNBC_ov1[[11]])
TCGA=exprs(TNBC_ov1[[12]])
setwd("D:\\SDQ\\TNBC\\NEW_data\\NEW_data\\predict drug senstivity\\DataFiles\\Training Data\\")
GDSC2_Expr = readRDS(file='GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS(file = "GDSC2_Res.rds")
# testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
# colnames(testExpr)=paste0('test',colnames(testExpr))

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = TCGA,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

calcPhenotype(trainingExprData = GDSC2_Expr,
                      trainingPtype = GDSC2_Res,
                      testExprData = meta,
                      batchCorrect = 'eb',  #   "eb" for ComBat  
                      powerTransformPhenotype = TRUE,
                      removeLowVaryingGenes = 0.2,
                      minNumSamples = 10, 
                      printOutput = TRUE, 
                      removeLowVaringGenesFrom = 'rawData' )



