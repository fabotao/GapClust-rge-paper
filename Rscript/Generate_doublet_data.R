workdir              = "/home/sam/Documents/RareCellDetection/Proj/Doublet/"    
# where you put the data and results

setwd(workdir)
library(pbapply)
library(future.apply)
library(ROCR)
library(Rcpp)
library(splatter)
library(rflann)
library(e1071)
library(svd)
library(RaceID)
library(minerva)
library(FiRE)
library(Seurat)
library(ineq)
source('../../Main/utils.R')
#sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')

df <- read.csv('ExprM_3_subtypes_filter.csv', header=T)
mat <- as.matrix(df[,2:dim(df)[2]])
dimnames(mat)[[1]] <- df$X

pbmc <- CreateSeuratObject(count = mat)
pbmc <- NormalizeData(object = pbmc, verbose=F)
pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(mat)[1], verbose = F)
vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)


CD56.id <- which(grepl('CD56', dimnames(mat)[[2]]))
CD19.id <- which(grepl('CD19', dimnames(mat)[[2]]))
CD14.id <- which(grepl('CD14', dimnames(mat)[[2]]))

set.seed(2020)
samplingRare<-c(sample(CD14.id,2))
set.seed(100)
seeds<-sample(c(1:10000),140)
#subsample to create two common cell types
numbersB<-c(1600,800,400,200,100,50,25)
numbersNK<-c(800,400,200,100,50,25,13)

for (q in 1:7){
  for (i in 1:20){
    set.seed(seeds[(q-1)*20+i])
    sampling<-c(sample(CD56.id,numbersNK[q]),sample(CD19.id,numbersB[q]))
    sampling<-c(samplingRare,sampling)
    ExprM.RawCounts.filter <- mat[, sampling]
    exprimentID<-paste("10X_rare",q,"_",i, '_ExprM.filter.RData',sep="")
    save(ExprM.RawCounts.filter, file=paste("data/", exprimentID, sep=""))
  }
}








