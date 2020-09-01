setwd('/home/sam/Documents/FBT/Single/')
library(Rcpp)
library(pbapply)
library(future.apply)
library(ROCR)
library(Rcpp)
library(splatter)
library(rflann)
library(e1071)
library(svd)
library(irlba)
#sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')
source('/home/sam/Documents/RareCellDetection/Main/utils.R')

data <- read.csv('ExprM_3_subtypes_filter.csv', header=T)
expr.mat <- as.matrix(data[,2:dim(data)[2]])
dimnames(expr.mat)[[1]] <- data[,1]

disp <- FastLogVMR(as(expr.mat, 'dgCMatrix'), F)

CD14.data <- expr.mat[order(disp, decreasing=T)[1:2000],(grepl('CD14', dimnames(expr.mat)[[2]]))]
CD19.data <- expr.mat[order(disp, decreasing=T)[1:2000],(grepl('CD19', dimnames(expr.mat)[[2]]))]

res.14 <- dbscan::hdbscan(dist(t(CD14.data)), minPts = 2)
res.19 <- dbscan::hdbscan(dist(t(CD19.data)), minPts = 2)

get.major <- function(cluster){
  cnt <- table(cluster)
  cnt <- cnt[names(cnt) != '0']
  return(names(cnt)[which.max(cnt)])
}

expr.data <- expr.mat[, c(dimnames(CD14.data)[[2]][res.14$cluster==get.major(res.14$cluster)], dimnames(CD19.data)[[2]][res.19$cluster==get.major(res.19$cluster)])]

CD14.cells <- which(grepl('CD14', dimnames(expr.data)[[2]]))
CD19.cells <- which(grepl('CD19', dimnames(expr.data)[[2]]))

set.seed(20200602)
CD19.cells.500 <- sample(CD19.cells, 500)

data1 <- expr.data[, c(CD14.cells, CD19.cells.500)]
group <- ifelse(grepl('CD14', dimnames(data1)[[2]]), 'Group1', 'Group2')

DE <- WilcoxDETest((data1), cells.1 = dimnames(data1)[[2]][group=='Group1'], 
                   cells.2 = dimnames(data1)[[2]][group=='Group2'], verbose = F)
logFC <- apply(data1, 1, function(x){log(mean(x[group=='Group1'])/mean(x[group=='Group2']))})
infinite.na.id <- which(is.infinite(logFC) | is.na(logFC))
data2 <- data1[-infinite.na.id,]
p_vals <- DE$p_val[-infinite.na.id]
logFC1 <- logFC[-infinite.na.id]

DE.down <- dimnames(data2)[[1]][logFC1 < -log2(2) & p.adjust(p_vals, 'fdr') < 0.05]
DE.up.total <- dimnames(data2)[[1]][p.adjust(p_vals, 'fdr') < 0.05 & logFC1 > log2(2)]
DE.up <- DE.up.total[order(logFC1[p.adjust(p_vals, 'fdr') < 0.05 & logFC1 > log2(2)], decreasing = T)[1:length(DE.down)]]

NDE <- dimnames(data2)[[1]][p_vals > 0.05]
data3 <- data2[c(DE.down, DE.up, NDE),]
dimnames(data3)[[1]] <- paste0(dimnames(data3)[[1]], c(rep('_DE_down', length(DE.down)), rep('_DE_up', length(DE.up)), rep('_NDE', length(NDE))))
save(data3, file='/home/sam/Documents/RareCellDetection/Proj/10X_DE_sensitivity/DE_NDE_genes.RData')


## 20 CD14 cells
CD14.num = 10
for (i in 1){
  set.seed(i)
  samplingRare<-c(sample(which(grepl('CD14', dimnames(data3)[[2]])), CD14.num))
  sampling<-c(samplingRare, which(grepl('CD19', dimnames(data3)[[2]])))
  CD14.20.cells <- data3[,sampling]
  save(CD14.20.cells, file=paste("/home/sam/Documents/RareCellDetection/Proj/10X_DE_sensitivity/CD14_10_cells_", i, "_ExprM.RData", sep=""))
}



## 5 CD14 cells
CD14.num = 5
for (i in 1){
  set.seed(i)
  samplingRare<-c(sample(which(grepl('CD14', dimnames(data3)[[2]])), CD14.num))
  sampling<-c(samplingRare, which(grepl('CD19', dimnames(data3)[[2]])))
  CD14.20.cells <- data3[,sampling]
  save(CD14.20.cells, file=paste("/home/sam/Documents/RareCellDetection/Proj/10X_DE_sensitivity/CD14_5_cells_", i, "_ExprM.RData", sep=""))
}


## 2 CD14 cells
CD14.num = 2
for (i in 1){
  set.seed(i)
  samplingRare<-c(sample(which(grepl('CD14', dimnames(data3)[[2]])), CD14.num))
  sampling<-c(samplingRare, which(grepl('CD19', dimnames(data3)[[2]])))
  CD14.20.cells <- data3[,sampling]
  save(CD14.20.cells, file=paste("/home/sam/Documents/RareCellDetection/Proj/10X_DE_sensitivity/CD14_2_cells_", i, "_ExprM.RData", sep=""))
}


