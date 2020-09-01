library(rflann)
library(Seurat)
library(ineq)
library(e1071)
library(irlba)
homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/Runtime/')  
# where you put the data and results

setwd(workdir)
source('../../Rscript/utils.R')

require('Matrix')
require('plyr')


find_elbow <- function(x, y){
  n <- length(x)
  firstPoint <- c(x[1], y[1])
  lineVec = c(x[n]-x[1], y[n]-y[1])
  lineVecNorm = lineVec/(sqrt(sum(lineVec^2)))
  vecFromFirst = cbind(x-x[1], y-y[1])
  scalaProd =rowSums(vecFromFirst * cbind(rep(lineVecNorm[1], n), rep(lineVecNorm[2], n)))
  vecFromFirstParallel = outer(scalaProd, lineVecNorm)
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2))
  idx = which.max(distToLine)
  return(x[idx])
}


load('../10X_full/data/PBMC_68K_normalized.RData')


cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 1:length(cell.num)){
  set.seed(cell.num[i] + i)
  data1 <- data4[,sample(1:dim(data4)[2], cell.num[i])]
  data1 <- data1[rowMeans(data1) > 0,]

  start.time <- Sys.time()

  res <- GapClust::GapClust(data1, k=200)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  save(time.taken, file=paste0('./Proj/runtime/OUR_runtime_', i, '.RData'))
  
  if(cell.num[i] == 68579){
    save(predictions, file=paste0('./Proj/runtime/OUR_68K.RData'))
  }
  rm(list=c('start.time', 'time.taken', 'end.time', 'predictions', 'data1', 'pbmc', 
            'tmp', 'den', 'rare.score', 'filtered', 'skew', 'ks', 'remain',
            'cliff.id', 'v1.k', 'v1', 'v2', 'knn.res', 'pca', 'features.vst', 'vst'))
  gc()
  
}


