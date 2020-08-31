library(rflann)
library(Seurat)
library(ineq)
library(e1071)
library(irlba)
workdir              = "/lustre/home/acct-clsyzs/clsyzs/btfa/Runtime"    
# where you put the data and results

setwd(workdir)

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


load('./Proj/10X_full/data/PBMC_68K_normalized.RData')


cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 1:length(cell.num)){
  set.seed(cell.num[i] + i)
  data1 <- data4[,sample(1:dim(data4)[2], cell.num[i])]
  data1 <- data1[rowMeans(data1) > 0,]

  start.time <- Sys.time()

  pbmc <- CreateSeuratObject(count = data1)
  pbmc <- NormalizeData(object = pbmc, verbose = F)
  
  ## Different Fano genes for clustering
  pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data1)[1], verbose = F)
  vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
  den <- density(vst)
  features.vst <- dimnames(data1)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
  
  ##################################################################################################
  tmp <- data1[dimnames(data1)[[1]] %in% (features.vst),]
  #pca <- calcul.pca(t(log2(tmp+1)), 50)
  pca <- irlba(t(log2(tmp+1)), nv=50) # More robust no error, contrast to calcul.pca
  pca$pca <-t(pca$d*t(pca$u))
  
  knn.res <- Neighbour(pca$pca, pca$pca, k=200)
  
  #filtered.all <- rep(0, dim(data2)[2])
  v1.k <- matrix(NA, dim(data2)[2], 199)
  #ginis <- c()
  skew <- c()
  for(j in 2:200){
    if(j==2){
      v <- (knn.res$distances[,j]) - knn.res$distances[,j-1]
    }else{
      v <- (knn.res$distances[,j]) - knn.res$distances[,j-1]
    }
    v1 <- v
    for(m in 1:length(v)){
      v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
    }
    v1.k[, j-1] <- v1
    
    #gini <- ineq(v1, 'gini')
    #filtered.all<- filtered.all  +   gini * v1
    #ginis <- c(ginis, gini)
    
    v2 <- v1
    v2 <- v1[order(v1, decreasing = T)[(j+1):length(v1)]]
    v2 <- c(v2, rep(mean(v1[order(v1, decreasing = T)[1:j]]), 2))
    skew <- c(skew, skewness(v2))
  }

  cliff.id <- which.max(skew)
  
  ks <- which(skew[1:cliff.id] >= mean(skew[1:cliff.id]))
  ks <- ks[ks > (floor(cliff.id/2))-1]
  rare.score <- apply(v1.k[,ks, drop=F], 1, sum)
  filtered <- order(rare.score, decreasing = T)[1:((cliff.id)*3)]  
  
  den <- density(rare.score[filtered])
  remain <- order(rare.score, decreasing = T)[1:cliff.id]
  
  predictions <- rep(0, dim(data1)[2])
  predictions[remain] <- 1

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


