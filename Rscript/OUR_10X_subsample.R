library(FiRE)
library(Rcpp)
library(rflann)
library(Seurat)
library(ineq)
library(e1071)
library(irlba)
homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/10X_subsample/')  
# where you put the data and results

setwd(workdir)
source('../../Rscript/utils.R')
#sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')


###############################################################
OUR.result <- list()
for (i in 2:100){
  #k=ks[j]
  exprimentID<-paste("10X_rare_",i,"_CD14_cells",sep="")
  load(paste0('results/', exprimentID, '_ExprM.filter.RData'))
  
  data <- ExprM.normCounts.filter[rowMeans(ExprM.normCounts.filter) > 0,]
  #data <- log2(data + 1)
  #disp <- FastLogVMR(as(data, 'dgCMatrix'), F)
  
  library(RaceID)
  #data2 <- log2(data +1)
  data2 <- data
  res <- GapClust::GapClust(data2)
  
  predictions <- rep(0, dim(data2)[2])
  predictions[res$rare_cell_indices[[which.max(res$skewness[as.integer(names(res$rare_cell_indices))])]]] <- 1
  OUR.result[[i-1]] <- predictions
  plot(1:length(res$skewness), res$skewness, cex=0.3, col=ifelse(1:197 %in% (i-1), 'red', 'black'))
  
}
save(OUR.result, file='results/OURClustering.RData')

