homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/10X_subsample/')  
# where you put the data and results

setwd(workdir)
library(data.table)
source('../../Rfunction/CellSIUS.R')
source('../../Rfunction/CellSIUS_final_cluster_assignment.R')
source('../../Rfunction/CellSIUS_GetResults.R')



CELLSIUS.result <- list()
for (i in 2:100){
  #k=ks[j]
  exprimentID<-paste("10X_rare_",i,"_CD14_cells",sep="")
  load(paste0('results/', exprimentID, '_ExprM.filter.RData'))
  
  data <- ExprM.normCounts.filter[rowMeans(ExprM.normCounts.filter) > 0,]
  km <- kmeans(t(data), centers=2)
  res <- CellSIUS(data, group_id=km$cluster, min_n_cells = 2)
  #CellSIUS_GetResults(res)
  predictions <- CellSIUS_final_cluster_assignment(CellSIUS.out=res, group_id=km$cluster, min_n_genes = 2)
  
  CELLSIUS.result[[i-1]] <- predictions
  
}
save(CELLSIUS.result, file='results/CELLSIUSClustering.RData')

