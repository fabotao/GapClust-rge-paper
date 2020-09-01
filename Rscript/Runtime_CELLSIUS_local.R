homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/Runtime/')  
# where you put the data and results

setwd(workdir)
source('../../Rscript/utils.R')

source('../../Rfunction/CellSIUS.R')
source('../../Rfunction/CellSIUS_final_cluster_assignment.R')
source('../../Rfunction/CellSIUS_GetResults.R')

load('../10X_full/data/PBMC_68K_normalized.RData')

df <- read.table('../10X_full/data/68k_pbmc_barcodes_annotation.tsv', header=T, sep='\t', stringsAsFactors = F)
cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 1:length(cell.num)){
  #k=ks[j]
  set.seed(cell.num[i] + i)
  ids <- sample(1:dim(data4)[2], cell.num[i])
  data1 <- data4[,ids]
  data1 <- data1[rowMeans(data1) > 0,]
  start.time <- Sys.time()

  res <- CellSIUS(data1, group_id=as.integer(factor(df$celltype[ids])), min_n_cells = 2)
  #CellSIUS_GetResults(res)
  predictions <- CellSIUS_final_cluster_assignment(CellSIUS.out=res, group_id=km$cluster, min_n_genes = 2)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  save(time.taken, file=paste0('CELLSIUS_runtime_', i, '.RData'))
  
  if(i == 68579){
    save(predictions, file=paste0('CELLSIUS_68K.RData'))
  }
  rm(list=c('predictions', 'data1', 'end.time', 'start.time', 'time.taken', 'res', 'ids'))
  gc()
  
  
}


