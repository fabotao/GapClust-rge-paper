library(RaceID)
workdir              = "/lustre/home/acct-clsyzs/clsyzs/btfa/Runtime"     
# where you put the data and results

setwd(workdir)
load('/Proj/10X_full/data/PBMC_68K_normalized.RData')

df <- read.table('/Proj/10X_full/data/68k_pbmc_barcodes_annotation.tsv', header=T, sep='\t', stringsAsFactors = F)
cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 1:length(cell.num)){
  #k=ks[j]
  set.seed(cell.num[i] + i)
  ids <- sample(1:dim(data4)[2], cell.num[i])
  data1 <- data4[,ids]
  data1 <- data1[rowMeans(data1) > 0,]
  start.time <- Sys.time()
  sc <- SCseq(data1)
  sc <- filterdata(sc,mintotal=min(colSums(data1))-1)
  fdata <- getfdata(sc)
  sc <- compdist(sc,metric="pearson")
  sc <- clustexp(sc)
  sc <- clustexp(sc,cln=length(unique(df$celltype[ids])),sat=FALSE)
  sc <- findoutliers(sc)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  save(time.taken, file=paste0('RaceID_runtime_', i, '.RData'))
  
  finalcls <- (sc@cpart)
  print(table(finalcls))
  if(i == 68579){
    save(finalcls, file=paste0('RaceID_68K.RData'))
  }
  rm(list=c('sc', 'time.taken', 'end.time', 'start.time', 'data1', 'finalcls'))
  gc()
}

