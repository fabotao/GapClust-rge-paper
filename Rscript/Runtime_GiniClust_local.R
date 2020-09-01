#fixed parameters for GiniClust2
minCellNum           = 3                                                # filtering, remove genes expressed in fewer than minCellNum cells
minGeneNum           = 2000                                             # filtering, remove cells expressed in fewer than minGeneNum genes
expressed_cutoff     = 1                                                # filtering, for raw counts
gini.bi              = 0                                                # fitting, default is 0, for qPCR data, set as 1. 
log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression   
log2.expr.cutoffh    = 20                                               # cutoff for range of gene expression 
Gini.pvalue_cutoff   = 0.0001                                           # fitting, Pvalue, control how many Gini genes chosen
Norm.Gini.cutoff     = 1                                                # fitting, NormGini, control how many Gini genes chosen, 1 means not used.
span                 = 0.9                                              # parameter for LOESS fitting
outlier_remove       = 0.75                                             # parameter for LOESS fitting
GeneList             = 1                                                # parameter for clustering, 1 means using pvalue, 0 means using HighNormGini
Gamma                = 0.9                                              # parameter for clustering
diff.cutoff          = 1                                                # MAST analysis, filter genes that don't have high a log2_foldchange to reduce gene num
lr.p_value_cutoff    = 1e-5                                             # MAST analysis, pvalue cutoff to identify differentially expressed genes
CountsForNormalized  = 100000                                           # if normalizing- by default not used
       
# where GiniClust2 R functions are stored

#dataset specific parameters:
MinPts               = 3                                                # parameter for DBSCAN
eps                  = 0.45                                             # parameter for DBSCAN
mycols               = c("grey50","greenyellow","red","blue","black","green","orange","purple","yellow","navy","magenta")
# color setting for tSNE plot
perplexity_G         = 30                                               # parameter for Gini tSNE
perplexity_F         = 30                                               # parameter for Fano tSNE
max_iter_G           = 1000                                             # parameter for Gini tSNE
max_iter_F           = 1000                                             # parameter for Fano tSNE
ks                   = c(2,2,2,3,3,3,3)                                 # a range of k's for k-means for subsampled data depending on rarity: use k=2 for rarer, k=3 for more common
gap_statistic        = FALSE                                            # whether the gap statistic should be used to determine k
K.max                = 10                                               # if using the gap statistic, highest k that should be considered
automatic_eps        = TRUE                                             # whether to determine eps using KNN- for consistency we use the same eps as full data set here
automatic_minpts     = TRUE                                             # whether to determine MinPts based on the size of the data set                                          
homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/Runtime/')  
# where you put the data and results

setwd(workdir)
Rfundir              = "../../Rfunction/"  
#dir.create(file.path(workdir, "results"), showWarnings = FALSE) #folder to save results
#dir.create(file.path(workdir, "figures"), showWarnings = FALSE) #folder to save figures
#load packages and functions

source(paste(Rfundir,"GiniClust2_packages.R",sep=""))
source(paste(Rfundir,"GiniClust2_functions.R",sep=""))

#generate 140 data sets, with different proportions
#source(paste(Rfundir,"Generate_10X_datasets.R",sep=""))

#for each of 99 data sets, run GiniClust2
#for each, plot a barplot comparing the reference and the GiniClust2 result
load('../10X_full/data/PBMC_68K_normalized.RData')

cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 2:100){
  set.seed(cell.num[i] + i)
  ids <- sample(1:dim(data4)[2], cell.num[i])
  data1 <- data4[,ids]
  data1 <- data1[rowMeans(data1) > 0,]
  start.time <- Sys.time()
  
  ExprM.normCounts.filter <- data1
  source(paste(Rfundir,"GiniClust2_fitting.R",sep=""))
  ExprM.RawCounts<-ExprM.RawCounts.filter
  source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep=""))

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  save(time.taken, file=paste0('Gini_runtime_', i, '.RData'))
  
  if(i == 68579){
    save(finalcls, file=paste0('RaceID_68K.RData'))
  }
  rm(list=c('data1', 'time.taken', 'end.time', 'start.time', 'ExprM.RawCounts', 'ExprM.normCounts.filter'))
  gc()

}

