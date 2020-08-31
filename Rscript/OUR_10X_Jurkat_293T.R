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
Rfundir              = "/home/sam/Documents/RareCellDetection/Rfunction/"     
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
automatic_minpts     = TRUE                                  


workdir              = "/home/sam/Documents/RareCellDetection/Proj/jurkat_two_species/"    
# where you put the data and results

setwd(workdir)
library(pbapply)
library(future.apply)
library(ROCR)
library(Rcpp)
#library(splatter)
library(rflann)
library(e1071)
library(svd)
library(RaceID)
#library(minerva)
library(FiRE)
library(ineq)
library(Seurat)
library(irlba)
source('../../Main/utils.R')
#sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')

# data <- splatSimulate(group.prob=c(0.99, 0.01),method='groups', verbose=F,
#                       batchCells=500,de.prob=c(0.4, 0.4), out.prob=0,
#                       de.facLoc=0.4, de.facScale=0.8, nGenes=5000)
# counts <- data@assays@data$counts
# simulate.data <- list(counts=counts, group=data$Group)
# save(simulate.data, file='DE_simulate_data.RData')

source(paste(Rfundir,"GiniClust2_packages.R",sep=""))
source(paste(Rfundir,"GiniClust2_functions.R",sep=""))


ranger_preprocess<-function(data_mat, ngenes_keep=1000, dataSave='./', optionToSave=F, minLibSize=0, verbose=T){
  
  ngenes_keep = 1000
  #write.csv(x = data_mat$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
  #l<-.normalize_by_umi(data_mat)
  l<-.normalize_by_umi_2(data_mat, dataSave, minLibSize, verbose)
  m_n<-l$m
  
  if (verbose){
    cat("Select variable Genes...\n")
  }
  
  df<- .get_variable_gene(m_n)
  gc()
  
  if (verbose){
    cat("Sort Top Genes...\n")
  }
  
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  
  if (verbose){
    cat("Cutoff Genes...\n")
  }
  
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  #system("rm genes",ignore.stderr = T)
  
  if (optionToSave){
    write.csv(features, file = paste(dataSave,"genes",sep=""), quote = F,row.names = F)
    write.csv(l$use_genes[features], file = paste(dataSave,"genes_used_all",sep=""), quote = F,row.names = F)
  }
  
  #genes = read.csv(file = "genes")
  #features = genes$x
  
  #Final Data
  m_n_68K<-m_n[,features]
  m_filt<-Matrix(log2(m_n_68K+1),sparse = T)
  
  if (verbose){
    cat(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2]))
  }
  #system("rm whole_matrix",ignore.stderr = T)
  if (optionToSave){
    writeMM(m_filt,file="whole_matrix")
  }
  
  list(preprocessedData=m_filt, selGenes=features)
}

# ---------------------------------------------
# normalize the gene barcode matrix by umi
# filter based on read count first
# ---------------------------------------------

.normalize_by_umi_2 <-function(x, dataSave, minLibSize, verbose){
  mat  = x$mat
  gene_symbols = x$gene_symbols
  
  #Filter Cells
  if (minLibSize > 0){
    
    keepCellIndex <- c()
    for (i in c(1:dim(mat)[1])){
      count = sum(mat[i,] > 0)
      if (count > minLibSize){
        keepCellIndex <- c(keepCellIndex, i)
      }
    }
    
    mat <- mat[keepCellIndex,]
    if (verbose){
      cat(paste("Dimensions of matrix after cell filtering : ",dim(mat),"\n"))
    }
    write.csv(keepCellIndex, file = paste(dataSave,"keepCellIndexes.csv",sep=""), quote = F,row.names = F)
  }
  
  #Filter Genes
  cs <- colSums(mat>2)
  x_use_genes <- which(cs > 3)
  
  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]
  if (verbose){
    cat("Dimensions os filtered Matrix:")
    cat(paste(dim(x_filt),"\n"))
  }
  
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=gene_symbols)
}


# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m){
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}


transform.gini <- function(){
  gini.res <- list()
  for(i in 1:length(Gini.result)){
    cnt <- table(Gini.result[[i]])
    cnt.top <- table(Gini.result[[i]][1:5])
    small <- names(cnt)[cnt <= sum(cnt) * 0.15] # 0.15 is for 100 cells or a little more
    small <- small[small != 'Singleton']
    small <- small[small %in% names(cnt.top)]
    if(length(small)>=1){
      gini.res[[i]] <- ifelse(Gini.result[[i]] %in% small, 1, 0)
      
    }else{
      gini.res[[i]] <- (rep(0, length(Gini.result[[i]])))
    }
  }
  return(gini.res)
}


transform.race <- function(){
  race.res <- list()
  for(i in 1:length(RaceID.result)){
    
    cnt <- table(RaceID.result[[i]])
    cnt.top <- table(RaceID.result[[i]][1:5])
    small <- names(cnt)[cnt <= sum(cnt) * 0.15] # 0.15 is for 100 cells or a little more
    small <- small[small %in% names(cnt.top)]
    if(length(small)>=1){
      race.res[[i]] <- ifelse(RaceID.result[[i]] %in% small, 1, 0)
      
    }else{
      race.res[[i]] <- (rep(0, length(RaceID.result[[i]])))
    }
    
  }
  return(race.res)
}




#
set.seed(2019)

RaceID.result <- list()
Gini.result <- list()
OUR.result <- list()
FiRE.result <- list()
jurkat.num.vec <- c(8, 16, 24, 32, 40, 48, 56, 65, 73, 82)
for(jurkat.num in jurkat.num.vec){
  #load('DE_simulate_data.RData')
  experment_id <- which(jurkat.num.vec %in% jurkat.num)
  exprimentID<-paste0('Jurkat_', jurkat.num, '_293_1580_cells.RData')
  load(paste0('./results/', exprimentID))
  
  data <- as.matrix(jurkat.293[rowMeans(jurkat.293) > 0,])
  zero.cnt <- apply(data, 1, function(x){length(x[x>0])})
  data <- data[zero.cnt>=2,]
  norm.data <- .normalize_by_umi(t(data), gene_symbols = dimnames(data)[[1]], minLibSize=0, verbose = F)
  data2 <- t(norm.data$m)
  data21 <- data2
  data2 <- log2(data2 + 1)
  #group <- simulate.data$group
  group <- ifelse(grepl('Jurkat', dimnames(data2)[[2]]), 'Group1', 'Group2')
  
  pbmc <- CreateSeuratObject(count = data21)
  pbmc <- NormalizeData(object = pbmc, verbose = F)
  
  ## Different Fano genes for clustering
  pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
  vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
  den <- density(vst)
  cut.cnt <- length(which(vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])))
  features.vst <- dimnames(data2)[[1]][order(vst, decreasing = T)[1:(cut.cnt)]]
  
  ##################################################################################################
  #tmp <- data2[dimnames(data2)[[1]] %in% union(features, features.vst),]
  #tmp <- data2[order(vst, decreasing = T)[1:(2*length(features.vst))],]
  tmp <- (data2[dimnames(data2)[[1]] %in% (features.vst),])
  #cell.sum <- colSums(tmp)
  #tmp <- t(t(tmp)/(cell.sum/mean(cell.sum)))
  #pca <- calcul.pca(t(tmp), 50)
  pca <- irlba(t((tmp)), nv=50) # More robust no error, contrast to calcul.pca
  pca$pca <-t(pca$d*t(pca$u))
  #pca <- calcul.pca(t(data[dimnames(data)[[1]] %in% features,]), 50)
  #pca <- prcomp(t(tmp))
  #pca$pca <- pca$x[,1:min(50, dim(pca$x)[2])]
  #pca <- calcul.pca(t(data[order(disp, decreasing = T)[1:2000],]), 50)
  #knn.res <- Neighbour(pca$pca, pca$pca, k=200)
  
  knn.res <- Neighbour(pca$pca, pca$pca, k=200)
  

  distance.diff <- knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE]
  
  diff.left <- distance.diff[, -1, drop = FALSE] - distance.diff[, -ncol(distance.diff), drop = FALSE]
  diff.both <- diff.left[, -ncol(diff.left), drop=FALSE] - diff.left[, -1, drop=FALSE]
  diff.both[,1] <- diff.both[,1] + distance.diff[,1]  # Very important due to distance variation to the first neighbor.
  
  v1.k <- matrix(NA, dim(data2)[2], 197)
  skew <- c()
  skew1 <- c()
  top.ave <- c()
  remain.ave <- c()
  for(j in 1:dim(diff.both)[2]){
    
    #v <- distance.diff[,j]
    #v <- pmax(distance.diff[,j], diff.both[,j-1])
    v <- diff.both[,j]
    v1 <- v
    for(m in 1:length(v)){
      #v1[m] <- (v[m] + sum(v[knn.res$indices[m,2:(1+ceiling(log2(j+1)))]]))/(1+length(2:(1+ceiling(log2(j+1)))))
      v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
    }
    v1.k[, j] <- (v1)
    
    v2 <- v1[order(v1, decreasing = T)[(j+2):length(v1)]]
    skew1 <- c(skew1, skewness(v2))
    #top.values <- v1[order(v1, decreasing = T)[1:(j)]]
    top.values <- v1[knn.res$indices[which.max(v1),1:(j+1)]]
    #v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(exp(sum(log(top.values[top.values>0]))/length(top.values)), (2)))
    v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
    skew <- c(skew, skewness(v2))
    
    v1.sort <- sort(v1, decreasing = T)
    top.ave <- c(top.ave, mean(top.values))
    remain.ave <- c(remain.ave, mean(v1.sort[j:length(v1)]))
  }
  
  ids <- which(skew > 2)
  col.mat <- matrix(0, length(ids), dim(tmp)[2])
  for(i in 1:length(ids)){
    top.cell <- which.max(v1.k[,(ids[i])])
    col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * mean(v1.k[, ids[i]])
  }
  
  id.max <- apply(col.mat, 2, which.max)
  max.val <- apply(col.mat, 2, max)
  id.max[max.val==0] <- 0
  cnt <- table(id.max)
  cnt <- cnt[names(cnt)!='0']
  id.max.match <- cnt[which(cnt == (ids[as.integer(names(cnt))] + 1))] - 1
  
  cls <- rep(0, dim(tmp)[2])
  for(id.match in id.max.match){
    cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
  }
  
  rare.cells <- list()
  for(id.match in id.max.match){
    rare.cells[[as.character(id.match)]] <- knn.res$indices[which.max(v1.k[,id.match]), 1:(id.match+1)]
  }
  
  print(str(rare.cells))
  plot(1:dim(diff.both)[2], skew, cex=0.3, col=ifelse(1:197 %in% (jurkat.num-1), 'red', 'black'))
  
}



