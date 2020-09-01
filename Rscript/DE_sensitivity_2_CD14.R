#fixed parameters for GiniClust2
minCellNum           = 2                                                # filtering, remove genes expressed in fewer than minCellNum cells
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
MinPts               = 2                                                # parameter for DBSCAN
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


homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/10X_DE_sensitivity/')  
Rfundir              = "../../Rfunction/"  
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
source('../../Rscript/utils.R')

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
    cnt.top <- table(Gini.result[[i]][1:2])
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
    cnt.top <- table(RaceID.result[[i]][1:2])
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


#load('DE_simulate_data.RData')
load('../10X_DE_sensitivity/CD14_2_cells_1_ExprM.RData')
#data <- simulate.data$counts[rowMeans(simulate.data$counts) > 0,]
data <- CD14.20.cells[rowMeans(CD14.20.cells) > 0,]
norm.data <- .normalize_by_umi(t(data), gene_symbols = dimnames(data)[[1]], minLibSize=0, verbose = F)
data3 <- t(norm.data$m)
#group <- simulate.data$group
group <- ifelse(grepl('CD14', dimnames(data3)[[2]]), 'Group1', 'Group2')
#disp <- FastLogVMR(as(data3, 'dgCMatrix'), F)
#DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][group=='Group1'], 
#                   cells.2 = dimnames(data3)[[2]][group=='Group2'], verbose = F)
#logFC <- apply(data3, 1, function(x){log(mean(x[group=='Group1'])/mean(x[group=='Group2']))})

#mic <- mine(t(data3), y=as.integer(factor(group))) #not work
#DE.ind <- which(p.adjust(DE$p_val, 'fdr') < 0.05 & abs(logFC)>log2(5))
#DE.ind <- which(abs(logFC) > log2(5) & p.adjust(DE$p_val, 'fdr') < 0.05)
#NDE.ind <- which(DE$p_val >= 0.05)
#names(DE.ind) <- NULL
#names(NDE.ind) <- NULL
DE.ind <- which(grepl('_DE', dimnames(data3)[[1]]))
NDE.ind <- which(grepl('_NDE', dimnames(data3)[[1]]))


#
set.seed(2019)
#DE.range <- 1:length(DE.ind)
DE.range <- 1:length(DE.ind)
auc.ave <- c()
auc.std <- c()
auc.ave.our <- c()
auc.std.our <- c()
auc.ave.gini <- c()
auc.std.gini <- c()
auc.ave.15 <- c()
auc.std.15 <- c()
auc.ave.10 <- c()
auc.std.10 <- c()
auc.ave.05 <- c()
auc.std.05 <- c()


for(num in DE.range){
  RaceID.result <- list()
  Gini.result <- list()
  auc.our <- c()
  auc.fire.15 <- c()
  auc.fire.10 <- c()
  auc.fire.05 <- c()
  for(k in 1:100){
    set.seed(num*10 + k)
    print(paste0(num,'_',k))
    DE.ind.rep <- sample(DE.ind, num)
    NDE.ind.rep <- sample(NDE.ind, num)
    data2 <- rbind(data3[setdiff(NDE.ind, NDE.ind.rep),], data3[DE.ind.rep,]) 
    dimnames(data2)[[1]] <- 1:dim(data2)[1]
    dimnames(data2)[[2]] <- 1:dim(data2)[2]
    
    #########################################################
    sc <- SCseq(data2)
    sc.temp <- sc
    sc <- filterdata(sc,mintotal=min(colSums(data2))-1, minexpr=min(data2[data2>0]-0.001), minnumber=1)
    fdata <- getfdata(sc)
    features <- sc@cluster$features
    
    
    sc <- compdist(sc,metric="pearson")
    sc <- clustexp(sc)
    sc <- clustexp(sc,cln=1,sat=FALSE)
    sc <- findoutliers(sc)
    
    finalcls <- (sc@cpart)
    print(table(finalcls))
    RaceID.result[[k]] <- finalcls
    rm('sc')
    gc()
    ##################################################################################################
    empty.id <- which(nchar(dimnames(data2)[[1]])==0)
    dimnames(data2)[[1]][empty.id] <- paste0('Define', empty.id)
    pbmc <- CreateSeuratObject(count = data2)
    pbmc <- NormalizeData(object = pbmc, verbose = F)
    
    ## Different Fano genes for clustering
    pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data2)[1], verbose = F)
    vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
    den <- density(vst)
    features.vst <- dimnames(data2)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
    tmp <- log2(data2[dimnames(data2)[[1]] %in% (features.vst),]+1)
    cell.sum <- apply(tmp, 2, sum)
    drop <- apply(tmp, 2, function(x){length(x[x>0])})
    zero.id <- c(which(drop < 3))

    if(length(zero.id)>=1){
      tmp <- tmp[, -zero.id]
    }else{
      tmp <- tmp
    }
    
    pca <- irlba(t(tmp), nv=min(c(50, dim(tmp)-1)))
    pca$pca <-t(pca$d*t(pca$u))
    
    sample.all <- ceiling(dim(tmp)[2] * 0.4)
    knn.res <- Neighbour(pca$pca, pca$pca, k=sample.all)
    
    distance.diff <- knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE]
    
    diff.left <- distance.diff[, -1, drop = FALSE] - distance.diff[, -ncol(distance.diff), drop = FALSE]
    diff.both <- diff.left[, -ncol(diff.left), drop=FALSE] - diff.left[, -1, drop=FALSE]
    diff.both[,1] <- diff.both[,1] + distance.diff[,1]  # Very important due to distance variation to the first neighbor.
    
    v1.k <- matrix(NA, dim(tmp)[2], dim(diff.both)[2])
    skew <- c()
    top.values.ave <- c()
    for(j in 1:dim(diff.both)[2]){
      v <- diff.both[,j]
      v1 <- v
      for(m in 1:length(v)){
        v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
      }
      v1.k[, j] <- (v1)
      v2 <- v1[order(v1, decreasing = T)[(j+2):length(v1)]]
      v2[is.na(v2)] <- 0
      top.values <- v1[knn.res$indices[which.max(v1),1:(j+1)]]
      v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
      skew <- c(skew, skewness(v2))
      top.values.ave <- c(top.values.ave, mean(top.values))
    }
    
    ids <- which(skew > 2)
    if(length(ids)==0){
      predictions <- c(rep(0, dim(tmp)[2]))
    }else{
      col.mat <- matrix(0, length(ids), dim(tmp)[2])
      for(i in 1:length(ids)){
        top.cell <- which.max(v1.k[,(ids[i])])
        col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * top.values.ave[ids[i]]
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
      
      id.match.max <- id.max.match[which.max(skew[id.max.match])]
      rare.cells <- knn.res$indices[which.max(v1.k[,id.match.max]), 1:(id.match.max+1)]
      
      predictions <- ifelse(1:dim(pca$pca)[1] %in% rare.cells, 1, 2)
      
    }
    if(length(zero.id)>=1){
      if(length(table(group[-zero.id]))>1){
        pred <- prediction(c(predictions, rep(2, length(zero.id))), c(group[-zero.id], rep('Group2', length(zero.id))))
        perf <- performance(pred, "auc")
        print(perf@y.values[[1]])
        auc.our <- c(auc.our, perf@y.values[[1]])
      }else{
        print(0.5)
        auc.our <- c(auc.our, 0.5)
      }
      
    }else{
      pred <- prediction(predictions, group)
      perf <- performance(pred, "auc")
      print(perf@y.values[[1]])
      auc.our <- c(auc.our, perf@y.values[[1]])
    }
    
    print(num)
    #########################################################################
    ExprM.normCounts.filter <- data2
    source(paste(Rfundir,"GiniClust2_fitting_DE.R",sep=""))
    ExprM.RawCounts<-ExprM.RawCounts.filter
    
    source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep=""))
    print(length(P_G))
    Gini.result[[k]] <- P_G
    #########################################################################
    
    
    data1 <- t(data2) #Samples * Features
    
    #Genes
    genes <- c(1:dim(data1)[2]) #It can be replaced with original gene names
    
    data_mat <- list(mat=data1, gene_symbols=genes)
    preprocessedList <- ranger_preprocess(data_mat)
    preprocessedData <- as.matrix(preprocessedList$preprocessedData)
    model <- new(FiRE::FiRE, 100, 50, 1017881, 5489, 0)
    model$fit(preprocessedData)
    score <- model$score(preprocessedData)
    
    #Apply IQR-based criteria to identify rare cells for further downstream analysis.
    q3 <- quantile(score, 0.75)
    iqr <- IQR(score)
    th1 <- q3 + (1.5*iqr)
    th2 <- q3 + (1*iqr)
    th3 <- q3 + (0.5*iqr)
    
    #Select indexes that satisfy IQR-based thresholding criteria.
    indIqr <- which(score >= th1)
    
    #Create a file with binary predictions
    predictions.15 <- integer(dim(data1)[1])
    predictions.15[which(score >= th1)] <- 1 #Replace predictions for rare cells with '1'.
    predictions.10 <- integer(dim(data1)[1])
    predictions.10[which(score >= th2)] <- 1 #Replace predictions for rare cells with '1'.
    predictions.05 <- integer(dim(data1)[1])
    predictions.05[which(score >= th3)] <- 1 #Replace predictions for rare cells with '1'.
    
    perf.15 <- performance(prediction(predictions.15, group), "auc")
    perf.10 <- performance(prediction(predictions.10, group), "auc")
    perf.05 <- performance(prediction(predictions.05, group), "auc")
    
    print(perf.15@y.values[[1]])
    print(perf.10@y.values[[1]])
    print(perf.05@y.values[[1]])
    
    auc.fire.15 <- c(auc.fire.15, perf.15@y.values[[1]])
    auc.fire.10 <- c(auc.fire.10, perf.10@y.values[[1]])
    auc.fire.05 <- c(auc.fire.05, perf.05@y.values[[1]])
    
    
    
    
  }
  race.res <- transform.race()
  auc.vec <- c()
  
  for(m in 1:length(RaceID.result)){
    
    predictions <- race.res[[m]]
    if(all(is.na(predictions))){
      auc.vec <- c(auc.vec, NA)
    }else{
      pred <- prediction(predictions, group)
      perf <- performance(pred, "auc")
      
      print(num)
      print(perf@y.values[[1]])
      auc.vec <- c(auc.vec, perf@y.values[[1]])
    }
    
  }
  
  gini.res <- transform.gini()
  auc.gini <- c()
  
  for(n in 1:length(gini.res)){
    
    predictions <- gini.res[[n]]
    if(length(predictions)!=length(group)){
      print(table(predictions))
      perf@y.values[[1]] <- NA
      print(n)
    }else{
      pred <- prediction(predictions, group)
      perf <- performance(pred, "auc")
    }
    
    print(num)
    print(perf@y.values[[1]])
    auc.gini <- c(auc.gini, perf@y.values[[1]])
  }
  
  
  auc.ave.na <- if(all(is.na(auc.vec))) NA else mean(auc.vec, na.rm=T) 
  auc.std.na <- if(all(is.na(auc.vec))) NA else sd(auc.vec, na.rm=T) 
  auc.ave <- c(auc.ave, auc.ave.na)
  auc.std <- c(auc.std, auc.std.na)
  auc.ave.our <- c(auc.ave.our, mean(auc.our))
  auc.std.our <- c(auc.std.our, sd(auc.our))
  auc.ave.gini <- c(auc.ave.gini, mean(auc.gini))
  auc.std.gini <- c(auc.std.gini, sd(auc.gini))
  
  auc.ave.15 <- c(auc.ave.15, mean(auc.fire.15))
  auc.std.15 <- c(auc.std.15, sd(auc.fire.15))
  auc.ave.10 <- c(auc.ave.10, mean(auc.fire.10))
  auc.std.10 <- c(auc.std.10, sd(auc.fire.10))
  auc.ave.05 <- c(auc.ave.05, mean(auc.fire.05))
  auc.std.05 <- c(auc.std.05, sd(auc.fire.05))
  
  rm(list=c('data2', 'pca', 'pp', 'ww', 'knn.res', 'auc.vec', 'auc.our', 'auc.gini', 'auc.fire.15', 'auc.fire.10', 'auc.fire.05',
            'pred', 'perf', 'Gini.result', 'predictions', 'data1', 'data_mat', 'RaceID.result', 'preprocessedList', 'preprocessedData',
            'model', 'ExprM.normCounts.filter', 'ExprM.RawCounts'))
  gc()
  res.num <- list(auc.ave, auc.std, auc.ave.our, auc.std.our, auc.ave.gini, auc.std.gini, auc.ave.15, auc.std.15,
                  auc.ave.10, auc.std.10, auc.ave.05, auc.std.05)
  #save(res.num, file=paste0('./AUC/AUC_', num, '.RData'))
  
}


auc.result <- list(auc.ave, auc.std,
                   auc.ave.our, auc.std.our,
                   auc.ave.gini, auc.std.gini,
                   auc.ave.15, auc.std.15,
                   auc.ave.10, auc.std.10,
                   auc.ave.05, auc.std.05)


save(auc.result, file='./results/AUC_result_ALL_2_cells.RData')





