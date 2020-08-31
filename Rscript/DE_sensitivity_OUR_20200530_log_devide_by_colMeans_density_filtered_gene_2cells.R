workdir              = "/home/sam/Documents/RareCellDetection/Proj/10X_subsample_A/"     
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
library(Seurat)
library(irlba)
library(scran)
source('../../Main/utils.R')
sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')

auc.mat <- matrix(NA, 100, 200)

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

## save DE and non-DE genes for plots
#DE_dat <- data3[DE.ind,]
#NDE_dat <- data3[NDE.ind,]
#save(DE_dat, file = 'DE_genes_20200530.RData')
#save(NDE_dat, file = 'NDE_genes_20200530.RData')



set.seed(2019)
DE.range <- 1:length(DE.ind)
auc.ave <- c()
auc.std <- c()

for(num in DE.range){
  auc.vec <- c()
  for(k in 1:100){
    set.seed(num*10 + k)
    print(paste0(num,'_',k))
    DE.ind.rep <- sample(DE.ind, num)
    NDE.ind.rep <- sample(NDE.ind, num)
    data <- rbind(data3[setdiff(NDE.ind, NDE.ind.rep),], data3[DE.ind.rep,]) 
    empty.id <- which(nchar(dimnames(data)[[1]])==0)
    dimnames(data)[[1]][empty.id] <- paste0('Define', empty.id)
    pbmc <- CreateSeuratObject(count = data)
    pbmc <- NormalizeData(object = pbmc, verbose = F)

    ## Different Fano genes for clustering
    pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
    vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
    den <- density(vst)
    features.vst <- dimnames(data)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
    #features.vst <- dimnames(data)[[1]][vst > mean(vst) + 1.96 * sd(vst)]

    tmp <- log2(data[dimnames(data)[[1]] %in% (features.vst),]+1)
    cell.sum <- apply(tmp, 2, sum)
    drop <- apply(tmp, 2, function(x){length(x[x>0])})
    zero.id <- c(which(drop < 3))
    #if(length(zero.id)>=1){
    #  tmp <- t(t(tmp[,-zero.id])/(cell.sum[-zero.id]/(sum(cell.sum)/(length(cell.sum)-length(zero.id)))))
    #}else{
    #  tmp <- t(t(tmp)/(cell.sum/(sum(cell.sum)/(length(cell.sum)))))
    #}    
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
    
    v1.k <- matrix(NA, dim(tmp)[2], sample.all-3)
    skew <- c()
    skew1 <- c()
    top.ave <- c()
    remain.ave <- c()
    top.val <- c()
    for(j in 1:(sample.all-3)){
      
      v <- diff.both[,j]
      v1 <- v
      for(m in 1:length(v)){
        #v1[m] <- (v[m] + sum(v[knn.res$indices[m,2:(1+ceiling(log2(j+1)))]]))/(1+length(2:(1+ceiling(log2(j+1)))))
        v1[m] <- (v[m] + (v[knn.res$indices[m,2]]))/2
        #v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
      }
      v1.k[, j] <- (v1)
      
      v2 <- v1[order(v1, decreasing = T)[(j+1):length(v1)]]
      skew1 <- c(skew1, skewness(v2))
      top.values <- v1[order(v1, decreasing = T)[1:(j)]]
      v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(exp(sum(log(top.values[top.values>0]))/length(top.values)), 2))
      #v2[v2<0] <- 0
      skew <- c(skew, skewness(v2))
      top.val <- c(top.val, exp(sum(log(top.values[top.values>0]))/length(top.values)))
      v1.sort <- sort(v1, decreasing = T)
      top.ave <- c(top.ave, mean(v1.sort[1:(j)]))
      remain.ave <- c(remain.ave, mean(v1.sort[j:length(v1)]))
    }
    #plot(1:197, skew, cex=0.3, col=ifelse(1:197 %in% (i-1), 'red', 'black'))
    skew.max <- which.max(skew)
    if(skew.max==1){
      top.cell <- which.max(v1.k[,(skew.max)])
    }else{
      top.cell <- which.max(v1.k[,(skew.max)])
    }
    
    rare.cells <- knn.res$indices[top.cell,1:(skew.max+1)]
    plot(1:length(skew), skew, col=c('red', rep('black', length(skew)-1)), cex=0.3)
    predictions <- ifelse(1:dim(pca$pca)[1] %in% rare.cells, 1, 2)
    if(length(zero.id)>=1){
      if(length(table(group[-zero.id]))>1){
        pred <- prediction(c(predictions, rep(2, length(zero.id))), c(group[-zero.id], rep('Group2', length(zero.id))))
        perf <- performance(pred, "auc")
        print(perf@y.values[[1]])
        auc.vec <- c(auc.vec, perf@y.values[[1]])
      }else{
        print(0.5)
        auc.vec <- c(auc.vec, 0.5)
      }

    }else{
      pred <- prediction(predictions, group)
      perf <- performance(pred, "auc")
      print(perf@y.values[[1]])
      auc.vec <- c(auc.vec, perf@y.values[[1]])
    }

    rm(list=c('zero.id'))
    gc()
  }
  auc.vec[auc.vec < 0.5] <- 0.5
  auc.ave <- c(auc.ave, mean(auc.vec))
  auc.std <- c(auc.std, sd(auc.vec))
  write.csv(auc.ave, file='auc_2.csv')
}


auc.our <- list(auc.ave, auc.std)


save(auc.our, file='./results/AUC_res_OUR_20200530_log_devide_by_colMeans_density_filtered_gene_2cells.RData')





