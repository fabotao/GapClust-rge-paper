.normalize_by_umi <-function(mat, gene_symbols, minLibSize, verbose){
   
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
    #write.csv(keepCellIndex, file = paste(dataSave,"keepCellIndexes.csv",sep=""), quote = F,row.names = F)
  }
  
  #Filter Genes
  cs <- colSums(mat>2)
  cs <- apply(mat, 2, function(x){length(x[x>0])})
  x_use_genes <- which(cs > 2)
  
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


# Differential expression using Wilcoxon Rank Sum
#
# Identifies differentially expressed genes between two groups of cells using
# a Wilcoxon Rank Sum test
#
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param verbose Print a progress bar
# @param ... Extra parameters passed to wilcox.test
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# features
#
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#
# @export
#
# @examples
# pbmc_small
# WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(x = group.info), drop = FALSE]
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
    }
  )
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}




calcul.pca = function(x,n){
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  ppk<-svd:::propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(pca=pca, use_genes=use_genes)
}


which.peaks <- function(x,partial=TRUE,decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else {
      which(diff(diff(x)>0)>0)+1
    }
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)+1
    }
  }
}

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


find.inflection <- function(v1, v2, k=10){
  d2 <- diff(v2)
  d2 <- d2>0
  d2 <- d2*2 -1 
  k=k
  cutoff <- 10
  scores <- sapply(k:(length(d2)-k), FUN=function(i){
    score <- abs(mean(-d2[ i-1:k ], na.rm=T) + mean(d2[ i+0:k ], na.rm=T))
  })
  
  scores <- sapply(k:(length(v2)-k), FUN=function(i){
    left <- (v2[sapply(i-1:k, max, 1) ]<v2[i])*2-1
    right <- (v2[sapply(i+1:k, min, length(v2)) ]<v2[i])*2-1
    
    score <- abs(sum(left) + sum(right))
  })
  
  inflections <- (k:(length(v2)-k))[scores>=cutoff]
  return(inflections)
}

#########################################
fitbackground <- function(x,mthr=-1){
  m <- apply(x,1,mean)
  v <- apply(x,1,var )
  
  ml <- log2(m)
  vl <- log2(v)
  f <- ml > -Inf & vl > -Inf
  ml <- ml[f]
  vl <- vl[f]
  mm <- -8
  repeat{
    fit <- lm(vl ~ ml + I(ml^2)) 
    if( coef(fit)[3] >= 0 | mm >= mthr){
      break
    }
    mm <- mm + .5
    f <- ml > mm
    ml <- ml[f]
    vl <- vl[f]
  }
  
  vln <- log2(v)  - log2(sapply(m,FUN=uvar,fit=fit))
  n <- names(vln)[vln>0]
  return(list(fit=fit,n=n))
}

lvar  <- function(x,fit) 2**(coef(fit)[1] + log2(x)*coef(fit)[2] + coef(fit)[3] * log2(x)**2)
lsize <- function(x,lvar,fit) x**2/(max(x + 1e-6,lvar(x,fit)) - x)
uvar  <- function(x,fit){
  err <- coef(summary(fit))[, "Std. Error"]
  2**(coef(fit)[1] + err[1] + log2(x)*(coef(fit)[2] + err[2]) + (coef(fit)[3] + err[3]) * log2(x)**2)
}


############################################################################
#fits a line to the relationship of log(variance) and log(mean) using
#local polynomial regression (loess). Then standardizes the feature values
#using the observed mean and expected variance (given by the fitted line).
#Feature variance is then calculated on the standardized values after clipping
#to a maximum (see clip.max parameter).
###########################################################################
fitvst <- function(object){
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as(object = object, Class = 'dgCMatrix')
  }
  hvf.info <- data.frame(mean = rowMeans(x = object))
  hvf.info$variance <- SparseRowVar2(
    mat = object,
    mu = hvf.info$mean,
    display_progress = F
  )
  hvf.info$variance.expected <- 0
  hvf.info$variance.standardized <- 0
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, ],
    span = loess.span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  # use c function to get variance after feature standardization
  hvf.info$variance.standardized <- SparseRowVarStd(
    mat = object,
    mu = hvf.info$mean,
    sd = sqrt(hvf.info$variance.expected),
    vmax = clip.max,
    display_progress = verbose
  )
  colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  
  
}
