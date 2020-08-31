library(FiRE)
library(Rcpp)
library(rflann)
library(Seurat)
library(ineq)
library(e1071)
library(irlba)
library(scran)
library(ggplot2)
library(pbapply)
library(future.apply)
workdir              = "/home/sam/Documents/RareCellDetection/Proj/Intestine/"     
# where you put the data and results

setwd(workdir)
source('../../Main/utils.R')
#sourceCpp('/home/sam/Documents/FBT/Single/package/Rare/src/utils.cpp')

matrix_dir = 'GSE123516_RAW/'
setwd(workdir)

library(Matrix)

read.data <- function(file_id){
  barcode.path <- paste0(matrix_dir, file_id, "_barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, file_id, "_genes.tsv.gz")
  matrix.path <- paste0(matrix_dir, file_id, ".mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  return(mat)
}

files <- c('GSM3308717_C04', 'GSM3308718_C05', 'GSM3308719_C06', 'GSM3308720_C07')

mat1 <- read.data(files[1])
mat2 <- read.data(files[2])
#mat3 <- read.data(files[3])
#mat4 <- read.data(files[4])


all(mat1@Dimnames[[1]] == mat2@Dimnames[[1]])
#all(mat1@Dimnames[[1]] == mat3@Dimnames[[1]])
#all(mat1@Dimnames[[1]] == mat4@Dimnames[[1]])
#all(mat2@Dimnames[[1]] == mat4@Dimnames[[1]])
#all(mat3@Dimnames[[1]] == mat4@Dimnames[[1]])
#all(mat2@Dimnames[[1]] == mat3@Dimnames[[1]])


###############################################################
OUR.result <- list()
#for (i in 2:100){
  #k=ks[j]

  
  #data <- ExprM.normCounts.filter[rowMeans(ExprM.normCounts.filter) > 0,]
  #data <- log2(data + 1)
  #disp <- FastLogVMR(as(data, 'dgCMatrix'), F)
  
  data <- as.matrix(mat2)
  drop <- apply(data, 1, function(x){length(x[x>0])})
  data <- data[drop>=2,]
  sample.expressed <- apply(data, 2, function(x){length(x[x>0])})
  data <- data[, sample.expressed>=200]
  norm.data <- .normalize_by_umi(t(data), gene_symbols = dimnames(data)[[1]], minLibSize=0, verbose = F)
  data2 <- t(norm.data$m)
  
  data2 <- data[dimnames(data2)[[1]],]
  sf <- computeSumFactors(data2)
  
  data2 <- t(t(data2)/sf)
  #data2 <- log2(data +1)
  #sc <- SCseq(data2)
  #sc.temp <- sc
  #sc <- filterdata(sc,mintotal=min(colSums(data2))-1, minexpr=min(data2[data2>0]-0.001), minnumber=1)
  #fdata <- getfdata(sc)
  #features <- sc@cluster$features
  
  pbmc <- CreateSeuratObject(count = data2)
  #pbmc <- NormalizeData(object = pbmc, verbose = F) # top vst features were same whether normalization is applied or not for this dataset, as is the same with 68 k PBMC data set.
  
  ## Different Fano genes for clustering
  pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data2)[1], verbose = F)
  vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
  den <- density(vst)
  features.vst <- dimnames(data2)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
  #features.vst <- dimnames(data2)[[1]][vst > mean(vst) + 1.96 * sd(vst)]
  #features.vst <- dimnames(data2)[[1]][order(vst, decreasing = T)[1:2000]]
  ##################################################################################################
  #tmp <- data2[dimnames(data2)[[1]] %in% union(features, features.vst),]
  tmp <- data2[dimnames(data2)[[1]] %in% (features.vst),]
  tmp <- log2(tmp+1)
  #tmp1 <- .normalize_by_umi(t(tmp), gene_symbols = dimnames(tmp)[[1]], minLibSize=0, verbose = F)
  cell.sum <- apply(tmp, 2, sum)
  tmp <- t(t(tmp)/(cell.sum/mean(cell.sum)))
  #tmp <- tmp1$m
  #pca <- calcul.pca(t(log2(tmp+1)), 50)
  pca <- irlba(t(tmp), nv=50) # More robust no error, contrast to calcul.pca
  pca$pca <-t(pca$d*t(pca$u))
  #pca <- calcul.pca(t(data[dimnames(data)[[1]] %in% features,]), 50)
  #pca <- prcomp(t(tmp))
  #pca$pca <- pca$x[,1:min(50, dim(pca$x)[2])]
  #pca <- calcul.pca(t(data[order(disp, decreasing = T)[1:2000],]), 50)
  #knn.res <- Neighbour(pca$pca, pca$pca, k=200)
  
  knn.res <- Neighbour(pca$pca, pca$pca, k=450)
  
  distance.diff <- (knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE])
  
  diff.left <- distance.diff[, -1, drop = FALSE] - distance.diff[, -ncol(distance.diff), drop = FALSE]
  diff.both <- diff.left[, -ncol(diff.left), drop=FALSE] - diff.left[, -1, drop=FALSE]
  diff.both[,1] <- diff.both[,1] + distance.diff[,1]  # Very important due to distance variation to the first neighbor.
  
  v1.k <- matrix(NA, dim(data2)[2], 447)
  skew <- c()
  skew1 <- c()
  top.ave <- c()
  remain.ave <- c()
  for(j in 1:dim(diff.both)[2]){
    #v <- (distance.diff[,j])
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
    top.values <- v1[knn.res$indices[which.max(v1),1:(j+1)]]
    #top.values <- v1[order(v1, decreasing = T)[1:(j)]]
    #v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(exp(mean(log(top.values[top.values>0]))), (2)))
    v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
    skew <- c(skew, skewness(v2))
    
    v1.sort <- sort(v1, decreasing = T)
    top.ave <- c(top.ave, mean(v1.sort[1:(j)]))
    remain.ave <- c(remain.ave, mean(v1.sort[(j+1):length(v1)]))
  }  
  

  
#}
  
ids <- which(skew > 2)
#col <- rep('cls', dim(tmp)[2])  

# for(id in ids){
#   top.cell <- which.max(v1.k[,(id)])
#   col[knn.res$indices[top.cell,1:(id+1)]] <- paste0(col[knn.res$indices[top.cell,1:(id+1)]], '_', id)
#   
# }

col.mat <- matrix(0, length(ids), dim(tmp)[2])
for(i in 1:length(ids)){
  top.cell <- which.max(v1.k[,(ids[i])])
  col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]]
}

#col.mat <- col.mat[, colSums(col.mat)>0]
id.max <- apply(col.mat, 2, which.max)
cnt <- table(id.max)
id.max.match <- as.integer(names(cnt)[which(cnt == ids[sort(unique(id.max))] + 1)])

cls <- rep(0, dim(tmp)[2])
for(id.match in id.max.match){
  cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
}

rare.cells <- list()
for(id.match in id.max.match){
  rare.cells[[as.character(ids[id.match])]] <- knn.res$indices[which.max(v1.k[,ids[id.match]]), 1:(ids[id.match]+1)]
}





tsne <- Rtsne::Rtsne(t(tmp), dims=2, perplexity=30)
#plot(tsne$Y[,1], tsne$Y[,2], cex=0.3, col=factor(col))
tsne.df <- data.frame(TSNE.1=tsne$Y[,1], TSNE.2=tsne$Y[,2], cluster=col)
#save(tsne.df, file='tsne_data.RData')
p <- ggplot(data=tsne.df, aes(x=TSNE.1, y=TSNE.2, colour=cluster)) + geom_point(size=0.3) +
  theme_bw() + theme(legend.position='none',  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.text = element_text(size = 16),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank()) + labs(x='tSNE-1', y='tSNE-2') +
  scale_color_manual(name=NULL, values = c("cls"='grey', "cls_130_149"="red1", "cls_149"="DeepSkyBlue", "cls_28" = "Magenta", 'cls_5'='green3'))
p



genes <- c('Agr2', 'Muc2', 'Tff3', 
           'Cd3e', 'Cd3d', 'Cd3g', 'Ptprc',
           'Lyz1', 'Defa17', 'Defa22', 'Defa24', 'Ang4', 
           'Lgr5', 'Ascl2', 'Axin2', 'Olfm4', 'Slc12a2',
           'Mki67', 'Cdk4', 'Mcm6', 'Mcm5', 'Pcna', 
           'Clu', 'Anxa1', 'Ly6d', 'Areg', 'Lamc2', 
           'Dclk1', 'Trpm5')


#            'Alas2', 'Hbb-bs', 'Hbb-bt', 'Hba-a1', 'Hba-a2',
# 'Cd8a', 'Cd8b1', 'Cd4', 'Cd19'


library(pheatmap)


col.rare <- col[col!= 'cls']
col.rare <- as.character(factor(col.rare, levels=c('cls_130_149', 'cls_149', 'cls_28', 'cls_5'), labels=c('1_1', '1_2', '2', '3')))
data1 <- data[, col!='cls']
annotation_col = data.frame(
  "Clusters" = factor(col.rare, levels=c('1_1', '1_2', '2', '3'))
)
rownames(annotation_col) = colnames(data1)
cols = list(Clusters = c("1_1" = "red1", "1_2"="DeepSkyBlue", '2'='Magenta', '3'='green3'))
plot.df <- log2(data1[(genes), c(which(col.rare=='1_1'), which(col.rare=='1_2'), which(col.rare=='2'), which(col.rare=='3'))]+1)
dimnames(plot.df)[[1]] <- paste0('    ', dimnames(plot.df)[[1]])
pheatmap(plot.df, 
         cluster_cols = F, cluster_rows = F, annotation_col = annotation_col, annotation_names_col=F, show_colnames = F, annotation_colors = cols,
         gaps_col = c(131, 150, 179), legend=F, annotation_legend = F)







data22 <- data2[!duplicated(dimnames(data2)[[1]]),]
data3 <- data22[, col %in% c('cls_149', 'cls_130_149', 'cls_28', 'cls_5')]
grp <-col[col %in% c('cls_149', 'cls_130_149', 'cls_28', 'cls_5')]
# 
vals <- sort(unique(c(data3)))
data3[data3==0] <- vals[2]

## cls_149
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_149'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_149'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_149'])/mean(x[grp!='cls_149']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC1 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE1 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df1 <- data.frame(logFC=logFC1, p_val=DE1, gene=dimnames(data33)[[1]])
df1 <- df1[p.adjust(df1$p_val, 'fdr') < 0.05,]
#gene1 <- as.character(df1$gene[order(df1$logFC, decreasing=T)[1:30]])
#pheatmap(t(log2(data3[unique(c(gene1)),c(which(col.rare=='1_1'), which(col.rare=='1_2'),  which(col.rare=='2'), which(col.rare=='3'))]+1)), 
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
#                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene1 <- c('Cd74', 'H2-Aa', 'H2-Ab1', 'Ly6d', 'Ebf1', 'Cd79a', 'Mef2c', 'H2-Eb1', 'Gm43603', 'Tnfrsf13c')

## cls_130_149
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_130_149'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_130_149'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_130_149'])/mean(x[grp!='cls_130_149']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC2 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE2 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df2 <- data.frame(logFC=logFC2, p_val=DE2, gene=dimnames(data33)[[1]])
df2 <- df2[p.adjust(df2$p_val, 'fdr') < 0.05,]
#gene2 <- as.character(df2$gene[order(df2$logFC, decreasing=T)[1:30]])
#pheatmap(t(log2(data3[unique(c(gene2)),c(which(col.rare=='1_1'), which(col.rare=='1_2'),  which(col.rare=='2'), which(col.rare=='3'))]+1)), 
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
#                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene2 <- c('Gzma', 'Ccl5', 'Gzmb', 'Rgs1', 'Nkg7', 'Cd7', 'Cd3g', 'Fcer1g', 'AW112010', 'Cd8a')

## cls_28
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_28'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_28'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_28'])/mean(x[grp!='cls_28']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC3 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE3 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df3 <- data.frame(logFC=logFC3, p_val=DE3, gene=dimnames(data33)[[1]])
df3 <- df3[p.adjust(df3$p_val, 'fdr') < 0.05,]
#gene3 <- as.character(df3$gene[order(df3$logFC, decreasing=T)[1:30]])
#pheatmap(t(log2(data3[unique(c(gene3)),c(which(col.rare=='1_1'), which(col.rare=='1_2'),  which(col.rare=='2'), which(col.rare=='3'))]+1)), 
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
#                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene3 <- c('Krt18', 'Cd24a', 'Adh1', 'Cystm1', 'Aldh2', 'Dclk1', 'Sh2d6', 'Rgs13', 'Hck', 'Trpm5')


### cls_5
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_5'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_5'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_5'])/mean(x[grp!='cls_5']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC4 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE4 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df4 <- data.frame(logFC=logFC4, p_val=DE4, gene=dimnames(data33)[[1]])
df4 <- df4[p.adjust(df4$p_val, 'fdr') < 0.05,]
gene4 <- as.character(df4$gene[order(df4$logFC, decreasing=T)[1:30]])
#pheatmap(t(log2(data3[unique(c(gene4)),c(which(col.rare=='1_1'), which(col.rare=='1_2'),  which(col.rare=='2'), which(col.rare=='3'))]+1)), 
#                 cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
#                 annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene4 <- c('Hbb-bs', 'Hba-a1', 'Hba-a2', 'Hbb-bt', 'Alas2')

col.rare <- col[col!= 'cls']
col.rare <- as.character(factor(col.rare, levels=c('cls_130_149', 'cls_149', 'cls_28', 'cls_5'), labels=c('1_1', '1_2', '2', '3')))
annotation_row = data.frame(
  "Clusters" = factor(col.rare, levels=c('1_1', '1_2', '2', '3'))
)
rownames(annotation_row) = colnames(data3)
cols = list(Clusters = c("1_1" = "red1", "1_2"="DeepSkyBlue", '2'='Magenta', '3'='green3'))

pheatmap(t(log2(data3[unique(c(gene2, gene1, gene3, gene4)),c(which(col.rare=='1_1'), which(col.rare=='1_2'),  which(col.rare=='2'), which(col.rare=='3'))]+1)), 
        cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
        annotation_names_row=F, legend=T, annotation_legend = T, angle_col = "45")






