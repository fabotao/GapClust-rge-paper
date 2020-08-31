workdir              = "/home/sam/Documents/RareCellDetection/Proj/jurkat_two_species/"    
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
library(Seurat)
library(scran)
library(Matrix)
source('../../Main/utils.R')

matrix_dir <- 'data/'
read.data <- function(){
  barcode.path <- paste0(matrix_dir, "barcodes.tsv")
  features.path <- paste0(matrix_dir, "genes.tsv")
  matrix.path <- paste0(matrix_dir, "matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  return(mat)
}

mat <- read.data()
cluster <- read.csv('data/clusters.csv', header=T, stringsAsFactors = F)

data <- mat[rowMeans(mat) > 0,]
tmp <- data
tmp@x <- rep(1, length(tmp@x))

data1 <- data[rowSums(tmp)>=2,]

sf <- computeSumFactors(data1)
data1 <- t(t(data1)/sf)

all(cluster$Barcode == data1@Dimnames[[2]])

# pbmc <- CreateSeuratObject(count = data1)
# pbmc <- NormalizeData(object = pbmc, verbose=F)
# pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data1)[1], verbose = F)
# vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
# den <- density(vst)
# features.vst <- dimnames(data1)[[1]][vst > (den$x[which.max(den$y)] - min(den$x)) * 2 + min(den$x)]
# 
# data2 <- as.matrix(data1[features.vst,])
# 
# Heatmap(cbind(data2[, cluster$Cluster==1], data2[, cluster$Cluster==2]), cluster_columns = F)

ref.data <- read.table('data/preprocessedData_jurkat_two_species_1580.txt', header=F, sep=' ')
label <- read.csv('data/labels_jurkat_two_species_1580.txt', header=F)
genes <- read.csv('data/genesSel_jurkat_two_species_1580.txt', header=F)
ref.data.raw <- read.table('data/jurkat_two_species_1580.txt', header=F, sep=' ')

data.all <- as.matrix(mat[genes$V1,])
data.ref <- t(as.matrix(ref.data))
data.ref.raw <- as.matrix(ref.data.raw[genes$V1,])

type.ref1 <- rowMeans(data.all[, cluster$Cluster==1])
type.ref2 <- rowMeans(data.all[, cluster$Cluster==2])

corr1 <- cor((data.ref), type.ref1)
corr2 <- cor((data.ref), type.ref2)
plot(corr1, corr2, cex=0.2)  # cluster == 2 means Jurkat cells (Rare cell type)

cell.293t <- data.all[, cluster$Cluster==1]
cell.jurkat <- data.all[, cluster$Cluster==2]


# Find the 1540 293 cell from raw data
cell.293t.all <- mat[, cluster$Cluster==1]
cell.293.ref <- data.ref.raw[, label$V1==1]
ids <- c()
for(i in 1:dim(cell.293.ref)[2]){
  
  ids <- c(ids, which.max(cor((cell.293t), cell.293.ref[,i])))
  
}

#cell.293.1580 <- cell.293t.all[, ids]

gene <- dimnames(mat)[[1]][genes$V1]
tsne <- Rtsne::Rtsne(t(as.matrix(data1[dimnames(data1)[[1]] %in% gene,])), dims=2, perplexity=30)

corr1 <- (cor(log2(as.matrix(data1[dimnames(data1)[[1]] %in% gene,]) + 1)))
corr2 <- cor(corr1, cluster$Cluster)

den1 <- density(-corr2[corr2>0])
den2 <- density(corr2[corr2<0])
elbow1 <- find_elbow(den1$x[which.max(den1$y):length(den1$x)], den1$y[which.max(den1$y):length(den1$y)])
elbow2 <- find_elbow(den2$x[which.max(den2$y):length(den2$x)], den2$y[which.max(den2$y):length(den2$y)])

cell.293t.id <- (order(corr2)[1:1540])
cell.jurkat.id <- order(corr2, decreasing = T)[1:300]

# # Find the 1540 293 cell from raw data
# cell.jurkat.all <- mat[, cluster$Cluster==2]
# cell.jurkat.ref <- data.ref.raw[, label$V1==2]
# ids1 <- c()
# for(j in 1:dim(cell.jurkat.ref)[2]){
#   
#   ids1 <- c(ids1, which.max(cor((cell.jurkat), cell.jurkat.ref[,j])))
#   
# }
# 
set.seed(2019)
#ids2 <- sample(setdiff(1:dim(cell.jurkat.all)[2], ids1), 82-40)
#cell.jurkat.1580 <- cell.jurkat.all[, union(ids1, ids2)]

cell.293.1580 <- as.matrix(data1[, cell.293t.id])
cell.jurkat.filtered <- as.matrix(data1[, cell.jurkat.id])


#> 1540 * 0.05/(1-0.05)
#[1] 81.05263


dilut.ratios <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)

for(ratio in dilut.ratios){
  jurkat.num <- ceiling(1540 * ratio/(1-ratio))
  jurkat.293 <- cbind(cell.293.1580, cell.jurkat.filtered[, sample(1:dim(cell.jurkat.filtered)[2], jurkat.num)])
  jurkat.293 <- jurkat.293[rowSums(jurkat.293) > 0,]
  dimnames(jurkat.293)[[2]] <- paste0(dimnames(jurkat.293)[[2]], ifelse(1:dim(jurkat.293)[2] %in% 1:1540, '-293T', '-Jurkat'))
  save(jurkat.293, file=paste0('./results/Jurkat_', jurkat.num, '_293_1580_cells.RData'))
}

