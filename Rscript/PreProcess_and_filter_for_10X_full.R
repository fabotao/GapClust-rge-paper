#Preprocess and filter 10X data using 10X code: https://github.com/10XGenomics/single-cell-3prime-paper

wd<-"/home/sam/Documents/RareCellDetection/Proj/10X_full/single-cell-3prime-paper-master/pbmc68k_analysis/"
setwd(wd)
source(paste(wd,file.path('util.R'),sep=""))
pbmc_68k <- readRDS(paste("../../data/",file.path('pbmc68k_data.rds'),sep=""))
all_data <- pbmc_68k$all_data
m<-all_data[[1]]$hg19$mat
genes <- all_data[[1]]$hg19$genes
genesymbols<-all_data[[1]]$hg19$gene_symbols
barcodes <- all_data[[1]]$hg19$barcodes

#ExprM.RawCounts<-as.matrix(t(m))
#rownames(ExprM.RawCounts)<-genes



l<-.normalize_by_umi(m) 
genes_used<-genes[l$use_genes]


data1 <- (t(l$m))
rownames(data1) <- genes_used
colnames(data1) <- barcodes

data4 <- m[, which(genes %in% genes_used)]
data4 <- t(data4)
rownames(data4) <- genes_used
colnames(data4) <- barcodes

library(scran)
sf <- computeSumFactors(data4)

data4 <- t(t(data4)/sf)

gene.symbols <- genesymbols[which(genes %in% genes_used)]

save(data4, file='../../data/PBMC_68K_normalized.RData')
save(gene.symbols, file='../../data/PBMC_68K_normalized_gene_symbols.RData')

#save(ExprM.RawCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
#save(ExprM.RawCounts.filter, file=paste("results/", exprimentID, "_ExprM.filter.RData", sep=""))
