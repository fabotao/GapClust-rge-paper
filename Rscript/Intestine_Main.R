library(FiRE)
library(rflann)
library(Seurat)
library(ineq)
library(e1071)
library(irlba)
library(scran)
library(ggplot2)
library(pbapply)
library(future.apply)
homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/Intestine/')  
# where you put the data and results

setwd(workdir)
source('../../Rscript/utils.R')

library(Matrix)

load('Intestine_GSM3308718_C05.RData')
###############################################################
data <- as.matrix(mat2)
norm.data <- .normalize_by_umi(t(data), gene_symbols = dimnames(data)[[1]], minLibSize=199, verbose = F)
data2 <- t(norm.data$m)

#data2 <- data2[dimnames(data2)[[1]],]
#sf <- computeSumFactors(data2)

#data2 <- t(t(data2)/sf)

res <- GapClust::GapClust(data2)


pbmc <- CreateSeuratObject(counts = data2)
pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data2)[1], verbose = F)
vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
den <- density(vst)
features.vst <- dimnames(data2)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
tmp <- data2[features.vst,]
tmp <- log2(as.matrix(tmp)+1)

col <- rep('cls', dim(data2)[2])
for(i in 1:length(res$rare_cell_indices)){
  col[res$rare_cell_indices[[i]]] <- paste0('cls_', i)
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
  scale_color_manual(name=NULL, values = c("cls"='grey', "cls_1"="red1", "cls_2"="DeepSkyBlue", "cls_3" = "Magenta", 'cls_4'='green3',
                                           "cls_5"='blue', "cls_6"="orange"))
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
col.rare <- as.character(factor(col.rare, levels=c('cls_1', 'cls_2', 'cls_3', 'cls_4', 'cls_5', 'cls_6'), labels=c('1', '2', '3', '4', '5', '6')))
data1 <- data[, col!='cls']
annotation_col = data.frame(
  "Clusters" = factor(col.rare, levels=c('1', '2', '3', '4', '5', '6'))
)
rownames(annotation_col) = colnames(data1)
cols = list(Clusters = c("1" = "red1", "2"="DeepSkyBlue", '3'='Magenta', '4'='green3', '5'='blue', '6'='orange'))
plot.df <- log2(data2[(genes), c(which(col=='cls_1'), which(col=='cls_2'), which(col=='cls_3'), which(col=='cls_4'), 
                                 which(col=='cls_5'), which(col=='cls_6'))]+1)
dimnames(plot.df)[[1]] <- paste0('    ', dimnames(plot.df)[[1]])
pheatmap(plot.df, 
         cluster_cols = F, cluster_rows = F, annotation_col = annotation_col, annotation_names_col=F, show_colnames = F, annotation_colors = cols,
         gaps_col = c(6, 13, 27, 54, 83), legend=F, annotation_legend = F)







data22 <- data2[!duplicated(dimnames(data2)[[1]]),]
data3 <- data22[, col %in% c('cls_1', 'cls_2', 'cls_3', 'cls_4', 'cls_5', 'cls_6')]
grp <- col[col %in% c('cls_1', 'cls_2', 'cls_3', 'cls_4', 'cls_5', 'cls_6')]
# 
vals <- sort(unique(c(data3)))
data3[data3==0] <- vals[2]

## cls_1
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_1'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_1'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_1'])/mean(x[grp!='cls_1']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC1 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE1 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df1 <- data.frame(logFC=logFC1, p_val=DE1, gene=dimnames(data33)[[1]])
df1 <- df1[p.adjust(df1$p_val, 'fdr') < 0.05,]
gene1 <- as.character(df1$gene[order(df1$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene1)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene1 <- c('Hba-a1', 'Hbb-bs', 'Hba-a2', 'Hbb-bt', 'Alas2')

#gene1 <- c('Cd74', 'H2-Aa', 'H2-Ab1', 'Ly6d', 'Ebf1', 'Cd79a', 'Mef2c', 'H2-Eb1', 'Gm43603', 'Tnfrsf13c')

## cls_2
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_2'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_2'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_2'])/mean(x[grp!='cls_2']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC2 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE2 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df2 <- data.frame(logFC=logFC2, p_val=DE2, gene=dimnames(data33)[[1]])
df2 <- df2[p.adjust(df2$p_val, 'fdr') < 0.05,]
gene2 <- as.character(df2$gene[order(df2$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene2)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene2 <- c('Spink4', 'Retnlb', 'Mptx1', 'Mmp7', "Ccl6")


#gene2 <- c('Gzma', 'Ccl5', 'Gzmb', 'Rgs1', 'Nkg7', 'Cd7', 'Cd3g', 'Fcer1g', 'AW112010', 'Cd8a')

## cls_3
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_3'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_3'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_3'])/mean(x[grp!='cls_3']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC3 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE3 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df3 <- data.frame(logFC=logFC3, p_val=DE3, gene=dimnames(data33)[[1]])
df3 <- df3[p.adjust(df3$p_val, 'fdr') < 0.05,]
gene3 <- as.character(df3$gene[order(df3$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene3)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                    cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                    annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene3 <- c('Afp', 'Chga', 'Tac1', 'Reg4', 'Tph1', 'Chgb')

# Ref: Single-cell messenger RNA sequencing reveals rare intestinal cell types


#gene3 <- c('Krt18', 'Cd24a', 'Adh1', 'Cystm1', 'Aldh2', 'Dclk1', 'Sh2d6', 'Rgs13', 'Hck', 'Trpm5')


### cls_4
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_4'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_4'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_4'])/mean(x[grp!='cls_4']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC4 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE4 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df4 <- data.frame(logFC=logFC4, p_val=DE4, gene=dimnames(data33)[[1]])
df4 <- df4[p.adjust(df4$p_val, 'fdr') < 0.05,]
gene4 <- as.character(df4$gene[order(df4$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene4)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                 cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                 annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene4 <- c('Dclk1', 'Rgs13', 'Hck', 'Gng13', 'Trpm5')


#gene4 <- c('Hbb-bs', 'Hba-a1', 'Hba-a2', 'Hbb-bt', 'Alas2')


### cls_5
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_5'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_5'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_5'])/mean(x[grp!='cls_5']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC5 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE5 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df5 <- data.frame(logFC=logFC5, p_val=DE5, gene=dimnames(data33)[[1]])
df5 <- df5[p.adjust(df4$p_val, 'fdr') < 0.05,]
gene5 <- as.character(df5$gene[order(df5$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene5)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                 cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                 annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene5 <- c('Gip', 'Cck', 'Fabp5', 'Rbp4', 'Scg2')


#gene4 <- c('Hbb-bs', 'Hba-a1', 'Hba-a2', 'Hbb-bt', 'Alas2')

### cls_6
DE <- WilcoxDETest((data3), cells.1 = dimnames(data3)[[2]][grp=='cls_6'], 
                   cells.2 = dimnames(data3)[[2]][grp!='cls_6'], verbose = F)
logFC <- apply(data3, 1, function(x){log(mean(x[grp=='cls_6'])/mean(x[grp!='cls_6']))})

data33 <- data3[!is.na(DE$p_val) & !is.na(logFC),]
logFC6 <- logFC[!is.na(DE$p_val) & !is.na(logFC)]
DE6 <- DE$p_val[!is.na(DE$p_val) & !is.na(logFC)]
df6 <- data.frame(logFC=logFC6, p_val=DE6, gene=dimnames(data33)[[1]])
df6 <- df6[p.adjust(df6$p_val, 'fdr') < 0.05,]
gene6 <- as.character(df6$gene[order(df6$logFC, decreasing=T)[1:30]])
pheatmap(t(log2(data3[unique(c(gene6)),c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))]+1)), 
                 cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_colors = cols,
                 annotation_names_row=F, legend=T, annotation_legend = F, angle_col = "45")
gene6 <- c('Gzma', 'Ccl5', 'Rgs1', 'Cd3g', 'Cd7', 'Gzmb')



#gene4 <- c('Hbb-bs', 'Hba-a1', 'Hba-a2', 'Hbb-bt', 'Alas2')







col.rare <- col[col!= 'cls']
col.rare <- as.character(factor(col.rare, levels=c('cls_1', 'cls_2', 'cls_3', 'cls_4', 'cls_5', 'cls_6'), labels=c('1', '2', '3', '4', '5', '6')))
annotation_row = data.frame(
  "Clusters" = factor(col.rare, levels=c('1', '2', '3', '4', '5', '6'))
)
rownames(annotation_row) = colnames(data3)
cols = list(Clusters = c("1" = "red1", "2"="DeepSkyBlue", '3'='Magenta', '4'='green3', '5'='blue', '6'='orange'))
ind <-  c(which(grp=='cls_1'), which(grp=='cls_2'),  which(grp=='cls_3'), which(grp=='cls_4'), which(grp=='cls_5'), which(grp=='cls_6'))
pheatmap(t(log2(data3[unique(c(gene1, gene2, gene3, gene4, gene5, gene6)), ind]+1)), 
        cluster_cols = F, cluster_rows = F, show_rownames = F, annotation_row = annotation_row, annotation_colors = cols,
        annotation_names_row=F, legend=T, annotation_legend = T, angle_col = "45",
        gaps_row = c(6, 13, 27, 54, 83), gaps_col=c(length(gene1), length(gene1) + length(gene2), length(gene1) + length(gene2) + length(gene3),
                                                    length(gene1) + length(gene2) + length(gene3) + length(gene4), 
                                                    length(gene1) + length(gene2) + length(gene3) + length(gene4) + length(gene5)))






