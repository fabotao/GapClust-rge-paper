library(ggplot2)
library(ComplexHeatmap)
library(splatter)
workdir              = "/home/sam/Documents/RareCellDetection/Proj/10X_subsample_A/"     
# where you put the data and results

setwd(workdir)


load('results/RaceIDClustering.RData') # RaceID.result
load('results/OURClustering.RData') # GapClust.result
load('results/GiniClustering.RData') # Gini.result
load('results/FiREClustering.RData') # FiRE.result
load('results/FiREClustering10.RData') # FiRE.result10
load('results/FiREClustering05.RData') # FiRE.result05
load('results/CELLSIUSClustering.RData') # CELLSIUS.result

true.res <- RaceID.result
for (i in 2:100){
  #k=ks[j]
  exprimentID<-paste("10X_rare_",i,"_CD14_cells",sep="")
  load(paste0('results/', exprimentID, '_ExprM.filter.RData'))
  true.res[[i-1]] <- ifelse(grepl('CD14', dimnames(ExprM.normCounts.filter)[[2]]), 1, 0)
  rm(list=c('ExprM.normCounts.filter', 'exprimentID'))
}

Gini.result <- sapply(Gini.result, as.character)

transform.gini <- function(){
  gini.res <- list()
  for(i in 1:length(Gini.result)){
    cnt <- table(Gini.result[[i]])
    cnt.top <- table(Gini.result[[i]][1:(i+1)])
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
    cnt.top <- table(RaceID.result[[i]][1:(i+1)])
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


transform.cell <- function(){
  cell.res <- list()
  for(i in 1:length(CELLSIUS.result)){
    cell.res[[i]] <- (rep(0, length(CELLSIUS.result[[i]])))
    cnt.top <- table(CELLSIUS.result[[i]][1:(i+1)])
    cnt.top <- cnt.top[cnt.top > 0 & grepl('_', names(cnt.top))]
    cnt <- table(CELLSIUS.result[[i]])
    if(length(cnt.top)>=1){
      for(j in 1:length(cnt.top)){
        label <- names(cnt.top)[j]
        if(cnt[label] < i+1){
          cell.res[[i]][CELLSIUS.result[[i]] == label] <- 1
        }
      }
    }else{
      cell.res[[i]] <- (rep(0, length(CELLSIUS.result[[i]])))
    }
  }
  return(cell.res)
}


gini.res <- transform.gini()
race.res <- transform.race()
cell.res <- transform.cell()

for(i in 1:99){
  print(table(RaceID.result[[i]]))
   print(table(race.res[[i]]))
}



calcul.f1 <- function(res){
  f1.vec <- c()
  for(i in 1:length(res)){
    TP <- length(which(res[[i]] == 1 & true.res[[i]] == 1))
    FP <- length(which(res[[i]] == 1 & true.res[[i]] == 0))
    FN <- length(which(res[[i]] == 0 & true.res[[i]] == 1))
    if(TP ==0 & FP ==0){
      P = 0
    }else{
      P= TP / (TP + FP)
    }
    if(TP == 0 & FN == 0){
      R = 0
    }else{
      R = TP / (TP + FN)
    }
    if(P == 0 & R == 0){
      f1.vec <- c(f1.vec, 0)
    }else{
      f1.vec <- c(f1.vec, 2*P*R/(P+R))
    }
  }
  return(f1.vec)
}


f1 <- calcul.f1(gini.res)
f2 <- calcul.f1(OUR.result)
f3 <- calcul.f1(race.res)
f4 <- calcul.f1(FiRE.result)
f5 <- calcul.f1(FiRE.result10)
f6 <- calcul.f1(FiRE.result05)
f7 <- calcul.f1(cell.res)

df <- data.frame(cells=rep(2:100, 4), f1_score=c(f1, f3, f7, f2), method=c(rep('GiniClust', 99), rep('RaceID', 99), rep('CellSIUS', 99), rep('GapClust', 99)))
df$method <- factor(df$method, levels=c('GiniClust', 'RaceID', 'CellSIUS', 'GapClust'), ordered = T)
sp <- ggplot(df, aes(x=cells, y=f1_score, fill=factor(method), color=factor(method))) + geom_point(size=2.5)
sp + facet_grid(. ~ method) + theme(panel.grid.minor.x=element_line(size = 0.3, linetype = 'solid',colour = "grey90"), 
                                    panel.grid.minor.y=element_line(size = 0.3, linetype = 'solid',colour = "grey90"), 
                                    panel.background = element_rect(fill = "white", linetype = "solid"),
                                    panel.grid.major = element_line(size = 0.8, linetype = 'solid',colour = "grey90"),
                                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                                    legend.position="none",
                                    strip.text.x = element_text(size = 20, colour = "black", angle = 0),
                                    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                                    axis.text.x = element_text(size=14, colour='black'),
                                    axis.text.y = element_text(size=14, colour='black'),
                                    axis.title.x = element_text(size=18),
                                    axis.title.y = element_text(size=16),
                                    axis.ticks.x = element_line(size=1),
                                    axis.ticks.y = element_line(size=1)) + labs(x='# of rare cells', y=expression(paste(F["1 "], score))) +
  scale_color_manual(values = c("GiniClust" = "DeepSkyBlue","RaceID" = "red1","CellSIUS" = "Magenta", 'GapClust' = 'green3')) +
  scale_x_continuous(breaks=c(2, 20, 40, 60, 80, 100), labels=c('2', '20', '40', '60', '80', '100'))


df <- data.frame(cells=rep(2:100, 4), f1_score=c(f4, f5, f6, f2), method=c(rep('FiRE(1.5×IQR)', 99), rep('FiRE(1.0×IQR)', 99), rep('FiRE(0.5×IQR)', 99), rep('GapClust', 99)))
df$method <- factor(df$method, levels=c('FiRE(1.5×IQR)', 'FiRE(1.0×IQR)', 'FiRE(0.5×IQR)', 'GapClust'), ordered = T)
sp <- ggplot(df, aes(x=cells, y=f1_score, fill=factor(method), color=factor(method))) + geom_point(size=2.5)
sp + facet_grid(. ~ method) + theme(panel.grid.minor.x=element_line(size = 0.3, linetype = 'solid',colour = "grey90"), 
                                    panel.grid.minor.y=element_line(size = 0.3, linetype = 'solid',colour = "grey90"), 
                                    panel.background = element_rect(fill = "white", linetype = "solid"),
                                    panel.grid.major = element_line(size = 0.8, linetype = 'solid',colour = "grey90"),
                                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                                    legend.position="none",
                                    strip.text.x = element_text(size = 20, colour = "black", angle = 0),
                                    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                                    axis.text.x = element_text(size=14),
                                    axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=16),
                                    axis.title.y = element_text(size=16),
                                    axis.ticks.x = element_line(size=1),
                                    axis.ticks.y = element_line(size=1)) + labs(x='# of rare cells', y=expression(paste(F["1 "], score))) +
  scale_color_manual(values = c("FiRE(1.5×IQR)" = "DeepSkyBlue","FiRE(1.0×IQR)" = "red1","FiRE(0.5×IQR)" = "Magenta", 'GapClust' = 'green3')) +
  scale_x_continuous(breaks=c(2, 20, 40, 60, 80, 100), labels=c('2', '20', '40', '60', '80', '100'))






load('results/OURClustering_order.RData')  # FiRE.result.order
load('results/FiREClustering_order.RData')  # GapClust.result.order

mat1 <- matrix(NA, 100, 99)
mat2 <- matrix(NA, 100, 99)


for(m in 1:length(OUR.result.order)){
  score <- OUR.result.order[[m]]
  top100 <- score[order(score, decreasing = T)[1:100]]
  mat1[,m] <- ifelse(grepl('CD14', names(top100)), 1, 0)
  
}


for(n in 1:length(FiRE.result.order)){
  score <- FiRE.result.order[[n]]
  top100 <- score[order(score, decreasing = T)[1:100]]
  mat2[,n] <- ifelse(grepl('CD14', names(top100)), 1, 0)
  
}

dimnames(mat2)[[1]] <- rep(' ',  dim(mat1)[1])
dimnames(mat2)[[1]][c(2, 20, 40, 60, 80, 99)] <- c('1', '20', '40', '60', '80', '100')

dimnames(mat2)[[2]] <- rep(' ',  dim(mat1)[1])
dimnames(mat2)[[2]][c(2, 20, 40, 60, 80, 99)] <- c('1', '20', '40', '60', '80', '100')

dimnames(mat1)[[2]] <- rep(' ',  dim(mat1)[2])
dimnames(mat1)[[2]][c(2, 19, 39, 59, 79, 99)] <- c('2', '20', '40', '60', '80', '100')

dimnames(mat2)[[2]] <- rep(' ',  dim(mat2)[2])
dimnames(mat2)[[2]][c(2, 19, 39, 59, 79, 99)] <- c('2', '20', '40', '60', '80', '100')

colors = structure(c('grey90', 'grey90'), names = c("0", "1"))
ht1 <- Heatmap(mat1, col = colors, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, column_title = "GapClust",
               row_names_side = c("left"), border = T,
               column_names_gp = gpar(fontsize=18), row_names_gp = gpar(fontsize=18), column_title_gp = gpar(fontsize = 24),
        rect_gp = gpar(col = "white", lwd = 0.5), column_labels = colnames(mat1), column_names_rot=0, column_gap = unit(6, 'mm'),
        layer_fun = function(j, i, x, y, w, h, fill) {
          ind_mat = restore_matrix(j, i, x, y)
          for(ir in seq_len(nrow(ind_mat))) {
            # start from the second column
            for(ic in seq_len(ncol(ind_mat))) {
              ind = ind_mat[ir, ic]   # current column
              v = mat1[i[ind], j[ind]]
              if(v > 0) { # if they have the same sign
                col = "steelblue3"
                grid.points(x[c(ind)], y[c(ind)], 
                            pch = 19, gp = gpar(col = col), size = unit(1.4, "mm"))
              }
            }
          }
        })
ht2 <- Heatmap(mat2, col = colors, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, column_title = "FiRE",
               row_names_side = c("left"), row_title_gp = gpar(fontsize = 20), row_title = "Rank", border = T,
               column_names_gp = gpar(fontsize=18), row_names_gp = gpar(fontsize=18), column_title_gp = gpar(fontsize = 24),
        rect_gp = gpar(col = "white", lwd = 0.5), column_labels = colnames(mat1), column_names_rot=0, column_gap = unit(6, 'mm'),
        layer_fun = function(j, i, x, y, w, h, fill) {
          ind_mat = restore_matrix(j, i, x, y)
          for(ir in seq_len(nrow(ind_mat))) {
            # start from the second column
            for(ic in seq_len(ncol(ind_mat))) {
              ind = ind_mat[ir, ic]   # current column
              v = mat2[i[ind], j[ind]]
              if(v > 0) { # if they have the same sign
                col = "steelblue3"
                grid.points(x[c(ind)], y[c(ind)], 
                            pch = 19, gp = gpar(col = col), size = unit(1.4, "mm"))
              }
            }
          }
        })


ht_list = ht2 + ht1
ComplexHeatmap::draw(ht_list, ht_gap = unit(0.5, "cm"), column_title='# of rare cells', column_title_gp = gpar(fontsize = 22), column_title_side = c("bottom"))



mat.lgd <- matrix(c(1, 0), 2, 1)
dimnames(mat.lgd)[[1]] <- c('Rare cell', 'Abundant cell')
ht.legend <- Heatmap(mat.lgd, col = colors, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, 
               row_names_side = c("right"),show_row_names=T,
               row_names_gp = gpar(fontsize=25), 
               rect_gp = gpar(col = "white", lwd = 5),
               layer_fun = function(j, i, x, y, w, h, fill) {
                 ind_mat = restore_matrix(j, i, x, y)
                 for(ir in seq_len(nrow(ind_mat))) {
                   # start from the second column
                   for(ic in seq_len(ncol(ind_mat))) {
                     ind = ind_mat[ir, ic]   # current column
                     v = mat.lgd[i[ind], j[ind]]
                     if(v > 0) { # if they have the same sign
                       col = "steelblue3"
                       grid.points(x[c(ind)], y[c(ind)], 
                                   pch = 19, gp = gpar(col = col), size = unit(10, "mm"))
                     }
                   }
                 }
               })

ComplexHeatmap::draw(ht.legend)




## Draw 68K PBMC tsne Plot for visualization
tsne <- read.table('../10X_full/data/68k_pbmc_barcodes_annotation.tsv', header=T, stringsAsFactors = F, sep='\t')
labels <- read.csv('../10X_full/single-cell-3prime-paper-master/pbmc68k_analysis/PureClusters.csv', header=F, stringsAsFactors = F)
tsne$cluster <- labels$V1
## Make the large plot without legend
p <- ggplot(data=tsne, aes(x=TSNE.1, y=TSNE.2, colour=cluster)) + geom_point(size=0.001) +
  theme_bw() + theme(legend.position='none',  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.text = element_text(size = 16),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank()) + labs(x='tSNE1', y='tSNE2') +
  scale_color_manual(name=NULL, values = c("CD14+ Monocyte"='LightPink', "CD19+ B"="Chocolate1", "CD34+"="Yellow", "CD4+ T Helper2" = "Magenta",
                                           "CD4+/CD25 T Reg"='Snow3', "CD4+/CD45RA+/CD25- Naive T"= "RoyalBlue1", "CD4+/CD45RO+ Memory"= "black",
                                           "CD56+ NK"= 'Orange4', "CD8+ Cytotoxic T"= "Firebrick2", "CD8+/CD45RA+ Naive Cytotoxic"= "CadetBlue2", "Dendritic"="Purple"))
p

# Make the small plot with three cell types without legend
p <- ggplot(data=tsne, aes(x=TSNE.1, y=TSNE.2, colour=cluster)) + geom_point(size=0.001) +
  theme_bw() + theme(legend.position='none',  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.text = element_text(size = 16),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank()) + labs(x=NULL, y=NULL) +
  scale_color_manual(name=NULL, values = c("CD14+ Monocyte"='LightPink', "CD19+ B"="Chocolate1", "CD34+"="white", "CD4+ T Helper2" = "white",
                                           "CD4+/CD25 T Reg"='white', "CD4+/CD45RA+/CD25- Naive T"= "white", "CD4+/CD45RO+ Memory"= "white",
                                           "CD56+ NK"= 'Orange4', "CD8+ Cytotoxic T"= "white", "CD8+/CD45RA+ Naive Cytotoxic"= "white", "Dendritic"="white"))
p


# Make the legend with very large point size
p <- ggplot(data=tsne, aes(x=TSNE.1, y=TSNE.2, colour=celltype)) + geom_point(size=4) +
  theme_bw() + theme(legend.position='right',  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.text = element_text(size = 12),
                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank()) + labs(x='tSNE1', y='tSNE2') +
  scale_color_manual(name=NULL, values = c("CD14+ Monocyte"='LightPink', "CD19+ B"="Chocolate1", "CD34+"="Yellow", "CD4+ T Helper2" = "Magenta",
                                           "CD4+/CD25 T Reg"='Snow3', "CD4+/CD45RA+/CD25- Naive T"= "RoyalBlue1", "CD4+/CD45RO+ Memory"= "black",
                                           "CD56+ NK"= 'Orange4', "CD8+ Cytotoxic T"= "Firebrick2", "CD8+/CD45RA+ Naive Cytotoxic"= "CadetBlue2", "Dendritic"="Purple"))
p














###################################################################################################
## DE_sensitivity_plot
load('results/AUC_result_ALL_2_cells.RData')
load('results/AUC_res_OUR_20200530_log_devide_by_colMeans_density_filtered_gene_2cells.RData')
total <- length(auc.our[[1]])
df <- data.frame(x=rep(1:total, 6),
                 y=c(1-auc.result[[1]], auc.our[[1]], 1-auc.result[[5]], 1- auc.result[[7]], 1-auc.result[[9]], 1-auc.result[[11]]),
                 group=c(rep('RaceID', total), rep('GapClust', total), rep('GiniClust', total), 
                         rep('FiRE(1.5 x IQR)', total), rep('FiRE(1.0 x IQR)', total), rep('FiRE(0.5 x IQR)', total)))

#df <- data.frame(x=rep(1:total, 6),
#                 y=c(1-res.num[[1]], res.num[[3]], 1-res.num[[5]], 1- res.num[[7]], 1-res.num[[9]], 1-res.num[[11]]),
#                 group=c(rep('RaceID', total), rep('GapClust', total), rep('GiniClust', total), 
#                         rep('FiRE(1.5 x IQR)', total), rep('FiRE(1.0 x IQR)', total), rep('FiRE(0.5 x IQR)', total)))
df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)', 'FiRE(1.5 x IQR)'), ordered = T)
sp1 <- ggplot(data=df) + geom_point(aes(x=x, y=y, color=group), size=1)
sp1 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position=c(0.05, 0.85), 
                         legend.text = element_text(size = 14, face = "bold"),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='# of DE genes', y='Area under curve') +
  scale_color_manual(name=NULL, values = c("RaceID"='red', "GapClust"="green4", "GiniClust"="DeepSkyBlue",
                                           "FiRE(1.5 x IQR)" = "steelblue3","FiRE(1.0 x IQR)" = "navy","FiRE(0.5 x IQR)" = "magenta"))

df1 <- df[!df$group %in% c('FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)'),]
df1$group <- as.character(df1$group)
df1$group <- factor(df1$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(1.5 x IQR)'),
                    labels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE'), ordered = T)
sp2 <- ggplot(data=df1) + geom_point(aes(x=x, y=y, color=group), size=1.5)
sp2 + theme_bw() + theme(legend.position=c(0.15, 0.80),  # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.text = element_text(size = 16, face = "bold"), 
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='# of DE genes', y='Area under curve') +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="green4", "GiniClust"="DeepSkyBlue", "FiRE" = "navy")) +
  scale_x_continuous(breaks=c(0, 50, 100), labels=c('0', '50', '100'))




## 5 cells##
load('results/AUC_result_ALL_5_cells.RData')
load('results/AUC_res_OUR_20200530_log_devide_by_colMeans_density_filtered_gene_5cells.RData')
total <- length(auc.our[[1]])
df <- data.frame(x=rep(1:total, 6),
                 y=c(1-auc.result[[1]], auc.our[[1]], 1-auc.result[[5]], 1- auc.result[[7]], 1-auc.result[[9]], 1-auc.result[[11]]),
                 group=c(rep('RaceID', total), rep('GapClust', total), rep('GiniClust', total), 
                         rep('FiRE(1.5 x IQR)', total), rep('FiRE(1.0 x IQR)', total), rep('FiRE(0.5 x IQR)', total)))
df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)', 'FiRE(1.5 x IQR)'), ordered = T)

df1 <- df[!df$group %in% c('FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)'),]
df1$group <- as.character(df1$group)
df1$group <- factor(df1$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(1.5 x IQR)'),
                    labels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE'), ordered = T)
sp2 <- ggplot(data=df1) + geom_point(aes(x=x, y=y, color=group), size=1.5)
sp2 + theme_bw() + theme(legend.position='none', # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.text = element_text(size = 16, face = "bold"), 
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='# of DE genes', y='Area under curve') +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="green4", "GiniClust"="DeepSkyBlue", "FiRE" = "navy")) +
  scale_x_continuous(breaks=c(0, 50, 100), labels=c('0', '50', '100'))


## 10 cells ##
load('results/AUC_result_ALL_10_cells.RData')
load('results/AUC_res_OUR_20200530_log_devide_by_colMeans_density_filtered_gene.RData')
total <- length(auc.our[[1]])
df <- data.frame(x=rep(1:total, 6),
                 y=c(1-auc.result[[1]], auc.our[[1]], 1-auc.result[[5]], 1- auc.result[[7]], 1-auc.result[[9]], 1-auc.result[[11]]),
                 group=c(rep('RaceID', total), rep('GapClust', total), rep('GiniClust', total), 
                         rep('FiRE(1.5 x IQR)', total), rep('FiRE(1.0 x IQR)', total), rep('FiRE(0.5 x IQR)', total)))
df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)', 'FiRE(1.5 x IQR)'), ordered = T)

df1 <- df[!df$group %in% c('FiRE(0.5 x IQR)', 'FiRE(1.0 x IQR)'),]
df1$group <- as.character(df1$group)
df1$group <- factor(df1$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE(1.5 x IQR)'),
                    labels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE'), ordered = T)
sp2 <- ggplot(data=df1) + geom_point(aes(x=x, y=y, color=group), size=1.5)
sp2 + theme_bw() + theme(legend.position='none', # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.text = element_text(size = 16, face = "bold"), 
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='# of DE genes', y='Area under curve') +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="green4", "GiniClust"="DeepSkyBlue", "FiRE" = "navy")) +
  scale_x_continuous(breaks=c(0, 50, 100), labels=c('0', '50', '100'))






# Heatmap for DE and NDE genes
library(pheatmap)
load('DE_NDE_genes.RData')
DE.up <- data3[grepl('DE_up', dimnames(data3)[[1]]),]
DE.down <- data3[grepl('DE_down', dimnames(data3)[[1]]),]
NDE_dat <- data3[grepl('NDE', dimnames(data3)[[1]]),]

drop.NDE <- apply(NDE_dat, 1, function(x){length(x[x>0])/length(x)})
drop.DE.up <- apply(DE.up, 1, function(x){length(x[x>0])/length(x)})
drop.DE.down <- apply(DE.down, 1, function(x){length(x[x>0])/length(x)})
NDE.data <- NDE_dat[order(drop.NDE, decreasing = T)[1:(160 * 3)],]
DE.up <- DE.up[order(drop.DE.up, decreasing = T),]
DE.down <- DE.down[order(drop.DE.down, decreasing = T),]

NDE.data <- NDE.data[,c(97:596, 1:96)]
DE.up <- DE.up[,c(97:596, 1:96)]
DE.down <- DE.down[,c(97:596, 1:96)]

plot.data <- log2(rbind(DE.down, DE.up, NDE.data)+1)

annotation_col = data.frame(
  "CellType" = factor(c(rep("CD19+ B", 500), rep("CD14+ Monocyte", 96)))
)
rownames(annotation_col) = colnames(plot.data)
annotation_row = data.frame(
  "GeneType" = factor(c(rep("DE gene", 160), rep("non DE gene", 480)))
)
rownames(annotation_row) = c(rownames(plot.data))

cols = list(CellType = c("CD19+ B" = "gray80", "CD14+ Monocyte"="firebrick"),
            GeneType = c("DE gene"="green4", "non DE gene"="Chocolate3"))
pheatmap(plot.data, cluster_rows = F, cluster_cols = F, gaps_row = 160, annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_colnames = F, show_rownames = F, annotation_names_col=F, annotation_names_row=F, annotation_colors = cols,
         annotation_legend=F)

NDE.data1 <- NDE.data[,c(1:500, 501:510)]
DE.up1 <- DE.up[,c(1:500, 501:510)]
DE.down1 <- DE.down[,c(1:500, 501:510)]

plot.data1 <- log2(rbind(DE.down1[1,], NDE.data1[2:dim(NDE.data1)[1],]) + 1)
pheatmap(plot.data1, cluster_rows = F, cluster_cols = F, annotation_col = annotation_col,
         show_colnames = F, show_rownames = F, annotation_names_col=F, annotation_names_row=F, annotation_colors = cols,
         annotation_legend=F)

plot.data2 <- log2(rbind(DE.down1[1:40,], DE.up1[1:40,], NDE.data1[81:dim(NDE.data1)[1],]) + 1)
pheatmap(plot.data2, cluster_rows = F, cluster_cols = F, annotation_col = annotation_col,
         show_colnames = F, show_rownames = F, annotation_names_col=F, annotation_names_row=F, annotation_colors = cols,
         annotation_legend=F)

plot.data3 <- log2(rbind(DE.down1[1:80,], DE.up1[1:80,], NDE.data1[161:dim(NDE.data1)[1],]) + 1)
pheatmap(plot.data3, cluster_rows = F, cluster_cols = F, annotation_col = annotation_col,
         show_colnames = F, show_rownames = F, annotation_names_col=F, annotation_names_row=F, annotation_colors = cols,
         annotation_legend=F)





