library(ggplot2)
library(ComplexHeatmap)
library(splatter)
workdir              = "/home/sam/Documents/RareCellDetection/Proj/Doublet/"     
# where you put the data and results

setwd(workdir)


load('results/RaceIDClustering.RData') # RaceID.result
load('results/OURClustering.RData') # GapClust.result
load('results/GiniClustering.RData') # Gini.result
load('results/FiREClustering.RData') # FiRE.result
load('results/CELLSIUSClustering.RData') # FiRE.result

true.res <- RaceID.result
for (q in 1:7){
  for (i in 1:20){
    #set.seed(seeds[(q-1)*20+i])
    exprimentID<-paste("10X_rare",q,"_",i, '_ExprM.filter.RData',sep="")
    load(paste("data/", exprimentID, sep=""))
    true.res[[((q-1) * 20) + i]] <- ifelse(grepl('CD14', dimnames(ExprM.RawCounts.filter)[[2]]), 1, 0)
  }
}


Gini.result <- sapply(Gini.result, as.character)

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


transform.cell <- function(){
  cell.res <- list()
  for(i in 1:length(CELLSIUS.result)){
    cell.res[[i]] <- (rep(0, length(CELLSIUS.result[[i]])))
    cnt.top <- table(CELLSIUS.result[[i]][1:2])
    cnt.top <- cnt.top[cnt.top > 0 & grepl('_', names(cnt.top))]
    cnt <- table(CELLSIUS.result[[i]])
    if(length(cnt.top)>=1){
      for(j in 1:length(cnt.top)){
        label <- names(cnt.top)[j]
        if(cnt[label] < sum(cnt) * 0.15){
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


for(i in 1:10){
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
f5 <- calcul.f1(cell.res)

group <- c(rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20), rep(5, 20), rep(6, 20), rep(7, 20))
f1.ave <- tapply(f1, group, mean)
f2.ave <- tapply(f2, group, mean)
f3.ave <- tapply(f3, group, mean)
f4.ave <- tapply(f4, group, mean)
f5.ave <- tapply(f5, group, mean)

df <- data.frame(x=as.numeric(as.character(c(names(f1.ave), names(f2.ave), names(f3.ave), names(f4.ave), names(f5.ave)))),
                 y=c(f1.ave, f2.ave, f3.ave, f4.ave, f5.ave),
                 group=c(rep('GiniClust', 7), rep('GapClust', 7), rep('RaceID', 7),  rep('FiRE', 7), rep('CellSIUS', 7)))

df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'CellSIUS', 'FiRE'), ordered = T)
sp1 <- ggplot(data=df) + geom_point(aes(x=x, y=y, color=group), size=3) + geom_line(aes(x=x, y=y, color=group)) 
sp1 + theme_bw() + theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position='right', panel.border = element_rect(colour = "black", fill=NA, size=1),
                         legend.text = element_text(size = 14, face = "bold"),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='% of rare cells', y=expression(paste(F["1 "], score))) +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="green4", "GiniClust"="DeepSkyBlue", "CellSIUS" = "Magenta", "FiRE" = "navy")) +
  scale_x_continuous(breaks=sort(unique(df$x)), labels=c('0.08', '0.16', '0.33', '0.66', '1.32', '2.60', '5.00'))


