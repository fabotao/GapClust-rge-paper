library(ggplot2)
library(ComplexHeatmap)
library(splatter)
workdir              = "/home/sam/Documents/RareCellDetection/Proj/jurkat_two_species/"     
# where you put the data and results

setwd(workdir)


load('results/RaceIDClustering.RData') # RaceID.result
load('results/OURClustering.RData') # GapClust.result
load('results/GiniClustering.RData') # Gini.result
load('results/FiREClustering.RData') # FiRE.result

true.res <- RaceID.result
jurkat.num.vec <- c(8, 16, 24, 32, 40, 48, 56, 65, 73, 82)
for (jurkat.num in jurkat.num.vec){
  #k=ks[j]
  exprimentID<-paste("Jurkat_",jurkat.num,"_293_1580_cells",sep="")
  load(paste0('results/', exprimentID, '.RData'))
  true.res[[which(jurkat.num.vec %in% jurkat.num)]] <- ifelse(grepl('Jurkat', dimnames(jurkat.293)[[2]]), 1, 0)
  rm(list=c('jurkat.293', 'exprimentID'))
}

Gini.result <- sapply(Gini.result, as.character)

transform.gini <- function(){
  gini.res <- list()
  for(i in 1:length(Gini.result)){
    cnt <- table(Gini.result[[i]])
    cnt.top <- table(Gini.result[[i]][1541:(1540 + jurkat.num.vec[i])])
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
    #cnt.top <- table(RaceID.result[[i]][1541:(1540 + jurkat.num.vec[i])])
    #small <- names(cnt)[cnt <= sum(cnt) * 0.15] # 0.15 is for 100 cells or a little more
    #small <- small[small %in% names(cnt.top)]
    major <- names(cnt)[which.max(cnt)]
    race.res[[i]] <- ifelse(RaceID.result[[i]] %in% major, 0, 1)
    # if(length(small)>=1){
    #   #race.res[[i]] <- ifelse(RaceID.result[[i]] %in% small, 1, 0)
    #   race.res[[i]] <- ifelse(RaceID.result[[i]] %in% major, 0, 1)
    # }else{
    #   race.res[[i]] <- (rep(0, length(RaceID.result[[i]])))
    # }
  }
  return(race.res)
}

gini.res <- transform.gini()
race.res <- transform.race()

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


df <- data.frame(x=rep(c(8, 16, 24, 32, 40, 48, 56, 65, 73, 82), 4),
                 y=c(f3, f2, f1, f4),
                 group=c(rep('RaceID', 10), rep('GapClust', 10), rep('GiniClust', 10), rep('FiRE', 10)))

df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'FiRE'), ordered = T)
sp1 <- ggplot(data=df) + geom_point(aes(x=x, y=y, color=group), size=3) + geom_line(aes(x=x, y=y, color=group)) 
sp1 + theme_bw() + theme(#panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                         legend.position='right', panel.border = element_rect(colour = "black", fill=NA, size=1),
                         legend.text = element_text(size = 14, face = "bold"),
                         axis.text.x = element_text(size=14),
                         axis.text.y = element_text(size=14),
                         axis.title.x = element_text(size=16),
                         axis.title.y = element_text(size=16),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_blank()) + labs(x='% of rare cells', y=expression(paste(F["1 "], score))) +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="forestgreen", "GiniClust"="DeepSkyBlue", "FiRE" = "navy")) +
  scale_x_continuous(breaks=c(8, 16, 24, 32, 40, 48, 56, 65, 73, 82), labels=c('0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '4.5', '5.0'))



