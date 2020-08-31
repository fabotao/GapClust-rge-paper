library(RaceID)
homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/10X_subsample/')  
# where you put the data and results

setwd(workdir)

RaceID.result <- list()
for (i in 2:100){
  #k=ks[j]
  exprimentID<-paste("10X_rare_",i,"_CD14_cells",sep="")
  load(paste0('results/', exprimentID, '_ExprM.filter.RData'))

  data1 <- ExprM.normCounts.filter[rowMeans(ExprM.normCounts.filter) > 0,]
  sc <- SCseq(data1)
  sc <- filterdata(sc,mintotal=min(colSums(data1))-1)
  fdata <- getfdata(sc)
  sc <- compdist(sc,metric="pearson")
  sc <- clustexp(sc)
  #plotsaturation(sc,disp=FALSE)
  #plotsaturation(sc,disp=TRUE)
  #plotjaccard(sc)
  sc <- clustexp(sc,cln=2,sat=FALSE)
  sc <- findoutliers(sc)
  #plotbackground(sc)
  #plotsensitivity(sc)
  #plotoutlierprobs(sc)
  #clustheatmap(sc)
  finalcls <- (sc@cpart)
  print(table(finalcls))
  RaceID.result[[i-1]] <- finalcls
  
}

save(RaceID.result, file='results/RaceIDClustering.RData')


