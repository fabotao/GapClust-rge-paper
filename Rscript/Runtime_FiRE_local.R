homedir = '/home/sam/Documents/FBT/Single/package/GapClust_manuscript'
workdir              = paste0(homedir, '/Rproj/Runtime/')  
# where you put the data and results

setwd(workdir)

require('Matrix')
require('plyr')

ranger_preprocess<-function(data_mat, ngenes_keep=1000, dataSave='./', optionToSave=F, minLibSize=0, verbose=T){
  
  ngenes_keep = 1000
  #write.csv(x = data_mat$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
  #l<-.normalize_by_umi(data_mat)
  l<-.normalize_by_umi_2(data_mat, dataSave, minLibSize, verbose)
  m_n<-l$m
  
  if (verbose){
    cat("Select variable Genes...\n")
  }
  
  df<- .get_variable_gene(m_n)
  gc()
  
  if (verbose){
    cat("Sort Top Genes...\n")
  }
  
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  
  if (verbose){
    cat("Cutoff Genes...\n")
  }
  
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  #system("rm genes",ignore.stderr = T)
  
  if (optionToSave){
    write.csv(features, file = paste(dataSave,"genes",sep=""), quote = F,row.names = F)
    write.csv(l$use_genes[features], file = paste(dataSave,"genes_used_all",sep=""), quote = F,row.names = F)
  }
  
  #genes = read.csv(file = "genes")
  #features = genes$x
  
  #Final Data
  m_n_68K<-m_n[,features]
  m_filt<-Matrix(log2(m_n_68K+1),sparse = T)
  
  if (verbose){
    cat(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2]))
  }
  #system("rm whole_matrix",ignore.stderr = T)
  if (optionToSave){
    writeMM(m_filt,file="whole_matrix")
  }
  
  list(preprocessedData=m_filt, selGenes=features)
}

# ---------------------------------------------
# normalize the gene barcode matrix by umi
# filter based on read count first
# ---------------------------------------------

.normalize_by_umi_2 <-function(x, dataSave, minLibSize, verbose){
  mat  = x$mat
  gene_symbols = x$gene_symbols
  
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
    write.csv(keepCellIndex, file = paste(dataSave,"keepCellIndexes.csv",sep=""), quote = F,row.names = F)
  }
  
  #Filter Genes
  cs <- colSums(mat>2)
  x_use_genes <- which(cs > 3)
  
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


# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m){
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}


load('../10X_full/data/PBMC_68K_normalized.RData')


cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

for (i in 1:length(cell.num)){
  set.seed(cell.num[i] + i)
  data1 <- data4[,sample(1:dim(data4)[2], cell.num[i])]
  data1 <- data1[rowMeans(data1) > 0,]

  data <- t(data1) #Samples * Features
  
  start.time <- Sys.time()
  #Genes
  genes <- c(1:dim(data)[2]) #It can be replaced with original gene names
  
  data_mat <- list(mat=data, gene_symbols=genes)
  preprocessedList <- ranger_preprocess(data_mat)
  preprocessedData <- as.matrix(preprocessedList$preprocessedData)
  model <- new(FiRE::FiRE, 100, 50, 1017881, 5489, 0)
  model$fit(preprocessedData)
  score <- model$score(preprocessedData)
  
  #Apply IQR-based criteria to identify rare cells for further downstream analysis.
  q3 <- quantile(score, 0.75)
  iqr <- IQR(score)
  th <- q3 + (1.5*iqr)
  
  #Select indexes that satisfy IQR-based thresholding criteria.
  indIqr <- which(score >= th)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  save(time.taken, file=paste0('FiRE_runtime_', i, '.RData'))
  
  
  #Create a file with binary predictions
  predictions <- integer(dim(data)[1])
  predictions[indIqr] <- 1 #Replace predictions for rare cells with '1'.
  
  if(i == 68579){
    save(predictions, file=paste0('RaceID_68K.RData'))
  }
  rm(list=c('start.time', 'time.taken', 'end.time', 'predictions', 'data1', 'data', 
            'indIqr', 'q3', 'iqr', 'score', 'model', 'preprocessedData', 'preprocessedList',
            'data_mat', 'genes'))
  gc()
  
}


