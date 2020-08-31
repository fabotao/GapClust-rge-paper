setwd('/home/sam/Documents/RareCellDetection/Proj/Runtime/')

library(ggplot2)

## Memory usage plot
memory.dat <- data.frame(method=c('GiniClust', 'RaceID', 'CellSIUS','FiRE', 'GapClust'),
                         size=c(117.93, 78.58, 168.16, 19.47, 2.79),
                         col=c('DeepSkyBlue', 'red1', 'Magenta', 'navy', 'green4'))
memory.dat$method <- factor(memory.dat$method, levels=c('GapClust', 'RaceID', 'GiniClust', 'CellSIUS','FiRE'), ordered = T)

col=c('green4', 'red1', "DeepSkyBlue", "Magenta", "navy")

p1 <- ggplot(memory.dat,aes(method, size, fill=method)) + geom_bar(stat='identity', width=0.5) + 
  scale_fill_manual(name=NULL, values = col) + theme_bw() +
  theme(legend.position="none", panel.border = element_rect(color='black', fill=NA, size=1),
        axis.text.x = element_text(size=14, colour = 'black'),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)
  ) + labs(x=' ', y='Memory usage/GB')

#write.csv(memory.dat, file='memory.csv', fileEncoding = 'utf-8')

## Computation time


cell.num <- c(1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 68579)

cellsius.time <- c()
for(i in 1:14){
  load(paste0('CELLSIUS_runtime_', i, '.RData'))
  cellsius.time <- c(cellsius.time, as.numeric(time.taken, units='mins'))
}


fire.time <- c()
for(i in 1:16){
  load(paste0('FiRE_runtime_', i, '.RData'))
  fire.time <- c(fire.time, as.numeric(time.taken, units='mins'))
}


our.time <- c()
for(i in 1:16){
  load(paste0('OUR_runtime_', i, '.RData'))
  our.time <- c(our.time, as.numeric(time.taken, units='mins'))
}

gini.time <- c()
for(i in 1:11){
  load(paste0('Gini_runtime_', i, '.RData'))
  gini.time <- c(gini.time, as.numeric(time.taken, units='mins'))
}



race.time <- c()
for(i in 1:9){
  load(paste0('RaceID_runtime_', i, '.RData'))
  race.time <- c(race.time, as.numeric(time.taken, units='mins'))
}


df <- data.frame(x=c(cell.num[1:14], cell.num[1:16], cell.num[1:16], cell.num[1:11], cell.num[1:9]),
                 y=c(cellsius.time, fire.time, our.time, gini.time, race.time),
                 group=c(rep('CellSIUS', 14), rep('FiRE', 16), rep('GapClust', 16), rep('GiniClust', 11), rep('RaceID', 9)))
df$logx <- log10(df$x)
df$logy <- log10(df$y)

df$group <- factor(df$group, levels=c('GapClust', 'RaceID', 'GiniClust', 'CellSIUS','FiRE'), ordered = T)
sp1 <- ggplot(data=df) + geom_point(aes(x=logx, y=logy, color=group), size=3) + geom_line(aes(x=logx, y=logy, color=group))
sp1 + theme_bw() + theme(#panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
  legend.position=c(0.15, 0.83), panel.border = element_rect(colour = "black", fill=NA, size=1),
  legend.text = element_text(size = 14, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.text.y = element_text(size=14),
  axis.title.x = element_text(size=16),
  axis.title.y = element_text(size=16),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank()) + labs(x=expression(paste(Log['10 '], '(#samples)')), y=expression(paste(Log["10 "], '(minutes)'))) +
  scale_color_manual(name=NULL, values = c("RaceID"='red1', "GapClust"="forestgreen", "GiniClust"="DeepSkyBlue", "FiRE" = "navy", 'CellSIUS'='magenta'))


#write.csv(df, file='runtime.csv', fileEncoding = 'utf-8')





