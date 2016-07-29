################################################
####  PLOT DELLE ROH DEL PROGRAMMA BISCARINI ###
################################################

#plot 
library(ggplot2)
razza<-read.table(file='detectRUNS.ROHom.csv',header=T,sep=',')
head(razza)

names(razza) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

pdf('nome_plot3.pdf',height=12, width=20)

for (a in sort(unique(razza$CHROMOSOME))){
  
  print(paste('Chromosome: ',a))
  
  #primo subset

  cromo<-subset(razza,CHROMOSOME==a)
  #secondo subset
  sotto<-cromo[,c(5,6,2,1)]
  sotto=sotto[order(sotto$POPULATION),]
  #creo un numero progressivo per animale
  newID=seq(1,length(unique(sotto$IND)))
  id=unique(sotto$IND)
  sotto$NEWID=newID[match(sotto$IND,id)]
  
  #colore delle razze
  colore=as.factor(sotto$POPULATION)
  
  #lughezza in mb
  sotto[,1]<-sotto[,1]/1000000
  sotto[,2]<-sotto[,2]/1000000
  
  #PRIMO PLOT
  grafico=ggplot(sotto,aes(colour=sotto$POPULATION))  + 
    geom_segment(data=sotto,aes(x = START, y = IND, xend = END, yend = IND),alpha=1,size=1)+
    xlim(0, max(sotto$END)) + 
    ggtitle(paste('Cromosome:',a))
  print(grafico)
  
  #SECONDO PLOT
  grafico1=ggplot(sotto,aes(colour=sotto$POPULATION))  + 
    geom_segment(data=sotto,aes(x = START, y = NEWID, xend = END, yend = NEWID),alpha=1,size=1)+
    xlim(0, max(sotto$END)) + 
    ggtitle(paste('Cromosome:',a))
  
  print(grafico1)
  
  #TERZO PLOT
  sotto$JitterSpecies <- ave(as.numeric(sotto$POPULATION), sotto$POPULATION,FUN = function(x) x + rnorm(length(x), sd = .1))
  grafico2 =  ggplot(sotto, aes(x = START, xend = END, y = JitterSpecies, yend = JitterSpecies)) +
    geom_segment() + xlim(0, max(sotto$END)) + 
    scale_y_continuous("Species", breaks = seq(unique(sotto$POPULATION)), labels = levels(sotto$POPULATION)) +
    ggtitle(paste('Cromosome:',a))
  
  
  print(grafico2)
  
}

dev.off()
