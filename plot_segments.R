################################################
####  PLOT DELLE ROH DEL PROGRAMMA BISCARINI ###
################################################

#plot 
library(ggplot2)
razza<-read.table(file='step1',header=T,sep=';')
head(razza)

pdf('nome_plot3.pdf',height=12, width=20)

for (a in sort(unique(razza$CHROMOSOME))){
  
  print(paste('Chromosome: ',a))
  
  #primo subset
  cromo<-subset(razza,CHROMOSOME==a)
  #secondo subset
  sotto<-cromo[,c(5,6,2,1)]
  sotto=sotto[order(sotto$BREED),]
  #creo un numero progressivo per animale
  newID=seq(1,length(unique(sotto$ANIMAL)))
  id=unique(sotto$ANIMAL)
  sotto$NEWID=newID[match(sotto$ANIMAL,id)]
  
  #colore delle razze
  colore=as.factor(sotto$BREED)
  
  #lughezza in mb
  sotto[,1]<-sotto[,1]/1000000
  sotto[,2]<-sotto[,2]/1000000
  
  #PRIMO PLOT
  grafico=ggplot(sotto,aes(colour=sotto$BREED))  + 
    geom_segment(data=sotto,aes(x = START, y = NEWID, xend = END, yend = NEWID),alpha=1,size=0.5)+
    xlim(0, max(sotto$END)) + 
    ggtitle(paste('Cromosome:',a))
  
  print(grafico)
  
  #SECONDO PLOT
  sotto$JitterSpecies <- ave(as.numeric(sotto$BREED), sotto$BREED,FUN = function(x) x + rnorm(length(x), sd = .1))
  grafico1 =  ggplot(sotto, aes(x = START, xend = END, y = JitterSpecies, yend = JitterSpecies)) +
    geom_segment()+
    scale_y_continuous("Species", breaks = seq(unique(sotto$BREED)), labels = levels(sotto$BREED)) +
    ggtitle(paste('Cromosome:',a))
  
  
  print(grafico1)
  
}

dev.off()
