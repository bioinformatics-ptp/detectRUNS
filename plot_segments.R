################################################
####  PLOT DELLE ROH DEL PROGRAMMA BISCARINI ###
################################################

#plot 
library("ggplot2")
razza<-read.table(file='detectRUNS.ROHet.csv',header=T,sep=',')
head(razza)

names(razza) <- c("ANIMAL","CHROMOSOME","START","END","COUNT","LENGTH")
razza$BREED <- rep("CHIL",nrow(razza))
razza <- razza[,c(7,1,2,3,4,5,6)]

#divisione delle razze
#scelta razza e cromosoma DEVI CAMBIARE IL NOME DEL FILE PDF SE LO VUOI E LA RAZZA
colore=as.factor(razza$BREED)

#pdf('reads_ROHet_BISCARINI.pdf',height=12, width=20)

for (a in seq(as.factor(table(razza$CHROMOSOME)))){
  cromo<-subset(razza,CHROMOSOME==a)
  sotto<-cromo[,c(4,5,2)]
  
  #lughezza in mb
  sotto[,1]<-sotto[,1]/1000000
  sotto[,2]<-sotto[,2]/1000000
  
  grafico=ggplot()  + geom_segment(data=sotto,aes(x = START, y = ANIMAL, xend = END, yend = ANIMAL),colour="red",alpha=1,size=4)+
      xlim(0, max(sotto$END)) + ggtitle("Biscarini's method")
  
  print(grafico)
  
}

#dev.off()



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


#############################################
####  PLOT DELLE ROH DEL PROGRAMMA MARRAS ###
#############################################
#preparazione dei dati
library(ggplot2)

razza<-read.table(file='example',header=T,sep=';')
head(razza)

#divisione delle razze
#scelta razza e cromosoma DEVI CAMBIARE IL NOME DEL FILE PDF SE LO VUOI E LA RAZZA
#controllare se funziona!!!!!
colore=as.factor(razza$BREED)


#pdf('reads_ROHet_MARRAS.pdf',height=12, width=20)

seq(as.factor(table(razza$CHROMOSOME)))
for (a in seq(as.factor(table(razza$CHROMOSOME)))){
  cromo<-subset(razza,CHROMOSOME==a)
  sotto<-cromo[,c(5,6,2)]
  
  #lughezza in mb
  sotto[,1]<-sotto[,1]/1000000
  sotto[,2]<-sotto[,2]/1000000
  
  grafico=ggplot()  + geom_segment(data=sotto,aes(x = START, y = ANIMAL, xend = END, yend = ANIMAL),colour="red",alpha=1,size=4)+
    xlim(0, max(sotto$END)) + ggtitle("Marras's method")
  
  print(grafico)
  
}

#dev.off()

