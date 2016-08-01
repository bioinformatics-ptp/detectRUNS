###############################
#### PLOT DELLE RUNS VICINE ###
###############################

#READ a FILE
dati<-read.table(file='detected.ROHom.csv',header=T,sep=';')
dati<-read.table(file='detected.ROHet.csv',header=T,sep=';')
head(dati)

names(dati) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

pdf('name_file1.pdf',height=12, width=20)

#select a POPULATION
for (bre in unique(dati$POPULATION)){
  print(paste('Subset for ',bre))
  razza=subset(dati,dati$POPULATION==bre)
  
  #select a chromosome
  for (a in sort(unique(razza$CHROMOSOME))){ 
    
    print (paste('CHR: ',a))
    cromo<-subset(razza,CHROMOSOME==a)
    sorted=cromo[order(cromo$START),]
    
    #start the order 
    yread <- c(); #keeps track of the x space that is used up by segments 
    
    # get x axis limits
    minstart <- min(sorted$START);
    maxend <- max(sorted$END);
    
    # initialise yread
    yread[1] <- minstart - 1;
    ypos <- c(); #holds the y pos of the ith segment
    
    #sorted<-subset(sorted,CHROMOSOME==1)
    
    for (r in 1:nrow(sorted)){
      read <- sorted[r,];   
      start <- read$START;
      placed <- FALSE;
      
      # iterate through yread to find the next availible
      # y pos at this x pos (start)
      y <- 1;
      while(!placed){
        
        if(yread[y] < start){
          ypos[r] <- y;
          yread[y] <- read$END;
          placed <- TRUE;
        } 
        
        # current y pos is used by another segment, increment
        y <- y + 1;
        # initialize another y pos if we're at the end of the list
        if(y > length(yread)){
          yread[y] <- minstart-1;
        }
      }
    } 
    
    maxy <- length(yread);
    sorted$ypos <- ypos;
    head(sorted)
    
    #PLOT TIME
    grafico=ggplot()  + 
      geom_segment(data=sorted,aes(x = START/1000000, y = ypos, xend = END/1000000, yend = ypos),colour="blue",alpha=1,size=1)+
      xlim(0, max(sorted$END/1000000)+10) + 
      ylab('n Runs') + xlab('Chromosome position (Mbps)') +
      ggtitle(paste("POPULATION: ",bre,'\nCromosome:',a))
    
    print(grafico)
    
  }
}

dev.off()