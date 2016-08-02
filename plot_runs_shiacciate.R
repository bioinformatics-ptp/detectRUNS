###############################
#### PLOT DELLE RUNS VICINE ###
###############################

#READ a FILE
runs<-read.table(file='detected.ROHom.csv',header=T,sep=';')
dati<-read.table(file='detected.ROHet_x.csv',header=T,sep=';')
head(dati)

names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

pdf('name_file1.pdf',height=12, width=20)

#select a POPULATION
for (rasse in unique(runs$POPULATION)){
  print(paste('Subset for ',rasse))
  teilsatz <- subset(runs,runs$POPULATION==rasse)
  
  #select a chromosome
  for (chrom in sort(unique(teilsatz$CHROMOSOME))){ 
    
    print (paste('CHR: ',chrom))
    krom <- subset(teilsatz,CHROMOSOME==chrom)
    krom <- krom[order(krom$START),]
    
    #start the order 
    yread <- c(); #keeps track of the x space that is used up by segments 
    
    # get x axis limits
    minstart <- min(krom$START);
    maxend <- max(krom$END);
    
    # initialise yread
    yread[1] <- minstart - 1;
    ypos <- c(); #holds the y pos of the ith segment
    
    #sorted<-subset(sorted,CHROMOSOME==1)
    
    for (r in 1:nrow(sorted)){
      read <- krom[r,];   
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
    krom$ypos <- ypos;
    head(krom)
    
    #PLOT TIME
    grafico=ggplot()  + 
      geom_segment(data=krom,aes(x = START/1000000, y = ypos, xend = END/1000000, yend = ypos),colour="blue",alpha=1,size=1)+
      xlim(0, max(sorted$END/1000000)+10) + 
      ylab('n Runs') + xlab('Chromosome position (Mbps)') +
      ggtitle(paste("POPULATION: ",bre,'\nCromosome:',a))
    
    print(grafico)
    
  }
}

dev.off()