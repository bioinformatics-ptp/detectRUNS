################################################
####  PLOT DELLE ROH DEL PROGRAMMA BISCARINI ###
################################################

#plot
library("ggplot2")
runs <- read.table(file='detected.ROHom.csv',header=T,sep=';')
head(runs)

names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

runs$POPULATION <- as.character(runs$POPULATION)
runs[runs$IND %in% sample(unique(runs$IND),6),"POPULATION"] <- "SWS"
write.table(runs,file="detected.ROHom.csv",quote=FALSE,row.names=FALSE,col.names=TRUE,sep=";")

plotRuns <- function(runsFile = 'detected.ROHet.csv', suppressInds = FALSE) {
  
  
  runs <- read.table(file=runsFile, header=TRUE, sep=';')
  
  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")
  
  for (chrom in sort(unique(runs$CHROMOSOME))) {
    
    #subset by chromosome
    krom <- subset(runs,CHROMOSOME==chrom)
    
    #rearrange subset
    teilsatz <- krom[,c(5,6,2,1)]
    teilsatz <- teilsatz[order(teilsatz$POPULATION),]
    
    #new progressive numerical ID
    newID <- seq(1,length(unique(teilsatz$IND)))
    id <- unique(teilsatz$IND)
    teilsatz$NEWID=newID[match(teilsatz$IND,id)]
    
    optionen <- scale_y_discrete("IDs",limits=unique(teilsatz$IND))
    alfa <- 1
    grosse <- 1
    
    if (length(id) > 50) {
      
      optionen <- theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
      alfa <- 0.75
      grosse <- 0.25
    }
    
    if (suppressInds) optionen <- theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
      
    #lughezza in mb
    teilsatz$START <- (teilsatz$START/(10^6))
    teilsatz$END <- (teilsatz$END/(10^6))
    
    row.names(teilsatz) <- NULL
    
    teilsatz$IND <- as.factor(teilsatz$IND)
    teilsatz$IND <- factor(teilsatz$IND, levels = unique(teilsatz$IND[order(teilsatz$NEWID)]))
    
    titel <- paste(unlist(strsplit(runsFile,"\\."))[2],"chromosome",chrom,sep="_")
    
    p <- ggplot(teilsatz)
    p <- p + geom_segment(data=teilsatz,aes(x = START, y = IND, xend = END, yend = IND,colour=as.factor(POPULATION)), alpha=alfa, size=grosse)
    p <- p + xlim(0, max(teilsatz$END)) + ggtitle(paste('Chromosome:',chrom)) 
    p <- p + guides(colour=guide_legend(title="Population")) + xlab("Mbps")
    p <- p + optionen
    
    pdf(paste(titel,".pdf",sep=""),height=8,width=10)
    print(p)
    dev.off()
  }
}

plotRuns("detected.ROHom_x.csv")
plotRuns("detected.ROHom.csv")
plotRuns("detected.ROHet.csv")

p <- ggplot(teilsatz)
p <- p + geom_segment(data=teilsatz,aes(x = START, y = IND, xend = END, yend = IND,colour=as.factor(POPULATION)), alpha=1, size=1)
p <- p + xlim(0, max(teilsatz$END))
p <- p + ggtitle(paste('Chromosome:',chrom)) 
p <- p + guides(colour=guide_legend(title="Population")) + xlab("Mbps")
# p <- p + scale_y_discrete("IDs",limits=unique(teilsatz$IND))
p <- p + theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
print(p)


titel <- "detected.ROHet.csv"
