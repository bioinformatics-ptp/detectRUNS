####################################
## FUNCTIONS TO MAKE PLOTS FROM RUNS
####################################

# three functions for three different plots:
# i) Runs per animal (animals on the y-axis)
# ii) stacked runs per animal
# iii) n. of times a SNP is in a run in the population

#required external packages
#library("ggplot2")

#' Function to plot runs per animal (see Williams et al. 2016, Animal Genetics)
#' IDs on the y-axis, bps on the x-axis: plots run (TRUE) / no run (FALSE)
#'
#' @param runs output file with runs per animal (breed, id, chrom, nSNP, start, end, length) #defaults to detectRUNS.ROHet.csv
#' @param suppressInds shall we suppress individual IDs on the y-axis? (defaults to FALSE)
#'
#' @return plot of runs by chromosome (pdf files)
#' @export
#'
#' @examples #not yet
#'
#'
#plot

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





