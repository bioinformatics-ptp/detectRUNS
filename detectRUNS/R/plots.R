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
#' @param ROHet shall we detect ROHet or ROHom? (TRUE/FALSE)
#'
#' @return converted vector of genotypes
#' @export
#'
#' @examples #not yet
#'
#'
#plot

plotRuns <- function(runsFile = 'detectRUNS.ROHet.csv', ROHet = TRUE) {


  runs <- read.table(file=runsFile, header=TRUE, sep=',')

  names(runs) <- c("POPULATION","IND","CHROMOSOME","START","END","COUNT","LENGTH")

  #divisione delle razze
  #scelta razza e cromosoma DEVI CAMBIARE IL NOME DEL FILE PDF SE LO VUOI E LA RAZZA
  kleur <- as.factor(runs$POPULATION)

  for (krom in seq(as.factor(table(runs$CHROMOSOME)))) {

    cromo <- subset(runs,CHROMOSOME==krom)
    sottoInsieme <- cromo[,c(4,5,2)]

    #lughezza in mb
    sottoInsieme[,1] <- sottoInsieme[,1]/1000000
    sottoInsieme[,2] <- sottoInsieme[,2]/1000000

    runType <- ifelse(ROHet,"ROHom","ROHet")
    titel <- paste(runType,"chromosome",krom,sep="_")

    p <- ggplot() + geom_segment(data=sottoInsieme, aes(x = START, y = IND, xend = END, yend = IND), colour="red", alpha=1 ,size=4) +
      xlim(0, max(sottoInsieme$END)) + ggtitle(titel)

    pdf(paste(titel,".pdf",sep=""),height=8,width=10)
    print(p)
    dev.off()
  }
}




