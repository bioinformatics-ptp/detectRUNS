##################################
### Compute ROH in R (a la Plink)
##################################

#' Run functions to detect RUNS (ROHom/ROHet)
#'
#' This is the main function of RUNS and would probably require
#' some good documentation.
#'
#' @param gegevens genotype dataset in the Plink .raw format
#' @param mappa Plink map file
#' @param window the size of sliding window
#' @param drempel the threshold of overlapping windows of the same state (homozygous/heterozygous) to call a SNP in a RUN
#' @param minSNP minimum n. of SNP in a RUN
#' @param ROHet should we look for ROHet or ROHom?
#'
#' @return n. of individuals for which runs have been written out
#' @export
#'
#' @import plyr
#' @import itertools
#' @import ggplot2
#'
#' @examples x <- RUNS.run(gegevens, mappa, windowSize = 20, drempel = 0.1, minSNP = 5, ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1)
#'

#library("plyr")

RUNS.run <- function(gegevens, mappa, windowSize = 15, drempel = 0.1, minSNP = 3, ROHet = TRUE,
                     maxOppositeGenotype = 1, maxMiss = 1) {

  #gegevens <- read.table("RoHet/DATA/subsetChillingham.raw",header=TRUE)
  #remove unnecessary fields from the .raw file
  gegevens <- gegevens[,-c(3,4,5,6)]

  report_filename <- paste("detectRUNS",ifelse(ROHet,"ROHet","ROHom"),"csv",sep=".")
  if(file.exists(report_filename)) system2("rm",report_filename)

  zustand <- vector()

  # require "plyr"
  staat <- daply(gegevens,"IID",function(x) {

    y <- schiebeFenster(as.integer(x[-c(1,2)]),windowSize,step=1,ROHet=ROHet,maxOppositeGenotype,maxMiss);
    snpRun <- snpInRun(y,windowSize,drempel)
    dRUN <- createRUNdf(snpRun,mappa,minSNP)
    zustand <- schreibRUN(as.character(x$IID),dRUN,ROHet,as.character(x$FID))
    return(zustand)
  })

  return(sum(staat))
}






