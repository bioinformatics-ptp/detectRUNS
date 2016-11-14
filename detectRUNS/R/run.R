##################################
### Compute ROH in R (a la Plink)
##################################

#' Run functions to detect RUNS (ROHom/ROHet)
#'
#' This is the main function of RUNS and would probably require
#' some good documentation.
#'
#' @param genotype genotype dataset in the Plink .raw format
#' @param mapFile Plink map file (either the file path/name or the R data.frame)
#' @param windowSize the size of sliding window
#' @param threshold the threshold of overlapping windows of the same state (homozygous/heterozygous) to call a SNP in a RUN
#' @param minSNP minimum n. of SNP in a RUN
#' @param ROHet should we look for ROHet or ROHom?
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#'
#' @return n. of individuals for which runs have been written out
#' @export
#'
#' @import plyr
#' @import itertools
#' @import ggplot2
#' @import itertools
#' @import utils
#'
#' @examples
#' data(chillingam)
#' x <- RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000, minDensity = 1/10)
#'

RUNS.run <- function(genotype, mapFile, windowSize = 15, threshold = 0.1, minSNP = 3, ROHet = TRUE,
                     maxOppositeGenotype = 1, maxMiss = 1, maxGap = 10^6, minLengthBps = 1000, minDensity = 1/10) {

  if(!is.data.frame(genotype)) {

    if(file.exists(genotype)){
      genotype <- utils::read.table(genotype,header=TRUE)
    }
  }

  if(!is.data.frame(mapFile)) {

    if(file.exists(mapFile)){
      mapFile <- utils::read.table(mapFile)
    }
  }

  names(mapFile) <- c("Chrom","SNP","cM","bps")

  #remove unnecessary fields from the .raw file
  genotype <- genotype[,-c(3,4,5,6)]

  ## write out populations/individuals for further plots (snpInRuns)
  utils::write.table(genotype[,c(1,2)],file="genotype.raw",quote=FALSE,row.names=FALSE,col.names=TRUE)

  report_filename <- paste("detected",ifelse(ROHet,"ROHet","ROHom"),"csv",sep=".")
  if(file.exists(report_filename)) file.remove(report_filename)

  # a vector to record if ROH if is written or not (FALSE, TRUE)
  is_run <- vector()

  # require "plyr"
  n_of_individuals <- daply(genotype,"IID",function(x) {

    gaps <- diff(mapFile$bps)
    y <- slidingWindow(as.integer(x[-c(1,2)]),gaps,windowSize,step=1,ROHet=ROHet,maxOppositeGenotype,maxMiss,maxGap);
    snpRun <- snpInRun(y,windowSize,threshold)
    dRUN <- createRUNdf(snpRun,mapFile,minSNP,minLengthBps,minDensity)

    # this function will write ROH on report_filename (defined inside writeRUN)
    is_run <- writeRUN(as.character(x$IID),dRUN,ROHet,as.character(x$FID))
    return(is_run)
  })

  return(sum(n_of_individuals))
}






