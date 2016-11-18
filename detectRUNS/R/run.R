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
#' @return a dataframe with RUNs of Homozygosity or Heterozygosity
#' @export
#'
#' @import plyr
#' @import itertools
#' @import ggplot2
#' @import itertools
#' @import data.table
#' @import utils
#' @importFrom bigmemory read.big.matrix
#'
#' @examples
#' data(chillingam)
#' x <- RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000, minDensity = 1/10)
#'

# TODO add output file parameter
RUNS.run <- function(genotype, mapFile, windowSize = 15, threshold = 0.1, minSNP = 3, ROHet = TRUE,
                     maxOppositeGenotype = 1, maxMiss = 1, maxGap = 10^6, minLengthBps = 1000, minDensity = 1/10) {

  # TODO: check for output_file existance

  # debug
  if (ROHet == TRUE) {
    message("Analysing Runs of Heterozygosity (ROHet)")
  } else if (ROHet == FALSE) {
    message("Analysing Runs of Homozygosity (ROHom)")
  } else {
    stop(paste("Unknown ROHet value:",ROHet, ". It MUST be only TRUE/FALSE (see documentation)"))
  }

  message(paste("Window size:", windowSize))
  message(paste("Threshold for calling SNP in a Run:", threshold))

  if(!is.data.frame(genotype)) {

    if(file.exists(genotype)){
      genotype_path <- genotype
      # using bigmemory to read data
      genotype <- read.big.matrix(genotype_path, sep = " ", header = T, type = "short")
      colClasses <- c(
        rep("character", 2),
        rep("NULL", ncol(genotype)-2)
      )

      animals <- fread(genotype_path, sep = " ", header = T, colClasses = colClasses)
    } else {
      stop(paste("file", genotype, "doesn't exists"))
    }

  } else {
    # read animals properly
    animals <- genotype[ ,c(1,2)]
  }

  if(!is.data.frame(mapFile)) {

    if(file.exists(mapFile)){
      # using data.table to read data
      mapFile <- fread(mapFile, header = F)
    } else {
      stop(paste("file", mapFile, "doesn't exists"))
    }
  }

  # TODO: check that genotype columns and mapFile rows (+6) are identical

  # setting colnames
  names(mapFile) <- c("Chrom","SNP","cM","bps")

  #remove unnecessary fields from the .raw file
  genotype <- genotype[ ,-c(1:6)]

  # record all runs in a dataframe
  RUNs <- data.frame(breed=character(), id=character(), chrom=character(), nSNP=integer(),
                     from=integer(), to=integer(), lengthBps=integer())

  # calculate gaps
  gaps <- diff(mapFile$bps)

  # define an internal function to call other functions
  find_run <- function(x, animal) {
    # get individual and breed
    ind <- as.character(animal$IID)
    breed <- as.character(animal$FID)

    # use sliding windows
    y <- slidingWindowCpp(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);
    snpRun <- snpInRunCpp(y,windowSize,threshold)
    dRUN <- createRUNdf(snpRun,mapFile,minSNP,minLengthBps,minDensity)

    # manipulate dRUN to order columns
    dRUN$id <- rep(ind, nrow(dRUN))
    dRUN$breed <- rep(breed, nrow(dRUN))
    dRUN <- dRUN[,c(7,6,4,3,1,2,5)]

    # debug
    if(nrow(dRUN) > 0) {
      message(paste("N. of RUNS for individual", ind, "is:", nrow(dRUN)))
    } else {
      message(paste("No RUNs found for animal",ind))
    }

    #return RUNs to caller
    return(dRUN)
  }

  for (i in 1:nrow(genotype)) {
    a_run <- find_run(
      as.integer(genotype[i, ]),
      animal = animals[i, ]
    )

    # bind this run (if has rows) to others RUNs (if any)
    RUNs <- rbind(RUNs, a_run)
  }

  # fix row names
  row.names(RUNs) <- NULL

  # TODO: write output file if requested

  # return calculated runs (data.frame)
  return(RUNs)
}

