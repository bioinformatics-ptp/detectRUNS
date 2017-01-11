##################################
### Compute ROH in R (a la Plink)
##################################


#' Run functions to detect RUNS (ROHom/ROHet)
#'
#' This is the main function of RUNS and would probably require
#' some good documentation.
#'
#' @param genotype_path genotype (.ped) file location
#' @param mapfile_path map file (.map) file location
#' @param windowSize the size of sliding window
#' @param threshold the threshold of overlapping windows of the same state (homozygous/heterozygous) to call a SNP in a RUN
#' @param minSNP minimum n. of SNP in a RUN
#' @param ROHet should we look for ROHet or ROHom?
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#' @param method one of either 'slidingWindow' (a la Plink) or 'consecutiveRuns' (a la Marras)
#'
#' @return a dataframe with RUNs of Homozygosity or Heterozygosity
#' @export
#'
#' @import plyr
#' @import itertools
#' @import ggplot2
#' @import itertools
#' @import utils
#'
#' @examples
#' # getting map and ped paths
#' genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
#' mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
#' runs <- RUNS.run(genotype_path, mapfile_path, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  maxGap=10^6, minLengthBps = 1000,
#' minDensity = 1/10, method='slidingWindow')
#'

# TODO add output file parameter
RUNS.run <- function(genotype_path, mapfile_path, windowSize = 15, threshold = 0.1,
                     minSNP = 3, ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,
                     maxGap = 10^6, minLengthBps = 1000, minDensity = 1/10,
                     method=c('slidingWindow','consecutiveRuns')) {

  # debug
  if (ROHet == TRUE) {
    message("Analysing Runs of Heterozygosity (ROHet)")
  } else if (ROHet == FALSE) {
    message("Analysing Runs of Homozygosity (ROHom)")
  } else {
    stop(paste("Unknown ROHet value:",ROHet, ". It MUST be only TRUE/FALSE (see documentation)"))
  }

  # if genotype is file, open file
  if(file.exists(genotype_path)){
    conn  <- file(genotype_path, open = "r")
  } else {
    stop(paste("file", genotype_path, "doesn't exists"))
  }

  if(file.exists(mapfile_path)){
    # using data.table to read data
    mapFile <- data.table::fread(mapfile_path, header = F)
  } else {
    stop(paste("file", mapfile_path, "doesn't exists"))
  }

  # setting colnames
  colnames(mapFile) <- c("Chrom","SNP","cM","bps")

  # define an internal function to call other functions
  find_run <- function(x, animal) {
    # get individual and group
    ind <- as.character(animal$IID)
    group <- as.character(animal$FID)

    # use sliding windows
    y <- slidingWindowCpp(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);
    snpRun <- snpInRunCpp(y,windowSize,threshold)
    dRUN <- createRUNdf(snpRun,mapFile,minSNP,minLengthBps,minDensity)

    # manipulate dRUN to order columns
    dRUN$id <- rep(ind, nrow(dRUN))
    dRUN$group <- rep(group, nrow(dRUN))
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

  method <- match.arg(method)

  message(paste("You are using the method:", method))

  if(method=='slidingWindow') {

    # calculate gaps
    gaps <- diff(mapFile$bps)

    message(paste("Window size:", windowSize))
    message(paste("Threshold for calling SNP in a Run:", threshold))
  }


  RUNs <- data.frame(group=character(), id=character(), chrom=character(), nSNP=integer(),
                     from=integer(), to=integer(), lengthBps=integer())

  # read file line by line (http://stackoverflow.com/questions/4106764/what-is-a-good-way-to-read-line-by-line-in-r)
  while (length(oneLine <- readLines(conn, n = 1, warn = FALSE)) > 0) {
    genotype <- (strsplit(oneLine, " "))
    genotype <- as.character(genotype[[1]])

    # check that genotype columns and mapFile rows (+6) are identical
    if (length(genotype)-6 != nrow(mapFile)*2) {
      stop("Number of markers differ in mapFile and genotype: are those file the same dataset?")
    }

    # get animal
    animal <- list(FID=genotype[1], IID=genotype[2])

    # convert into genotype (use from 7th column to last column)
    genotype <- pedConvertCpp(genotype[7:length(genotype)])

    # find run for this genotype
    if(method=='slidingWindow') {

      a_run <- find_run(genotype, animal)
    } else {

      a_run <- consecutiveRuns(genotype, animal, mapFile=mapFile, ROHet=ROHet, minSNP=minSNP,
                               maxOppositeGenotype=maxOppositeGenotype,
                               maxMiss=maxMiss, maxGap = maxGap)
    }

    # bind this run (if has rows) to others RUNs (if any)
    RUNs <- rbind(RUNs, a_run)

  }

  # close input stream
  close(conn)

  # fix row names
  row.names(RUNs) <- NULL

  # return calculated runs (data.frame)
  return(RUNs)
}

