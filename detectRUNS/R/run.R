###########################################################
### Compute Genomic Runs in R (homozygosity/heterozygosity)
###########################################################


#' Main function to detect RUNS (ROHom/ROHet) using sliding windows (a la Plink)
#'
#' This is one of the main function of detectRUNS and is used to detect runs
#' (of homozygosity or heterozygosity)
#' in the genome (diploid) with the sliding-window method.
#' All parameters to detect runs (e.g. minimum n. of SNP, max n. of missing genotypes,
#' max n. of opposite genotypes etc.) are specified here.
#' Input data are in the ped/map
#' Plink format (https://www.cog-genomics.org/plink/1.9/input#ped)
#'
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param windowSize the size of sliding window (number of SNP loci) (default = 15)
#' @param threshold the threshold of overlapping windows of the same state
#' (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)
#' @param minSNP minimum n. of SNP in a RUN (default = 3)
#' @param ROHet should we look for ROHet or ROHom? (default = FALSE)
#' @param maxOppWindow max n. of homozygous/heterozygous SNP in the
#' sliding window (default = 1)
#' @param maxMissWindow max. n. of missing SNP in the sliding window (default = 1)
#' @param maxGap max distance between consecutive SNP to be still considered a
#' potential run (default = 10^6 bps)
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#' @param maxOppRun max n. of opposite genotype SNPs in the run (optional)
#' @param maxMissRun max n. of missing SNPs in the run (optional)
#'
#' @details
#' This function scans the genome (diploid) for runs using the sliding-window method.
#' This is a wrapper function for many component functions that handle the input data
#' (ped/map files), perform internal conversions, accept parameters specifications,
#' select whether runs of homozygosity (RoHom) or of heterozygosity (RoHet)
#' are looked for.
#'
#' In the ped file, the groups samples belong to can be specified (first column).
#' This is important if comparisons between human ethnic groups or between animal breeds
#' or plant varieties or biological populations are to be performed.
#' Also, if cases and controls are to be compared, this is the place where this
#' information needs to be specified.
#'
#' This function returns a data frame with all runs detected in the dataset.
#' This data frame can then be written out to a csv file.
#' The data frame is, in turn, the input for other functions of the detectRUNS package
#' that create plots and produce statistics from the results
#' (see plots and statistics functions in this manual,
#' and/or refer to the detectRUNS vignette).
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity in the analysed dataset.
#' The returned dataframe contains the following seven columns: "group", "id", "chrom",
#' "nSNP", "from", "to", "lengthBps" (group: population, breed, case/control etc.;
#' id: individual identifier; chrom: chromosome on which the run is located;
#' nSNP: number of SNPs in the run; from: starting position of the run, in bps;
#' to: end position of the run, in bps; lengthBps: size of the run)
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
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#' # calculating runs with sliding window approach
#' \dontrun{
#' # skipping runs calculation
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,
#' minSNP = 15, ROHet = FALSE,  maxOppWindow = 1, maxMissWindow = 1, maxGap=10^6,
#' minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'

slidingRUNS.run <- function(genotypeFile, mapFile, windowSize = 15, threshold = 0.05,
                            minSNP = 3, ROHet = FALSE, maxOppWindow = 1, maxMissWindow = 1,
                            maxGap = 10^6, minLengthBps = 1000, minDensity = 1/1000,
                            maxOppRun = NULL, maxMissRun = NULL) {

  message(paste("You are using the method: SLIDING WINDOWS"))

  # debug
  if (ROHet == TRUE) {
    message("Analysing Runs of Heterozygosity (ROHet)")
  } else if (ROHet == FALSE) {
    message("Analysing Runs of Homozygosity (ROHom)")
  } else {
    stop(paste("Unknown ROHet value:",ROHet, ". It MUST be only TRUE/FALSE (see documentation)"))
  }

  # if genotype is file, open file
  if(file.exists(genotypeFile)){
    conn  <- file(genotypeFile, open = "r")
  } else {
    stop(paste("file", genotypeFile, "doesn't exists"))
  }

  if(file.exists(mapFile)){
    # using data.table to read data
    mapFile <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }

  # setting colnames
  colnames(mapFile) <- c("Chrom","SNP","cM","bps")

  # collect all parameters in a variable
  parameters <- list(windowSize=windowSize, threshold=threshold, minSNP=minSNP,
                     ROHet=ROHet, maxOppWindow=maxOppWindow,
                     maxMissWindow=maxMissWindow, maxGap=maxGap, minLengthBps=minLengthBps,
                     minDensity=minDensity, maxOppRun=maxOppRun, maxMissRun=maxMissRun)

  # calculate gaps
  gaps <- diff(mapFile$bps)

  message(paste("Window size:", windowSize))
  message(paste("Threshold for calling SNP in a Run:", threshold))

  # initialize data.frame of results
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
    a_run <- slidingRuns(genotype, animal, mapFile, gaps, parameters)

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


#' Main function to detect genomic RUNS (ROHom/ROHet) using the consecutive method
#'
#' This is the main detectRUNS function to scan the genome for runs (of homozygosity or heterozygosity)
#' using the consecutive method (Marras et al. 2015, Animal Genetics 46(2):110-121).
#' All parameters to detect runs (e.g. minimum n. of SNP, max n. of missing genotypes,
#' max n. of opposite genotypes etc.) are specified here.
#' Input data are in the ped/map
#' Plink format (https://www.cog-genomics.org/plink/1.9/input#ped)
#'
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param ROHet should we look for ROHet or ROHom? (default = FALSE)
#' @param maxOppRun max n. of opposite genotype SNPs in the run (default = 0)
#' @param maxMissRun max n. of missing SNPs in the run (default = 0)
#' @param minSNP minimum n. of SNP in a RUN (default = 15)
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run (defaults to 10^6)
#'
#' @details
#' This function scans the genome (diploid) for runs using the consecutive method.
#' This is a wrapper function for many component functions that handle the input data (ped/map files), performs internal conversions,
#' accepts parameters specifications, selects the statistical method to detect runs (sliding windows, consecutive loci) and whether
#' runs of homozygosity (RoHom) or of heterozygosity (RoHet) are looked for.
#'
#' In the ped file, the groups samples belong to can be specified (first column). This is important if comparisons between
#' human ethnic groups or between animal breeds or plant varieties or biological populations are to be performed.
#' Also, if cases and controls are to be compared, this is the place where this information needs to be specified.
#'
#' This function returns a data frame with all runs detected in the dataset. This data frame can then be written out to a csv file.
#' The data frame is, in turn, the input for other functions of the detectRUNS package that create plots and produce statistics
#' of the results (see plot and statistic functions in this manual, and/or refer to the vignette of detectRUNS).
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity in the analysed dataset.
#' The returned dataframe contains the following seven columns: "group", "id", "chrom",
#' "nSNP", "from", "to", "lengthBps" (group: population, breed, case/control etc.;
#' id: individual identifier; chrom: chromosome on which the run is located;
#' nSNP: number of SNPs in the run; from: starting position of the run, in bps;
#' to: end position of the run, in bps; lengthBps: size of the run)
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
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#' # calculating runs with consecutive run approach
#' \dontrun{
#' # skipping runs calculation
#' runs <- consecutiveRUNS.run(genotypeFile, mapFile, minSNP = 15, ROHet = FALSE,
#' maxOppRun = 0, maxMissRun = 0, maxGap=10^6,
#' minLengthBps = 100000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.consecutive.csv", package="detectRUNS")
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'

consecutiveRUNS.run <- function(genotypeFile, mapFile, ROHet = FALSE,
                                maxOppRun = 0, maxMissRun = 0,
                                minSNP = 15, minLengthBps = 1000,
                                maxGap = 10^6) {

  message(paste("You are using the method: CONSECUTIVE RUNS"))

  # debug
  if (ROHet == TRUE) {
    message("Analysing Runs of Heterozygosity (ROHet)")
  } else if (ROHet == FALSE) {
    message("Analysing Runs of Homozygosity (ROHom)")
  } else {
    stop(paste("Unknown ROHet value:",ROHet, ". It MUST be only TRUE/FALSE (see documentation)"))
  }

  # if genotype is file, open file
  if(file.exists(genotypeFile)){
    conn  <- file(genotypeFile, open = "r")
  } else {
    stop(paste("file", genotypeFile, "doesn't exists"))
  }

  if(file.exists(mapFile)){
    # using data.table to read data
    mapFile <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }

  # setting colnames
  colnames(mapFile) <- c("Chrom","SNP","cM","bps")

  # initialize data.frame of results
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

    a_run <- consecutiveRunsCpp(genotype, animal, mapFile=mapFile, ROHet=ROHet, minSNP=minSNP,
                                maxOppositeGenotype=maxOppRun, maxMiss=maxMissRun,
                                minLengthBps=minLengthBps, maxGap = maxGap)

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
