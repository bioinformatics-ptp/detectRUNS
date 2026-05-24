#####################
## FUNCTIONS FOR RUNS
#####################


#' Convert 0/1/2 genotypes to 0/1
#'
#' This is a utility function, that convert 0/1/2 genotypes (AA/AB/BB) into 0/1
#' (either homozygous/heterozygous)
#'
#' @param x vector of 0/1/2 genotypes
#' @keywords internal
#' @return converted vector of genotypes (0/1)
#'

genoConvert <- function(x) {
  new <- c(0,1,0,NA)
  old <- c(0,1,2,NA)
  return(new[match(x,old)])
}


#' Read from a .map file locations and return a data.table object
#' 
#' This is an utility function which check for file existance, define 
#' colClasses and then returns the read data.table object
#' @param mapFile map file (.map) file path
#' @keywords internal
#' @return data.table object
#' 

readMapFile <- function(mapFile) {
  # define colClasses
  colClasses <- c("character", "character", "character", "numeric")
  
  if(file.exists(mapFile)){
    # using data.table to read data
    mappa <- data.table::fread(mapFile, header = FALSE, colClasses = colClasses)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }
  
  # set column names
  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL
  
  return(mappa)
}


#' Function to check whether a window is (loosely) homozygous or not
#'
#' This is a core function. Parameters on how to consider a window homozygous are here (maxHet, maxMiss)
#'
#' @param x vector of 0/1 genotypes (from genoConvert())
#' @param gaps vector of differences between consecutive positions (gaps) in bps
#' @param maxHet max n. of heterozygous SNP in a homozygous window
#' @param maxMiss max n. of missing in a window
#' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
#' @param i index along the genome (genome-vector for each individual)
#' @param windowSize size of window (n. of SNP)
#'
#' @return a list: i) TRUE/FALSE (whether a window is heterozygous or NOT); ii) indexes of "opposite and missing" genotype
#'

homoZygotTest <- function(x,gaps,maxHet,maxMiss,maxGap,i,windowSize) {

  nHet <- sum(x==1,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  oppositeAndMissingSNP <- c(-1,0,9)[match(x,c(0,1,NA))]
  oppositeAndMissingSNP <- oppositeAndMissingSNP[oppositeAndMissingSNP!=-1]
  indexSNP <- seq(i,i+windowSize-1)[which(x==1 | is.na(x))]
  names(oppositeAndMissingSNP) <- indexSNP

  windowStatus <- ifelse(!(nHet > maxHet | nMiss > maxMiss), TRUE, FALSE)
  return(list("windowStatus"=windowStatus,"oppositeAndMissingSNP"=oppositeAndMissingSNP))
}


#' Function to check whether a window is (loosely) heterozygous or not
#'
#' This is a core function within the sliding-window workflow. Parameters on how to consider a window heterozygous are here (maxHom, maxMiss)
#'
#' @param x vector of 0/1 genotypes (from genoConvert())
#' @param gaps vector of differences between consecutive positions (gaps) in bps
#' @param maxHom max n. of homozygous SNP in a heterozygous window
#' @param maxMiss max n. of missing in a window
#' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
#' @param i index along the genome (genome-vector for each individual)
#' @param windowSize size of window (n. of SNP)
#'
#' @return a list: i) TRUE/FALSE (whether a window is heterozygous or NOT); ii) indexes of "opposite and missing" genotype
#'

heteroZygotTest <- function(x,gaps,maxHom,maxMiss,maxGap,i,windowSize) {

  nHom <- sum(x==0,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  oppositeAndMissingSNP <- c(0,-1,9)[match(x,c(0,1,NA))]
  oppositeAndMissingSNP <- oppositeAndMissingSNP[oppositeAndMissingSNP!=-1]
  indexSNP <- seq(i,i+windowSize-1)[which(x==0 | is.na(x))]
  names(oppositeAndMissingSNP) <- indexSNP

  windowStatus <- ifelse(!(nHom > maxHom | nMiss > maxMiss), TRUE, FALSE)
  return(list("windowStatus"=windowStatus,"oppositeAndMissingSNP"=oppositeAndMissingSNP))
}


#' Function to slide a window over a vector (individual's genotypes)
#'
#' This is a core function. The functions to detect RUNS are slid over the genome
#'
#' @param data vector of 0/1/2 genotypes
#' @param gaps vector of differences between consecutive positions (gaps) in bps
#' @param windowSize size of window (n. of SNP)
#' @param step by which (how many SNP) is the window slid
#' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
#' @param ROHet shall we detect ROHet or ROHom?
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#'
#' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
#'

slidingWindow <- function(data, gaps, windowSize, step, maxGap, ROHet=TRUE, maxOppositeGenotype=1, maxMiss=1) {

  data_length <- length(data)
  spots <- seq(from = 1, to = (data_length - windowSize + 1), by = step)
  result <- vector(length = length(spots))
  oppositeAndMissingGenotypes <- array(character(0))
  y <- genoConvert(data)

  message(paste("Analysing",ifelse(ROHet,"Runs of Heterozygosity (ROHet)","Runs of Homozygosity (ROHom)"),sep=" "))

  if(ROHet) {

    for(i in 1:length(spots)){
      ret <- heteroZygotTest(y[spots[i]:(spots[i]+windowSize-1)],gaps[spots[i]:(spots[i]+windowSize-2)],maxOppositeGenotype,maxMiss,maxGap,i,windowSize)
      result[i] <- ret$windowStatus
      # this will append the returned value of heteroZygotTest to existsting
      # oppositeAndMissingGenotypes array. Since windows slide with sovraposition,
      # we append new values to oppositeAndMissingGenotypes. We may calculate
      # this array once.
      oppositeAndMissingGenotypes <- c(oppositeAndMissingGenotypes,ret$oppositeAndMissingSNP[!(names(ret$oppositeAndMissingSNP) %in% names(oppositeAndMissingGenotypes))])
    }
    # to include a shrinking sliding-window at the end of the chromosome/genome, uncomment the following line
    # for(i in (length(spots)+1):data_length) result[i] <- heteroZygotTest(y[seq(i,data_length)],maxOppositeGenotype,maxMiss)
  } else {

    for(i in 1:length(spots)){
      ret <- homoZygotTest(y[spots[i]:(spots[i]+windowSize-1)],gaps[spots[i]:(spots[i]+windowSize-2)],maxOppositeGenotype,maxMiss,maxGap,i,windowSize)
      result[i] <- ret$windowStatus
      oppositeAndMissingGenotypes <- c(oppositeAndMissingGenotypes,ret$oppositeAndMissingSNP[!(names(ret$oppositeAndMissingSNP) %in% names(oppositeAndMissingGenotypes))])
    }
    # to include a shrinking sliding-window at the end of the chromosome/genome, uncomment the following line
    #for(i in (length(spots)+1):data_length) result[i] <- homoZygotTest(y[seq(i,data_length)],maxOppositeGenotype,maxMiss)
  }

  # print(paste(
  #   "Length of homozygous windows overlapping SNP loci (should be equal to the n. of SNP in the file):",
  #   length(result),sep=" "))

  return(list("windowStatus"=result,"oppositeAndMissingGenotypes"=oppositeAndMissingGenotypes))

}


#' Function to return a vector of T/F for whether a SNP is or not in a RUN
#'
#' This is a core function. The function to determine whether a SNP is or not in a RUN.
#' The ratio between homozygous/heterozygous windows and total n. of windows is computed here
#'
#' @param RunVector vector of TRUE/FALSE (is a window homozygous/heterozygous?)
#' @param windowSize size of window (n. of SNP)
#' @param threshold threshold to call a SNP in a RUN
#'
#' @return vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
#'

snpInRun <- function(RunVector,windowSize,threshold) {

  RunVector_length <- length(RunVector)

  # print(paste("Length of input vector:",RunVector_length,sep=" "))
  # print(paste("Window size:",windowSize,sep=" "))
  # print(paste("Threshold for calling SNP in a Run:",threshold,sep=" "))

  # compute total n. of overlapping windows at each SNP locus (see Bjelland et al. 2013)
  nWin <- c(seq(1,windowSize), rep(windowSize,(RunVector_length-windowSize-1)), seq(windowSize,1))

  # compute n. of homozygous/heterozygous windows that overlap at each SNP locus (Bjelland et al. 2013)
  # create two sets of indices to slice the vector of windows containing or not a run (RunVector)
  ind1 <- c(rep(1L, windowSize - 1L), seq_len(RunVector_length))
  ind2 <- c(seq_len(RunVector_length), rep(RunVector_length, windowSize - 1L))
  hWin <- mapply(function(a, b) sum(RunVector[a:b]), ind1, ind2)

  # ratio between homozygous/heterozygous windows and total overlapping windows at each SNP
  quotient <- hWin/nWin


  #vector of SNP belonging to a ROH (>= matches PLINK --homozyg behaviour)
  snpRun <- ifelse(quotient >= threshold, TRUE, FALSE)
  # print(paste(
  #   "Lenght of output file:",
  #   length(snpRun),sep=" "))

  return(snpRun)
}


#' Function to create a dataframe of RUNS per individual animal
#' Requires a map file (other filename to read or R object)
#' Parameters on maximum number of missing and opposite genotypes in the run (not the window) are implemented here
#'
#'
#' @param snpRun vector of TRUE/FALSE (is the SNP in a RUN?)
#' @param mapFile Plink-like map file (data.frame)
#' @param minSNP minimum n. of SNP to call a RUN
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#' @param oppositeAndMissingSNP indexed array of missing and opposite genotypes (SNP order in the genome is the index)
#' @param maxOppRun max n. of opposite genotype SNPs in the run (not in the window!)
#' @param maxMissRun max n. of missing SNPs in the run (not in the window!)
#'
#' @return a data.frame with RUNS per animal
#'
#' @import utils
#' @importFrom stats na.omit
#'

createRUNdf <- function(snpRun, mapFile, minSNP = 3, minLengthBps = 1000,
                        minDensity = 1/10, oppositeAndMissingSNP, maxOppRun=NULL,
                        maxMissRun=NULL, maxGap=NULL) {

  bps_all <- mapFile$bps
  dd <- cbind.data.frame(snpRun,"Chrom"=mapFile$Chrom,"n"=seq(1,nrow(mapFile)))

  parts <- split(dd, dd$Chrom)
  dL <- do.call(rbind, lapply(parts, function(x) {

    # define where RUNs change states
    # cutPoints for "from" and "to" on the original snpRun vector
    cutPoints <- x[which(diff(sign(x$snpRun)) != 0),"n"]
    from <- c(x$n[1], cutPoints + 1)
    to <- c(cutPoints, x$n[length(x$snpRun)])
    # internal cutPoints for n. SNP calculations
    cutPoints <- which(diff(sign(x$snpRun)) != 0)
    from_bis <- c(1, cutPoints + 1)
    to_bis <- c(cutPoints, length(x$snpRun))
    # count SNPs in each run segment
    lengte <- mapply(function(a, b) sum(x$snpRun[a:b]), from_bis, to_bis)

    # Run-level gap split (PLINK-compatible): for each TRUE run segment,
    # split into sub-runs wherever two consecutive SNPs have gap > maxGap.
    # Both boundary SNPs are kept (no SNP loss, unlike the snpRun-vector trick).
    if (!is.null(maxGap) && maxGap > 0L && any(lengte > 0L)) {
      exp_from <- integer(0); exp_to <- integer(0); exp_nsnp <- integer(0)
      for (i in seq_along(from)) {
        f <- from[i]; t <- to[i]
        if (lengte[i] == 0L) {
          exp_from <- c(exp_from, f); exp_to <- c(exp_to, t); exp_nsnp <- c(exp_nsnp, 0L)
        } else {
          cur_f <- f
          for (j in seq(f, t - 1L)) {
            if (bps_all[j + 1L] - bps_all[j] > maxGap) {
              exp_from <- c(exp_from, cur_f); exp_to <- c(exp_to, j)
              exp_nsnp <- c(exp_nsnp, j - cur_f + 1L)
              cur_f <- j + 1L
            }
          }
          exp_from <- c(exp_from, cur_f); exp_to <- c(exp_to, t)
          exp_nsnp <- c(exp_nsnp, t - cur_f + 1L)
        }
      }
      from <- exp_from; to <- exp_to; lengte <- exp_nsnp
    }

    # get n of rows
    n_rows <- length(lengte)

    data.frame("from"=from,
               "to"=to,
               "nSNP"=lengte,
               "chrom"=character(n_rows),
               "lengthBps"=numeric(n_rows), stringsAsFactors=FALSE)
  }))

  # filter RUNs by minSNP
  dL <- dL[dL$nSNP>=minSNP, ]
  dL <- na.omit(dL)

  # return if all rows are filtered
  if (nrow(dL) == 0) {
    return(dL)
  }

  chroms <- mapFile[dL$from, "Chrom"]

  # debug
  # print(chroms)

  dL$from <- mapFile[dL$from, "bps"]
  dL$to <- mapFile[dL$to,"bps"]

  # setting other values
  dL$chrom <- as.character(chroms)
  dL$lengthBps <- (dL$to - dL$from + 1L)

  # filters on minimum run length and minimum SNP density
  dL <- dL[dL$lengthBps >= minLengthBps,]
  dL$SNPdensity <- (dL$nSNP/dL$lengthBps)*1000 # n. SNP per kbps
  dL <- dL[dL$SNPdensity >= minDensity, ]
  dL$SNPdensity <- NULL

  # return if all rows are filtered
  if (nrow(dL) == 0) {
    return(dL)
  }

  # filters on max heterozygotes and missing in a run
  if(!is.null(maxOppRun) || !is.null(maxMissRun)) {
    # Add map information to opposite and missing SNPs
    W <- cbind.data.frame(oppositeAndMissingSNP)
    W <- cbind.data.frame(W, mapFile[as.numeric(row.names(W)), ])

    # Add nOpp and nMiss columns to dataframe
    dL$nOpp  <- mapply(function(chr, f, t) sum(W$Chrom==chr & W$bps>=f & W$bps<=t & W$oppositeAndMissingSNP==0), dL$chrom, dL$from, dL$to)
    dL$nMiss <- mapply(function(chr, f, t) sum(W$Chrom==chr & W$bps>=f & W$bps<=t & W$oppositeAndMissingSNP==9), dL$chrom, dL$from, dL$to)

    if(!is.null(maxOppRun)) {
      # filter RUNs by opposite SNPs
      dL <- dL[dL$nOpp <= maxOppRun,]
    }

    if (!is.null(maxMissRun)) {
      # filter RUNs by missing SNPs
      dL <- dL[dL$nMiss <= maxMissRun,]
    }

    # remove nOpp and nMiss columns
    dL <- dL[,-c(6,7)]

  }

  # fix row.names: whitout it, dataframe will have row names like unfiltered one
  row.names(dL) <- NULL

  return(dL)
}


#' Function to write out RUNS per individual animal
#'
#'
#' @param ind ID of animals
#' @param dRUN data.frame with RUNS per animal
#' @param ROHet shall we detect ROHet or ROHom?
#' @param group group (factor): population, breed, ethnicity, case/control etc.
#' @param outputName output filename
#'
#' @return TRUE/FALSE if RUNS are written out or not
#'
#' @import utils
#'

writeRUN <- function(ind, dRUN, ROHet=TRUE, group, outputName) {

  dRUN$id <- rep(ind,nrow(dRUN))
  dRUN$group <- rep(group,nrow(dRUN))
  dRUN <- dRUN[,c(7,6,4,3,1,2,5)]

  if(nrow(dRUN) > 0) {
    append = FALSE
    headers = TRUE

    if(file.exists(outputName)){
      warning(paste("Appending to", outputName, "..."))
      append = TRUE
      headers = FALSE
    }
    utils::write.table(
      sep=';', dRUN, file=outputName, quote=FALSE,
      col.names=headers, row.names=FALSE, append=append
    )
    is_run <- TRUE
  } else {
    message(paste("No RUNs found for sample",ind,sep=" "))
    is_run <- FALSE
  }
  return(is_run)
}


#' Function to count number of times a SNP is in a RUN. Need to be called per
#' chromosome (not using a dataframe with all results on all chromosomes)
#'
#'
#' @param runsChrom R object (dataframe) with results per chromosome (column names:"POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")
#' @param mapChrom R object (dataframe) with SNP name and position per chromosome (map file) (column names: "CHR","SNP_NAME","POSITION")
#' @param sample_info data.frame with columns group and id (one row per individual)
#'
#' @return dataframe with counts per SNP in runs (per population)
#' @keywords internal
#'
#' @import utils
#'

snpInsideRuns <- function(runsChrom, mapChrom, sample_info) {

  unique_groups <- sort(unique(runsChrom$POPULATION))
  k             <- nrow(mapChrom)
  snp_pos       <- mapChrom$POSITION   # must be sorted ascending

  results <- vector("list", length(unique_groups))

  for (i in seq_along(unique_groups)) {
    grp    <- unique_groups[i]
    runs_g <- runsChrom[runsChrom$POPULATION == grp, ]
    nGroup <- sum(sample_info$group == as.character(grp))

    if (nrow(runs_g) == 0L || k == 0L) {
      count_g <- integer(k)
    } else {
      # Convert run bp coordinates to SNP indices (1-based).
      # from_v: index of first SNP whose position >= START
      # to_v  : index of last  SNP whose position <= END
      from_v <- findInterval(runs_g$START - 1L, snp_pos) + 1L
      to_v   <- findInterval(runs_g$END,         snp_pos)

      valid <- from_v <= to_v & from_v >= 1L & to_v <= k
      if (!any(valid)) {
        count_g <- integer(k)
      } else {
        fv <- from_v[valid]
        tv <- to_v[valid]
        # Sweep-line: delta[j]+=1 at run start, delta[j]-=1 after run end.
        # cumsum(delta)[1:k] gives the number of runs covering each SNP.
        delta   <- tabulate(fv, nbins = k + 1L) -
                   tabulate(tv + 1L, nbins = k + 1L)
        count_g <- cumsum(delta)[seq_len(k)]
      }
    }

    results[[i]] <- data.frame(
      SNP_NAME   = mapChrom$SNP_NAME,
      CHR        = mapChrom$CHR,
      POSITION   = snp_pos,
      COUNT      = count_g,
      BREED      = as.factor(rep(grp, k)),
      PERCENTAGE = (count_g / nGroup) * 100,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Function to detect runs using sliding window approach
#'
#' This is a core function not intended to be exported
#'
#' @param indGeno vector of 0/1/NAs of individual genotypes (0: homozygote; 1: heterozygote)
#' @param individual list of group (breed, population, case/control etc.) and ID of individual sample
#' @param mapFile Plink map file (for SNP position)
#' @param gaps distance between SNPs
#' @param parameters list of parameters
#' @param cpp use cpp functions or not (DEBUG)
#'
#' @details
#' This method uses sliding windows to detect RUNs. Checks on minimum n. of SNP, max n. of opposite and missing genotypes,
#' max gap between adjacent loci and minimum length of the run are implemented (as in the sliding window method).
#' Both runs of homozygosity (RoHom) and of heterozygosity (RoHet) can be search for (option ROHet: TRUE/FALSE)
#' NOTE: this methods is intended to not be exported
#'
#' @return A data frame of runs per individual sample
#' @keywords internal
#'

slidingRuns <- function(indGeno, individual, mapFile, gaps, parameters, cpp=TRUE) {
  # get individual and group
  ind <- as.character(individual$IID)
  group <- as.character(individual$FID)

  # use sliding windows (check cpp)
  if (cpp == TRUE) {
    res <- slidingWindowCpp(indGeno, gaps, parameters$windowSize, step=1,
                            parameters$maxGap, parameters$ROHet, parameters$maxOppWindow,
                            parameters$maxMissWindow);

    snpRun <- snpInRunCpp(res$windowStatus, parameters$windowSize, parameters$threshold)
  } else {
    res <- slidingWindow(indGeno, gaps, parameters$windowSize, step=1,
                            parameters$maxGap, parameters$ROHet, parameters$maxOppWindow,
                            parameters$maxMissWindow);

    snpRun <- snpInRun(res$windowStatus, parameters$windowSize, parameters$threshold)
  }

  # TODO: check arguments names
  dRUN <- createRUNdf(snpRun, mapFile, parameters$minSNP, parameters$minLengthBps,
                      parameters$minDensity, res$oppositeAndMissingGenotypes,
                      parameters$maxOppRun, parameters$maxMissRun,
                      maxGap = parameters$maxGap)

  # manipulate dRUN to order columns
  dRUN <- as.data.frame(dRUN)
  dRUN$id <- rep(ind, nrow(dRUN))
  dRUN$group <- rep(group, nrow(dRUN))
  dRUN <- dRUN[,c(7,6,4,3,1,2,5)]

  #return RUNs to caller
  return(dRUN)
}


#' Function to detect consecutive runs in a vector (individual's genotypes)
#'
#' This is a core function. It implements the consecutive method for detection of runs in diploid genomes
#' (see Marras et al. 2015)
#'
#' @param indGeno vector of 0/1/NAs of individual genotypes (0: homozygote; 1: heterozygote)
#' @param individual list of group (breed, population, case/control etc.) and ID of individual sample
#' @param mapFile Plink map file (for SNP position)
#' @param ROHet shall we detect ROHet or ROHom?
#' @param minSNP minimum number of SNP in a run
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#' @param minLengthBps min length of a run in bps
#' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
#'
#' @details
#' The consecutive method detect runs by consecutively scanning SNP loci along the genome.
#' No sliding windows are used. Checks on minimum n. of SNP, max n. of opposite and missing genotypes,
#' max gap between adjacent loci and minimum length of the run are implemented (as in the sliding window method).
#' Both runs of homozygosity (RoHom) and of heterozygosity (RoHet) can be search for (option ROHet: TRUE/FALSE)
#'
#' @return A data frame of runs per individual sample
#'
#' @keywords internal
#'

consecutiveRuns <- function(indGeno, individual, mapFile, ROHet=TRUE, minSNP=3,
                            maxOppositeGenotype=1, maxMiss=1, minLengthBps=1000,
                            maxGap=10^6) {

  # runs of heterozygosity or of homozygosity?
  typ <- ifelse(ROHet,1,0)

  # animal data like IID or FID
  ind <- as.character(individual$IID)
  group <- as.character(individual$FID)

  # initialize variables. First chromosome in ordered mapFile
  lastChrom <- mapFile$Chrom[1]
  lastPos <- mapFile$bps[1]

  # a function to initialize runData from a position
  initializeRun <- function(chrom, start) {
    return(list("nOpposite"=0, "nMiss"=0, "runH"=0, "lengte"=0,
                "start"=start, chrom=chrom, end=NULL))
  }

  # a function to update RUNs data
  updateRUNS <- function(res, runData) {
    new_res <- data.frame("group"=group, "id"=ind, "chrom"=as.character(runData$chrom),
                          "nSNP"=runData$runH, "from"=runData$start,"to"=runData$end,
                          "lengthBps"=runData$lengte, stringsAsFactors = FALSE)
    return(rbind(res, new_res))
  }


  # a flag to determine if I'm in a RUN or not
  flag_run <- FALSE
  runData <- NULL

  # initialize dataframe of results. Defining data types accordingly slinding window
  res <- data.frame("group"=character(0),"id"=character(0),"chrom"=character(0),"nSNP"=integer(0),
                    "from"=integer(0),"to"=integer(0),"lengthBps"=integer(0), stringsAsFactors = FALSE)

  ##########################################################################################
  for (i in seq_along(indGeno)) {
    # Check for Chromosome
    currentChrom <- mapFile$Chrom[i]

    # get current position
    currentPos <- mapFile$bps[i]

    # test if chromosome is changed
    if (currentChrom != lastChrom ) {
      # if I have run, write to file
      # message(paste("Chromosome change", currentChrom, lastChrom))

      if(flag_run == TRUE && runData$runH >= minSNP && runData$lengte >= minLengthBps) {
        # message("Update RUN: chromosome changed")
        res <- updateRUNS(res, runData)
      }

      # update chrom and positions
      lastChrom <- currentChrom
      lastPos <- currentPos

      # unset RUN flag. New runs with new chromosomes!!!
      flag_run <- FALSE
      runData <- NULL
    }

    # calculate gap between consecutive SNP (in the same chrom)
    gap <- (currentPos - lastPos)

    # check if current gap is larger than max allowed gap. No matter current SNP
    if (flag_run == TRUE && gap >= maxGap) {
      # message("Gap condition")

      if(runData$runH >= minSNP && runData$lengte >= minLengthBps) {
        res <- updateRUNS(res, runData)
      }

      # unset RUN flag
      flag_run <- FALSE
      runData <- NULL
    }

    # Start for ==. Is a new RUN or not? ensure that indGeno[i] is a number
    if (indGeno[i] == typ && !is.na(indGeno[i])){
      # initialize run if not yet initialized, or just written after a big GAP
      if (flag_run == FALSE) {
        # message("Run initialized")
        # runData is a list of attributes for the current RUN
        runData <- initializeRun(currentChrom, currentPos)
        flag_run <- TRUE
      }

      # update runData values
      runData$runH <- runData$runH+1
      runData$end <- currentPos
      runData$lengte <- (runData$end - runData$start + 1L)

    } # condition: the genotype I want

    # start if !=
    else if (indGeno[i] != typ && !is.na(indGeno[i])){
      # if not in a run, don't do anything
      if (flag_run == FALSE) {
        next
      }

      # update nOpposite genotypes (in a run)
      runData$nOpposite <- runData$nOpposite + 1

      # check if maxOppositeGenotype is reached
      if (runData$nOpposite <= maxOppositeGenotype){
        # update runData values. This opposite genotype is a part of the RUN
        runData$runH <- runData$runH+1
        runData$end <- currentPos
        runData$lengte <- (runData$end - runData$start + 1L)

      } else {
        # message("max opposite reached")
        if (runData$runH >= minSNP && runData$lengte >= minLengthBps) {
          res <- updateRUNS(res, runData)
        }

        # unset RUN flag
        flag_run <- FALSE
        runData <- NULL

      } # condition nOpposite greather than maxOppositeGenotype

    } # condition: opposite genotype

    # start if 'NA'
    else if (is.na(indGeno[i])) {
      # if not in a run, don't do anything
      if (flag_run == FALSE) {
        next
      }

      # update missing values
      runData$nMiss <- runData$nMiss + 1

      # check if maxMiss is reached
      if (runData$nMiss <= maxMiss){
        # update runData values. This missing genotype is a part of the RUN
        runData$runH <- runData$runH+1
        runData$end <- currentPos
        runData$lengte <- (runData$end - runData$start + 1L)

      }
      else {
        # message("max missing reached")
        if (runData$runH >= minSNP && runData$lengte >= minLengthBps) {
          res <- updateRUNS(res, runData)
        }

        # unset RUN flag
        flag_run <- FALSE
        runData <- NULL

      } # condition: nMissing greather than permitted

    } # condition: missing genotype

    # update positions
    lastPos <- currentPos

  } # cicle: for i in genotypes

  # last snp if it is in a run
  if (flag_run == TRUE ) {
    # message("last SNP")
    if (runData$runH >= minSNP && runData$lengte >= minLengthBps) {
      res <- updateRUNS(res, runData)
    }

    # unset RUN flag
    flag_run <- FALSE
    runData <- NULL
  }

  return(res)
}

#' Read runs from external files
#'
#' Function to read in, from external files, the output of software for ROH:
#' \enumerate{
#' \item \code{detectRUNS}: output saved out to a file (e.g. write.table)
#' \item \code{Plink}: output from the \code{--homozyg} option (\code{.hom} files)
#' \item \code{BCFtools}: output from the \code{roh} option
#' }
#'
#' @param inputFile name of (path to) external file
#' @param program source program that produced the ROH file (one of \code{detectRUNS},
#' \code{Plink}, \code{BCFtools})
#'
#' @return dataframe in the correct format to be used with plots and statistics functions from \code{detectRUNS}
#' @export
#'
#' @examples
#' # getting map and ped paths
#' \dontrun{
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- scanRUNS(genotypeFile, method = "sliding", windowSize = 15, threshold = 0.1, minSNP = 15,
#' ROHet = FALSE, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#'
#' write.table(x= runs,file = 'Kijas2016_Sheep_subset.sliding.csv', quote=F, row.names = F)
#' }
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package = "detectRUNS")
#' newData=readExternalRuns(runsFile, program = 'detectRUNS')
#' 

readExternalRuns <- function(inputFile=NULL,program=c("plink","BCFtools","detectRUNS")) {

  # check method
  method <- match.arg(program)
  message(paste("Loading file from", method))

  # detectRUNS
  if (method == "detectRUNS"){
    FinalRuns <- read.table(textConnection(gsub("[,\\; \t]", "\t", readLines(inputFile))),
                            header=TRUE,stringsAsFactors = FALSE,
                            colClasses = c(rep("character", 3), rep("numeric", 4)))

    if(ncol(FinalRuns)!=7) stop(paste("Number of colums must be 7! Current n. of cols:",ncol(FinalRuns),sep=" "))
  }

  # plink
  if (method == "plink"){
    plinkDatei <- read.table(file=inputFile, header=TRUE,
                             colClasses = c(rep("character", 6), rep("numeric", 4),rep("character", 3)),
                             col.names =c("group","id","PHE","chrom","SNP1","SNP2","from","to","lengthBps","nSNP","DENSITY","PHOM","PHET"))
    # Subset
    FinalRuns <- plinkDatei[,c("group","id","chrom","nSNP","from","to","lengthBps")]

    #convert kbps to bps
    FinalRuns$lengthBps <- (FinalRuns$lengthBps*1000)
  }

  # BCFtools
  if (method == "BCFtools"){
    subsetBCF <- grep(pattern = "RG", x = readLines(inputFile),invert = FALSE,value = TRUE)
    BCFfinal <- read.table(text=gsub("\t", " ",subsetBCF),header = FALSE,
                           #colClasses = c("character","character","character","numeric","numeric","numeric","numeric"),
                           colClasses = c(rep("character", 3), rep("numeric", 4)),
                           col.names=c("group","id","chrom","from","to","lengthBps","nSNP","Quality")   )
    BCFfinal$chrom <- gsub("Chr", "",BCFfinal$chrom)

    # Final File
    FinalRuns = BCFfinal[,c("group","id","chrom","nSNP","from","to","lengthBps")]
  }

  return(FinalRuns)
}


#' Function to reorder data frames by CHROMOSOME
#'
#' The data frame will be reordered according to chromosome:
#' from 1 to n, then X, Y, XY, MT
#' The data frame needs to have a column with name "CHROMOSOME"
#'
#' @param dfx data frame to be reordered (with column "CHROMOSOME")
#'
#' @details
#' Reorder results based on chromosome
#'
#' @return A reordered data frame by chromosome
#' @export
#'

reorderDF <- function(dfx) {

  chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
  list_chr <- unique(dfx$chrom)
  chr_order <- chr_order[chr_order %in% list_chr]
  #order
  ordered_dfx <- dfx[match(chr_order,dfx$chrom),]

  return(ordered_dfx)
}


###########################################################
### PLINK binary file readers (BIM / FAM)
###########################################################


#' Read a PLINK BIM file
#'
#' Reads a PLINK .bim file and returns a \code{data.table} with one row per
#' SNP.  Useful for inspecting SNP metadata before calling \code{scanRUNS()}.
#'
#' @param bimFile Path to the .bim file.
#' @return A \code{data.table} with columns:
#'   \code{chrom}, \code{snp_id}, \code{cm}, \code{bp_pos}, \code{a1}, \code{a2}.
#' @export
#'
readBimFile <- function(bimFile) {
  if (!file.exists(bimFile))
    stop(paste("BIM file not found:", bimFile))

  bim <- data.table::fread(
    bimFile, header = FALSE,
    colClasses = c("character", "character", "numeric",
                   "integer",   "character", "character")
  )
  data.table::setnames(bim, c("chrom", "snp_id", "cm", "bp_pos", "a1", "a2"))
  return(bim)
}


#' Read a PLINK FAM file
#'
#' Reads a PLINK .fam file and returns a \code{data.table} with one row per
#' sample.  Useful for inspecting sample metadata before calling
#' \code{scanRUNS()}.
#'
#' @param famFile Path to the .fam file.
#' @return A \code{data.table} with columns:
#'   \code{fid}, \code{iid}, \code{pat}, \code{mat}, \code{sex}, \code{pheno}.
#' @export
#'
readFamFile <- function(famFile) {
  if (!file.exists(famFile))
    stop(paste("FAM file not found:", famFile))

  fam <- data.table::fread(
    famFile, header = FALSE,
    colClasses = c("character", "character", "character",
                   "character", "integer",   "numeric")
  )
  data.table::setnames(fam, c("fid", "iid", "pat", "mat", "sex", "pheno"))
  return(fam)
}
