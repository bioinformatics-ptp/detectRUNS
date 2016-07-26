#####################
## FUNCTIONS FOR RUNS
#####################


#required external packages
#library("itertools")

#' Function to convert 0/1/2 genotypes to 0/1 (either homozygous/heterozygous)
#'
#' This is a utility function.
#'
#' @param x vector of 0/1/2 genotypes
#'
#' @return converted vector of genotypes
#' @export
#'
#' @examples #not yet
#'
#'
genoUmschalten <- function(x) {

  neues <- c(0,1,0,NA)
  altes <- c(0,1,2,NA)
  return(neues[match(x,altes)])
}

#' Function to check whether a window is (loosely) homozygous or not
#'
#' This is a core function. Parameters on how to consider a window homozygous are here (maxHet, maxMiss)
#'
#' @param x vector of 0/1 genotypes (from genoUmschalten())
#' @param maxHet max n. of heterozygous SNP in a homozygous window
#' @param maxMiss max n. of missing in a window
#'
#' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'

homoZygotPrufen <- function(x,maxHet,maxMiss) {

  nHet <- sum(x==1,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  ifelse(!(nHet > maxHet | nMiss > maxMiss), TRUE,FALSE)
}

#' Function to check whether a window is (loosely) heterozygous or not
#'
#' This is a core function. Parameters on how to consider a window heterozygous are here (maxHom, maxMiss)
#'
#' @param x vector of 0/1 genotypes (from genoUmschalten())
#' @param maxHom max n. of homozygous SNP in a heterozygous window
#' @param maxMiss max n. of missing in a window
#'
#' @return vector of TRUE/FALSE (whether a window is heterozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'
heteroZygotPrufen <- function(x,maxHom,maxMiss) {

  nHom <- sum(x==0,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  ifelse(!(nHom > maxHom | nMiss > maxMiss), TRUE,FALSE)
}

#' Function to slide a window over a vector (individual's genotypes)
#'
#' This is a core function. The functions to detect RUNS are slidden over the genome
#'
#' @param x vector of 0/1/2 genotypes
#' @param window size of window (n. of SNP)
#' @param step by which (how many SNP) is the window slidden
#' @param ROHet shall we detect ROHet or ROHom?
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#'
#' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'

schiebeFenster <- function(data, window, step, ROHet=TRUE, maxOppositeGenotype=1, maxMiss=1) {

  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = total)
  y <- genoUmschalten(data)

  print(paste("Analysing",ifelse(ROHet,"Runs of Heterozygosity (ROHet)","Runs of Homozygosity (ROHom)"),sep=" "))

  if(ROHet) {

    for(i in 1:length(spots)){
      result[i] <- heteroZygotPrufen(y[spots[i]:(spots[i]+window-1)],maxOppositeGenotype,maxMiss)
    }
    for(i in (length(spots)+1):total) result[i] <- heteroZygotPrufen(y[seq(i,total)],maxOppositeGenotype,maxMiss)
  } else {

    for(i in 1:length(spots)){
      result[i] <- homoZygotPrufen(y[spots[i]:(spots[i]+window-1)],maxOppositeGenotype,maxMiss)
    }
    for(i in (length(spots)+1):total) result[i] <- homoZygotPrufen(y[seq(i,total)],maxOppositeGenotype,maxMiss)
  }

  print(paste(
    "Length of homozygous windows overlapping SNP loci (should be equal to the n. of SNP in the file):",
    length(result),sep=" "))

  return(result)
}

#' Function to return a vector of T/F for whether a SNP is or not in a RUN
#'
#' This is a core function. The function to determine whether a SNP is or not in a RUN.
#' The ratio between homozygous/heterozygous windows and total n. of windows is computed here
#'
#' @param RunVector vector of TRUE/FALSE (is a window homozygous/heterozygous?)
#' @param fenster size of window (n. of SNP)
#' @param schwelle threshold to call a SNP in a RUN
#'
#' @return vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
#' @export
#'
#' @examples #not yet
#'
#'
snpInRun <- function(RunVector,fenster,schwelle) {

  total <- length(RunVector)

  print(paste("Length of imput vector:",total,sep=" "))
  print(paste("Window size:",fenster,sep=" "))
  print(paste("Threshold for calling SNP in a Run:",schwelle,sep=" "))

  #requires itertools
  # compute total n. of overlapping windows at each SNP locus (see Bjelland et al. 2013)
  nWin <- c(seq(1,fenster),rep(fenster,(total-(2*fenster))),seq(fenster,1))

  # compute n. of homozygous/heterozygous windows that overlap at each SNP locus (Bjelland et al. 2013)
  iWin <- enumerate(nWin)
  hWin <- sapply(iWin, function(n) sum(RunVector[n$index:(n$index+n$value-1)]), simplify = TRUE)

  # ratio between homozygous/heterozygous windows and total overlapping windows at each SNP
  quotient <- hWin/nWin

  #vector of SNP belonging to a ROH
  snpRun <- ifelse(quotient>schwelle,TRUE,FALSE)

  print(paste(
    "Lenght of output file:",
    length(snpRun),sep=" "))

  return(snpRun)
}

#' Function to create a dataframe of RUNS per individual animal
#' Requires a map file (other filename to read or R object)
#'
#'
#' @param snpRun vector of TRUE/FALSE (is the SNP in a RUN?)
#' @param mapFile Plink map file (either the file path/name or the R data.frame)
#' @param minSNP minimun n. of SNP to call a RUN
#'
#' @return a data.frame with RUNS per animal
#' @export
#'
#' @examples #not yet
#'
#'

createRUNdf <- function(snpRun,mapFile,minSNP = 3) {

  #requires itertools

  mapa <- mapFile

  if(!is.data.frame(mapFile)) {

    if(file.exists(mapFile)){
      mapa <- read.table(mapFile)
    }
  }

  names(mapa) <- c("Chrom","SNP","cM","bps")

  cutPoints <- which(diff(sign(snpRun))!=0)
  von <- c(1,cutPoints+1)
  bis <- c(cutPoints,length(snpRun))

  iLaenge <- izip(a = von,b = bis)
  lengte <- sapply(iLaenge, function(n) sum(snpRun[n$a:n$b]))

  dL <- data.frame("von"=von,"bis"=bis,"nSNP"=lengte)
  dL <- dL[dL$nSNP>minSNP,]
  dL <- na.omit(dL)

  chroms <- mapa[dL$von,"Chrom"]
  dL$von <- mapa[dL$von,"bps"]
  dL$bis <- mapa[dL$bis,"bps"]
  dL$chrom <- as.integer(chroms)
  dL$lengthBps <- (dL$bis-dL$von)

  print(paste("N. of RUNS for this animal","is:",nrow(dL),sep=" "))
  return(dL)
}

#' Function to write out RUNS per individual animal
#'
#'
#' @param ind ID of animal
#' @param dRUN data.frame with RUNS per animal
#'
#' @return TRUE/FALSE if RUNS are written out or not
#' @export
#'
#' @examples #not yet
#'
#'

schreibRUN <- function(ind,dRUN,ROHet=TRUE,breed) {

  dRUN$id <- rep(ind,nrow(dRUN))
  dRUN$breed <- rep(breed,nrow(dRUN))
  dRUN <- dRUN[,c(7,6,4,3,1,2,5)]

  print(paste("N. of RUNS for individual",ind,"is:",nrow(dRUN),sep=" "))

  if(nrow(dRUN) > 0) {

    append = FALSE
    headers = TRUE

    report_filename <- paste("detectRUNS",ifelse(ROHet,"ROHet","ROHom"),"csv",sep=".")

    if(file.exists(report_filename)){
      append = TRUE
      headers = FALSE
    }
    write.table(
      sep=',', dRUN, file=report_filename, quote=FALSE,
      col.names=headers, row.names=FALSE, append=append
    )
    stato <- TRUE
  } else {

    print(paste("No RUNs found for animal",ind,sep=" "))
    stato <- FALSE
  }
  return(stato)
}
