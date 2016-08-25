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
#' @return TRUE/FALSE (whether a window is homozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'

homoZygotPrufen <- function(x,gaps,maxHet,maxMiss, maxGap) {

  nHet <- sum(x==1,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  ifelse(!(nHet > maxHet | nMiss > maxMiss | any(gaps > maxGap)), TRUE,FALSE)
}

#' Function to check whether a window is (loosely) heterozygous or not
#'
#' This is a core function. Parameters on how to consider a window heterozygous are here (maxHom, maxMiss)
#'
#' @param x vector of 0/1 genotypes (from genoUmschalten())
#' @param maxHom max n. of homozygous SNP in a heterozygous window
#' @param maxMiss max n. of missing in a window
#'
#' @return TRUE/FALSE (whether a window is heterozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'
heteroZygotPrufen <- function(x,gaps,maxHom,maxMiss,maxGap) {

  nHom <- sum(x==0,na.rm=TRUE)
  nMiss <- sum(is.na(x))
  ifelse(!(nHom > maxHom | nMiss > maxMiss | any(gaps > maxGap)), TRUE,FALSE)
}

#' Function to slide a window over a vector (individual's genotypes)
#'
#' This is a core function. The functions to detect RUNS are slidden over the genome
#'
#' @param x vector of 0/1/2 genotypes
#' @param gaps vector of differences between consecutive positions (gaps) in bps
#' @param window size of window (n. of SNP)
#' @param step by which (how many SNP) is the window slidden
#' @param ROHet shall we detect ROHet or ROHom?
#' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
#' @param maxMiss max. n. of missing SNP
#' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
#'
#' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
#' @export
#'
#' @examples #not yet
#'
#'

schiebeFenster <- function(data, gaps, fenster, step, ROHet=TRUE, maxOppositeGenotype=1, maxMiss=1, maxGap) {

  total <- length(data)
  spots <- seq(from = 1, to = (total - fenster + 1), by = step)
  result <- vector(length = length(spots))
  y <- genoUmschalten(data)

  print(paste("Analysing",ifelse(ROHet,"Runs of Heterozygosity (ROHet)","Runs of Homozygosity (ROHom)"),sep=" "))

  if(ROHet) {

    for(i in 1:length(spots)){
      result[i] <- heteroZygotPrufen(y[spots[i]:(spots[i]+fenster-1)],gaps[spots[i]:(spots[i]+fenster-2)],maxOppositeGenotype,maxMiss,maxGap)
    }
    # to include a shrinking sliding-window at the end of the chromosome/genome, uncomment the following line
    # for(i in (length(spots)+1):total) result[i] <- heteroZygotPrufen(y[seq(i,total)],maxOppositeGenotype,maxMiss)
  } else {

    for(i in 1:length(spots)){
      result[i] <- homoZygotPrufen(y[spots[i]:(spots[i]+fenster-1)],gaps[spots[i]:(spots[i]+fenster-2)],maxOppositeGenotype,maxMiss,maxGap)
    }
    # to include a shrinking sliding-window at the end of the chromosome/genome, uncomment the following line
    #for(i in (length(spots)+1):total) result[i] <- homoZygotPrufen(y[seq(i,total)],maxOppositeGenotype,maxMiss)
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
#' @param mapa Plink-like map file (the R data.frame from RUNS.run)
#' @param minSNP minimun n. of SNP to call a RUN
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#'
#' @return a data.frame with RUNS per animal
#' @export
#'
#' @examples #not yet
#'
#'

createRUNdf <- function(snpRun, mapa, minSNP = 3, minLengthBps = 1000, minDensity = 1/10) {

  #requires itertools

  ## write out map file for subsequent plots (snpInRuns)
  write.table(mapa,file="karte.map",quote=FALSE,row.names=FALSE,col.names=FALSE)

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

  #filters on minimum run length and minimum SNP density
  dL <- dL[dL$lengthBps > minLengthBps,]
  dL$SNPdensity <- (dL$nSNP/dL$lengthBps)*1000 # n. SNP per kbps
  dL <- dL[dL$SNPdensity > minDensity, ]
  dL$SNPdensity <- NULL

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

    report_filename <- paste("detected",ifelse(ROHet,"ROHet","ROHom"),"csv",sep=".")

    if(file.exists(report_filename)){
      append = TRUE
      headers = FALSE
    }
    write.table(
      sep=';', dRUN, file=report_filename, quote=FALSE,
      col.names=headers, row.names=FALSE, append=append
    )
    stato <- TRUE
  } else {

    print(paste("No RUNs found for animal",ind,sep=" "))
    stato <- FALSE
  }
  return(stato)
}


#' Function to count number of times a SNP is in a RUN
#'
#'
#' @param runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
#' @param mapKrom R map object with SNP per chromosome
#' @param popFile file with two columns POPULATION/ID (defaults to gegevens.raw, generated by RUNS.run())
#'
#' @return dataframe with counts per SNP in runs (per population)
#' @export
#'
#' @examples #not yet
#'
#'

snp_inside_ROH <- function(runs, mapKrom, popFile = "gegevens.raw") {

  pops <- read.table(popFile,header=TRUE)
  pops <- pops[,c(1,2)]
  names(pops) <- c("POP","ID")

  M <- data.frame("SNP_NAME"=character(),
                   "CHR"=integer(),
                   "POSITION"=numeric(),
                   "COUNT"=integer(),
                   "BREED"=factor(),
                   "PERCENTAGE"=numeric(),
                   stringsAsFactors=FALSE
  )

    for (ras in unique(runs$POPULATION)) {

      print(paste("Breed is:", ras))
      runsBreed <- runs[runs$POPULATION==ras,]
      nBreed <- nrow(pops[pops$POP==as.character(ras),])
      print(paste("N. of animals of breed",ras,nBreed,sep=" "))

      iPos <- ihasNext(mapKrom$POSITION)
      snpCount <- rep(NA,nrow(mapKrom))

      i <- 1
      while(hasNext(iPos)) {

        pos <- nextElem(iPos)
        inRun <- (pos >= runsBreed$START & pos <= runsBreed$END)
        snpCount[i] <- length(inRun[inRun==TRUE])
        i <- i + 1
      }

      mapKrom$COUNT <- snpCount
      mapKrom$BREED <- rep(ras,nrow(mapKrom))
      mapKrom$PERCENTAGE <- (snpCount/nBreed)*100
      mapKrom <- mapKrom[,c("SNP_NAME","CHR","POSITION","COUNT","BREED","PERCENTAGE")]
      M <- rbind.data.frame(M,mapKrom)
    }

  return(M)
}
