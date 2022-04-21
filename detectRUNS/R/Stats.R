#####################
## STATISTIC FOR RUNS
#####################

#' Function to found max position for each chromosome
#'
#'
#' @param mapFile Plink map file (for SNP position)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @return A data frame with the max position for chromosome
#' @keywords internal
#'

chromosomeLength <- function(mapFile){
  # read mapfile
  mappa <- readMapFile(mapFile)

  maps<-mappa[mappa$POSITION != 0, ] #delete chromosome 0

  # defining NULL variables to avoid warning messages
  CHR <- NULL
  POSITION <- NULL

  # find max value for chromosome
  LengthGenome=ddply(maps,.(CHR),summarize,max(POSITION))
  names<-c("CHROMOSOME","CHR_LENGTH")
  colnames(LengthGenome)<-names

  LengthGenome$CHR_LENGTH = as.numeric(as.vector(LengthGenome$CHR_LENGTH))

  # get total chromosome length
  message(paste("Total genome length:",sum(LengthGenome$CHR_LENGTH),sep=' '))

  return(LengthGenome)
}


#' Function to calculated Froh genome-wide or chromosome-wide
#'
#' This function calculates the individual inbreeding coefficients based on runs of
#' homozygosity (ROH), either per-chromosome (chromosome-wide) or based on the
#' entire genome (genome-wide). See details of calculations below
#'
#' @param runs R object (dataframe) with results on runs
#' @param mapFile Plink map file (to retrieve SNP position)
#' @param genome_wide vector of TRUE/FALSE (genome-wide or chromosome-wide;
#' defaults to TRUE/genome-wide)
#'
#' @details
#' Froh is calculated as:
#'
#' \eqn{ F_{ROH} = \frac{\sum ROH_{length}}{Length_{genome}} }
#'
#' Depending on whether genome-wide or chromosome-wide calculations are required,
#' the terms in the numerator and denominator will refer to the entire genome
#' or will be restricted to specific chromosomes.
#'
#' @import reshape2
#'
#' @return A data frame with the inbreeding coefficients of each individual sample
#'
#' @export
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' \dontrun{
#' # skipping runs calculation
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' Froh_inbreeding(runs = runsDF, mapFile = mapFile)
#' Froh_inbreeding(runs = runsDF, mapFile = mapFile, genome_wide=FALSE)
#'

Froh_inbreeding <- function(runs, mapFile, genome_wide=TRUE){

  LengthGenome=chromosomeLength(mapFile = mapFile)
  info_breed=unique(runs[c('group','id')])

  # Suppress warnings
  id <- NULL
  lengthBps <- NULL
  chrom <- NULL

  #sum of ROH for Sample
  if (genome_wide) {
    message("calculating Froh on all genome")

    # RESULTS!!!!!
    Froh <- ddply(runs,.(id),summarize,sum=sum(lengthBps))
    Froh$Froh_genome =  Froh$sum/sum(LengthGenome$CHR_LENGTH)

  } else {
    message("calculating Froh chromosome by chromosome")

    Froh_temp <- ddply(runs,.(id,chrom),summarize,sum=sum(lengthBps))
    Froh_temp=merge(Froh_temp,LengthGenome,by.y='CHROMOSOME',by.x='chrom')
    Froh_temp$Froh =  Froh_temp$sum/Froh_temp$CHR_LENGTH

    Froh=reshape2::dcast(Froh_temp,id ~ chrom ,value.var = "Froh")

    chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
    list_chr=unique(Froh_temp$chrom)
    new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))
    new_list_chr1=paste("Chr_",new_list_chr,sep="")
    new_list_chr=c("id",new_list_chr)

    # RESULTS!!!!!
    Froh <- Froh[new_list_chr]
    colnames(Froh) <- c('id',new_list_chr1)
  }

  Froh=merge(info_breed,Froh,by="id",all=TRUE)

  return(Froh)
}


#' Function to calculated Froh using a ROH-class
#'
#' This function calculates the individual inbreeding coefficients based on runs of
#' homozygosity (ROH) using only ROH of specific size classes.
#' The parameter \code{class} specify the size interval to split up calculations.
#' For example, if \code{class = 2} Froh based on ROH 0-2, 2-4, 4-8, 80-16, >16 Mbps long
#' will be calculated.
#'
#' @param runs R object (dataframe) with ROH results
#' @param mapFile Plink map file (for SNP position)
#' @param Class base ROH-length interval (in Mbps) (default: 0-2, 2-4, 4-8, 8-16, >16)
#'
#'
#' @return A data frame with individual inbreeding coefficients based on ROH-length of
#' specific size. The sum of ROH-length of specific size in each individual is
#' reported alongside
#' @export
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' \dontrun{
#' # skipping runs calculation
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' Froh_inbreedingClass(runs = runsDF, mapFile = mapFile, Class = 2)
#'

Froh_inbreedingClass <- function(runs, mapFile, Class=2){

  step_value=Class
  range_mb=c(0,0,0,0,0,99999)
  for (i in seq(from = 2 , to= length(range_mb)-1, by = 1) ){
    range_mb[i]=step_value
    step_value=step_value*2
  }

  #range_mb
  name_CLASS=c(paste(range_mb[1],"-",range_mb[2],sep=''),
               paste(range_mb[2],"-",range_mb[3],sep=''),
               paste(range_mb[3],"-",range_mb[4],sep=''),
               paste(range_mb[4],"-",range_mb[5],sep=''),
               paste(">",range_mb[5],sep=''),
               paste(">",range_mb[6],sep=''))

  # Creating the data frame
  runs$MB <- runs$lengthBps/1000000
  runs$CLASS=cut(as.numeric(runs$MB),range_mb)
  levels(runs$CLASS) = name_CLASS
  runs$CLASS=factor(runs$CLASS)
  table(runs$CLASS)

  LengthGenome=chromosomeLength(mapFile)

  # sum of ROH for Sample
  message("calculating Froh by Class")

  # Suppress warnings
  id <- NULL
  lengthBps <- NULL

  Froh_Class=unique(runs[c('group','id')])
  for (i in range_mb[1:5]){
    print(paste("Class used: >",i,sep=''))

    # subset ROHom/ROHet
    subset_roh <- runs[runs$MB >= i,]

    #if subset is empty (no runs for that class) skip/continue
    if(nrow(subset_roh)<1) next

    Froh_temp <- ddply(subset_roh,.(id),summarize,sum=sum(lengthBps))
    Froh_temp[[paste("Froh_Class_",i,sep="")]] =  Froh_temp$sum/sum(LengthGenome$CHR_LENGTH)
    colnames(Froh_temp)[2]<- paste("Sum_Class_",i,sep="")
    Froh_Class=merge(Froh_Class,Froh_temp,by="id",all=TRUE)
  }

  return(Froh_Class)

}


#' Summary statistics on detected runs
#'
#' This function processes the results from \code{slidingRUNS.run} and
#' \code{consecutiveRUNS.run} and produces a number of interesting descriptives
#' statistics on results.
#'
#' @param genotypeFile Plink ped file (for SNP position)
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results on detected runs
#' @param Class group of length (in Mbps) by class (default: 0-2, 2-4, 4-8, 8-16, >16)
#' @param snpInRuns TRUE/FALSE (default): should the function \code{snpInsideRuns} be
#' called to compute the proportion of times each SNP falls inside a run in the
#' group/population?
#'
#' @details
#' \code{summaryRuns} calculates: i) the number of runs per chromosome and group/population;
#' ii) the percent distribution of runs per chromosome and group; iii) the number of
#' runs per size-class and group; iv) the percent distribution of runs per size-class
#' and group; v) the mean length of runs per chromosome and group; vi) the mean
#' length of runs per size-class and group; vii) individual inbreeding coefficient
#' estimated from ROH; viii) individual inbreeding coefficient estimated from ROH
#' per chromosome; ix) individual inbreeding coefficient estimated from ROH per
#' size-class
#'
#' @return A list of dataframes containing the most relevant descriptives
#' statistics on detected runs. The list conveniently contains 9 dataframes that can
#' be used for further processing and visualization, or can be written out to text files
#' @export
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' \dontrun{
#' # skipping runs calculation
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' summaryRuns(runs = runsDF, mapFile = mapFile, genotypeFile = genotypeFile, Class = 2,
#' snpInRuns = FALSE)
#'

summaryRuns <- function(runs, mapFile, genotypeFile, Class=2, snpInRuns=FALSE){
  message("Checking files...")
  message(paste("Using class:",Class))

  # Avoid warnings
  group <- NULL
  CLASS <- NULL
  MB <- NULL
  chrom <- NULL

  n_class=Class

  result_Froh_genome_wide <- Froh_inbreeding(runs = runs,
                                             mapFile = mapFile,
                                             genome_wide = TRUE)
  result_Froh_chromosome_wide <- Froh_inbreeding(runs = runs,
                                                 mapFile = mapFile,
                                                 genome_wide = FALSE)
  result_Froh_class <- Froh_inbreedingClass(runs = runs,
                                                       mapFile = mapFile,
                                                       Class = n_class)


  runs$MB <- runs$lengthBps/1000000
  head(runs)
  #step_value=2

  range_mb <- c(0,0,0,0,0,99999)

  for (i in seq(from = 2 , to= length(range_mb)-1, by = 1) ){
    range_mb[i]=n_class
    n_class=n_class*2
  }

  #range_mb
  name_CLASS=c(paste(range_mb[1],"-",range_mb[2],sep=''),
               paste(range_mb[2],"-",range_mb[3],sep=''),
               paste(range_mb[3],"-",range_mb[4],sep=''),
               paste(range_mb[4],"-",range_mb[5],sep=''),
               paste(">",range_mb[5],sep=''),
               paste(">",range_mb[6],sep=''))

  message(paste("Class created:"  ,name_CLASS[0:5],sep=' '))
  runs$CLASS=cut(as.numeric(runs$MB),range_mb)
  levels(runs$CLASS) = name_CLASS
  runs$CLASS=factor(runs$CLASS)

  #RESULTS!!!!!
  summary_ROH_mean1 = ddply(runs,.(group,CLASS),summarize,sum=mean(MB))
  summary_ROH_mean_class = dcast(summary_ROH_mean1,CLASS ~ group ,value.var = "sum")
  levels(summary_ROH_mean_class$CLASS) = name_CLASS[0:5]

  #RESULTS!!!!!
  summary_ROH_mean_chr1 = ddply(runs,.(group,chrom),summarize,sum=mean(MB))
  summary_ROH_mean_chr = reorderDF(dcast(summary_ROH_mean_chr1,chrom ~ group ,value.var = "sum"))

  #RESULTS!!!!!
  summary_ROH_count =  ddply(runs,.(CLASS,group),nrow)
  summary_ROH_count1=dcast(summary_ROH_count, CLASS ~ group , value.var = "V1")
  rownames(summary_ROH_count1)=summary_ROH_count1$CLASS
  summary_ROH_count1$CLASS=NULL
  summary_ROH_count=summary_ROH_count1
  summary_ROH_percentage= as.data.frame(t(as.data.frame( t(summary_ROH_count)/colSums(summary_ROH_count,na.rm=TRUE))))
  summary_ROH_percentage$CLASS=row.names(summary_ROH_percentage)
  summary_ROH_percentage

  #RESULTS!!!!!
  summary_ROH_count_chr =  ddply(runs,.(chrom,group),nrow)
  summary_ROH_count_chr1=dcast(summary_ROH_count_chr, chrom ~ group , value.var = "V1")
  rownames(summary_ROH_count_chr1)=summary_ROH_count_chr1$chrom
  summary_ROH_count_chr1$chrom=NULL
  summary_ROH_count_chr=summary_ROH_count_chr1
  summary_ROH_percentage_chr= as.data.frame(t(as.data.frame( t(summary_ROH_count_chr)/colSums(summary_ROH_count_chr,na.rm=TRUE))))
  summary_ROH_percentage_chr$chrom=row.names(summary_ROH_percentage_chr)
  summary_ROH_percentage_chr

  result_summary <- list(summary_ROH_count_chr=summary_ROH_count_chr,
                          summary_ROH_percentage_chr=summary_ROH_percentage_chr,
                          summary_ROH_count=summary_ROH_count,
                          summary_ROH_percentage=summary_ROH_percentage,
                          summary_ROH_mean_chr=summary_ROH_mean_chr,
                          summary_ROH_mean_class=summary_ROH_mean_class,
                          result_Froh_genome_wide = result_Froh_genome_wide,
                          result_Froh_chromosome_wide = result_Froh_chromosome_wide,
                          result_Froh_class= result_Froh_class)

  if (snpInRuns){

    message("Calculating SNPs inside ROH")
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

    #read map file
    mappa <- readMapFile(mapFile)

    #Start calculation % SNP in ROH
    message("Calculation % SNP in ROH") #FILIPPO
    all_SNPinROH <- data.frame("SNP_NAME"=character(),
                               "CHR"=integer(),
                               "POSITION"=numeric(),
                               "COUNT"=integer(),
                               "BREED"=factor(),
                               "PERCENTAGE"=numeric(),
                               stringsAsFactors=FALSE)

    # create progress bar
    total <- length(unique(runs$CHROMOSOME))
    message(paste('Chromosome founds: ',total)) #FILIPPO
    n=0
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom <- mappa[mappa$CHR==chrom,]

      pops <- readPOPCpp(genotypeFile)
      snpInRuns <- snpInsideRunsCpp(runsChrom,mapKrom, pops)

      # remove Number column
      snpInRuns$Number <- NULL

      all_SNPinROH <- rbind.data.frame(all_SNPinROH,snpInRuns)
      n=n+1
      setTxtProgressBar(pb, n)
    }
    close(pb)

    result_summary=append(result_summary,list(SNPinRun = all_SNPinROH))
    message("Calculation % SNP in ROH finish") #FILIPPO
  }



  return(result_summary)
}
