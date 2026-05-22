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
  mappa <- as.data.frame(readMapFile(mapFile))

  maps <- mappa[mappa$POSITION != 0, ]  # delete chromosome 0

  # find max value for chromosome using base R
  chr_max <- tapply(maps$POSITION, maps$CHR, max)
  LengthGenome <- data.frame(
    CHROMOSOME = names(chr_max),
    CHR_LENGTH = as.numeric(chr_max),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # get total chromosome length
  message(paste("Total genome length:", sum(LengthGenome$CHR_LENGTH), sep=' '))

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

Froh_inbreeding <- function(runs, mapFile=NULL, genome_wide=TRUE){
  LengthGenome <- .get_chrom_lengths(runs, mapFile)
  info_breed   <- unique(.get_runs(runs)[c('group','id')])
  runs         <- .get_runs(runs)

  #sum of ROH for Sample
  if (genome_wide) {
    message("calculating Froh on all genome")

    # aggregate sum of lengthBps by individual id
    sum_by_id <- tapply(runs$lengthBps, runs$id, sum)
    Froh <- data.frame(id = names(sum_by_id), sum = as.numeric(sum_by_id),
                       stringsAsFactors = FALSE)
    Froh$Froh_genome <- Froh$sum / sum(LengthGenome$CHR_LENGTH)

  } else {
    message("calculating Froh chromosome by chromosome")

    key <- paste(runs$id, runs$chrom, sep = "\001")
    sum_by_id_chrom <- tapply(runs$lengthBps, key, sum)
    parts <- strsplit(names(sum_by_id_chrom), "\001", fixed = TRUE)
    Froh_temp <- data.frame(
      id    = vapply(parts, `[[`, "", 1L),
      chrom = vapply(parts, `[[`, "", 2L),
      sum   = as.numeric(sum_by_id_chrom),
      stringsAsFactors = FALSE
    )
    Froh_temp <- merge(Froh_temp, LengthGenome, by.y='CHROMOSOME', by.x='chrom')
    Froh_temp$Froh <- Froh_temp$sum / Froh_temp$CHR_LENGTH

    Froh <- as.data.frame(data.table::dcast(data.table::as.data.table(Froh_temp), id ~ chrom, value.var="Froh"))

    chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
    list_chr <- unique(Froh_temp$chrom)
    new_list_chr <- as.vector(sort(factor(list_chr, levels=chr_order, ordered=TRUE)))
    new_list_chr1 <- paste("Chr_", new_list_chr, sep="")
    new_list_chr <- c("id", new_list_chr)

    # RESULTS!!!!!
    Froh <- Froh[new_list_chr]
    colnames(Froh) <- c('id', new_list_chr1)
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

Froh_inbreedingClass <- function(runs, mapFile=NULL, Class=2){
  LengthGenome <- .get_chrom_lengths(runs, mapFile)
  runs         <- .get_runs(runs)

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

  # sum of ROH for Sample
  message("calculating Froh by Class")

  Froh_Class=unique(runs[c('group','id')])
  for (i in range_mb[1:5]){
    print(paste("Class used: >",i,sep=''))

    # subset ROHom/ROHet
    subset_roh <- runs[runs$MB >= i,]

    #if subset is empty (no runs for that class) skip/continue
    if(nrow(subset_roh)<1) next

    sum_by_id <- tapply(subset_roh$lengthBps, subset_roh$id, sum)
    Froh_temp <- data.frame(id = names(sum_by_id), sum = as.numeric(sum_by_id),
                            stringsAsFactors = FALSE)
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
#' @importFrom stats aggregate
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

summaryRuns <- function(runs, mapFile=NULL, genotypeFile=NULL, Class=2, snpInRuns=FALSE){
  runs_input <- runs   # keep original ROH object for metadata extraction
  runs <- .get_runs(runs)
  message("Checking files...")
  message(paste("Using class:",Class))

  n_class=Class

  result_Froh_genome_wide <- Froh_inbreeding(runs = runs_input, mapFile = mapFile, genome_wide = TRUE)
  result_Froh_chromosome_wide <- Froh_inbreeding(runs = runs_input, mapFile = mapFile, genome_wide = FALSE)
  result_Froh_class <- Froh_inbreedingClass(runs = runs_input, mapFile = mapFile, Class = n_class)


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
  mean1_agg <- aggregate(runs$MB, by=list(group=runs$group, CLASS=runs$CLASS), FUN=mean, na.rm=TRUE)
  colnames(mean1_agg)[3] <- "sum"
  summary_ROH_mean_class <- as.data.frame(data.table::dcast(data.table::as.data.table(mean1_agg), CLASS ~ group, value.var="sum"))
  levels(summary_ROH_mean_class$CLASS) <- name_CLASS[0:5]

  #RESULTS!!!!!
  mean_chr1_agg <- aggregate(runs$MB, by=list(group=runs$group, chrom=runs$chrom), FUN=mean, na.rm=TRUE)
  colnames(mean_chr1_agg)[3] <- "sum"
  summary_ROH_mean_chr <- reorderDF(as.data.frame(data.table::dcast(data.table::as.data.table(mean_chr1_agg), chrom ~ group, value.var="sum")))

  #RESULTS!!!!!
  count_agg <- aggregate(runs$MB, by=list(CLASS=runs$CLASS, group=runs$group), FUN=length)
  colnames(count_agg)[3] <- "V1"
  summary_ROH_count1 <- as.data.frame(data.table::dcast(data.table::as.data.table(count_agg), CLASS ~ group, value.var="V1"))
  rownames(summary_ROH_count1) <- summary_ROH_count1$CLASS
  summary_ROH_count1$CLASS <- NULL
  summary_ROH_count <- summary_ROH_count1
  summary_ROH_percentage <- as.data.frame(t(as.data.frame(t(summary_ROH_count)/colSums(summary_ROH_count, na.rm=TRUE))))
  summary_ROH_percentage$CLASS <- row.names(summary_ROH_percentage)

  #RESULTS!!!!!
  count_chr_agg <- aggregate(runs$MB, by=list(chrom=runs$chrom, group=runs$group), FUN=length)
  colnames(count_chr_agg)[3] <- "V1"
  summary_ROH_count_chr1 <- as.data.frame(data.table::dcast(data.table::as.data.table(count_chr_agg), chrom ~ group, value.var="V1"))
  rownames(summary_ROH_count_chr1) <- summary_ROH_count_chr1$chrom
  summary_ROH_count_chr1$chrom <- NULL
  summary_ROH_count_chr <- summary_ROH_count_chr1
  summary_ROH_percentage_chr <- as.data.frame(t(as.data.frame(t(summary_ROH_count_chr)/colSums(summary_ROH_count_chr, na.rm=TRUE))))
  summary_ROH_percentage_chr$chrom <- row.names(summary_ROH_percentage_chr)

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
    sample_info <- .get_sample_info(runs_input, genotypeFile)
    mappa       <- .get_snp_map(runs_input, mapFile)
    runs        <- runs[, 1:7, drop = FALSE]
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

    all_SNPinROH <- data.frame("SNP_NAME"=character(), "CHR"=integer(),
                               "POSITION"=numeric(), "COUNT"=integer(),
                               "BREED"=factor(), "PERCENTAGE"=numeric(),
                               stringsAsFactors=FALSE)

    total <- length(unique(runs$CHROMOSOME))
    message(paste('Chromosome founds: ', total))
    n <- 0
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom   <- mappa[mappa$CHR==chrom,]
      snp_result <- snpInsideRuns(runsChrom, mapKrom, sample_info)
      all_SNPinROH <- rbind.data.frame(all_SNPinROH, snp_result)
      n <- n + 1
      setTxtProgressBar(pb, n)
    }
    close(pb)

    result_summary <- append(result_summary, list(SNPinRun = all_SNPinROH))
    message("Calculation % SNP in ROH finish")
  }



  return(result_summary)
}


#' Function to retrieve most common runs in the population
#'
#' This function takes in input either the run results or the output from
#' the function \code{snpInsideRuns} (proportion of times a SNP is inside a run)
#' in the population/group, and returns a subset of the runs most commonly
#' found in the group/population. The parameter \code{threshold} controls the definition
#' of most common (e.g. in at least 50\%, 70\% etc. of the sampled individuals)
#'
#' @param genotypeFile Plink ped file (for SNP position)
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results on detected runs
#' @param threshold value from 0 to 1 (default 0.7) that controls the desired
#' proportion of individuals carrying that run (e.g. 70\%)
#' @param SnpInRuns dataframe with the proportion of times each SNP falls inside a
#' run in the population (output from \code{snpInsideRuns})
#'
#' @return A dataframe with the most common runs detected in the sampled individuals
#' (the group/population, start and end position of the run, chromosome and number of SNP
#' included in the run are reported in the output dataframe)
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
#' runsDF = readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' tableRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile, threshold = 0.5)
#'

tableRuns <- function(runs=NULL,SnpInRuns=NULL,genotypeFile=NULL, mapFile=NULL, threshold = 0.5) {
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1)
    stop("Threshold must be between 0 and 1")

  runs_input <- runs
  if (!is.null(runs)) runs <- .get_runs(runs)

  #set a threshold
  threshold_used=threshold*100
  message(paste('Threshold used:',threshold_used))

  # read map file
  mappa <- .get_snp_map(runs_input, mapFile)

  if(!is.null(runs) & is.null(SnpInRuns)){
    message('I found only Runs data frame. GOOD!')

    #change colnames in runs file
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

    #Start calculation % SNP in ROH
    message("Calculation % SNP in ROH")
    all_SNPinROH <- data.frame("SNP_NAME"=character(),
                               "CHR"=integer(),
                               "POSITION"=numeric(),
                               "COUNT"=integer(),
                               "BREED"=factor(),
                               "PERCENTAGE"=numeric(),
                               stringsAsFactors=FALSE)

    # create progress bar
    total <- length(unique(runs$CHROMOSOME))
    message(paste('Chromosome founds: ',total))
    n=0
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    sample_info <- .get_sample_info(runs_input, genotypeFile)

    #SNP in ROH
    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom <- mappa[mappa$CHR==chrom,]
      snpInRuns_result <- snpInsideRuns(runsChrom, mapKrom, sample_info)
      all_SNPinROH <- rbind.data.frame(all_SNPinROH, snpInRuns_result)
      n=n+1
      setTxtProgressBar(pb, n)
    }
    close(pb)
    message("Calculation % SNP in ROH finish")
  } else if (is.null(runs) & !is.null(SnpInRuns)) {
    message('I found only SNPinRuns data frame. GOOD!')
    all_SNPinROH=SnpInRuns
  } else{
    stop('You gave me Runs and SNPinRuns! Please choose one!')
  }

  #consecutive number
  all_SNPinROH$Number <- seq(1,length(all_SNPinROH$PERCENTAGE))

  #final data frame
  final_table <- data.frame("GROUP"=character(0),"Start_SNP"=character(0),"End_SNP"=character(0),
                            "chrom"=character(0),"nSNP"=integer(0),"from"=integer(0),"to"=integer(0))


  #vector of breeds
  group_list=as.vector(unique(all_SNPinROH$BREED))

  for (grp in group_list){
    message(paste('checking: ',grp))

    #create subset for group/thresold
    group_subset=as.data.frame(all_SNPinROH[all_SNPinROH$BREED %in% c(grp) & all_SNPinROH$PERCENTAGE > threshold_used,])

    #print(group_subset)

    #variable
    old_pos=group_subset[1,7]
    snp_pos1=group_subset[1,3]
    Start_SNP=group_subset[1,1]
    snp_count=0

    x=2
    while(x <= length(rownames(group_subset))) {

      snp_count = snp_count + 1
      new_pos=group_subset[x,7]
      old_pos=group_subset[x-1,7]
      chr_old=group_subset[x-1,2]
      chr_new =group_subset[x,2]

      diff=new_pos-old_pos

      if ((diff > 1) | (chr_new != chr_old) | x==length(rownames(group_subset))) {
        if (x==length(rownames(group_subset))){
          end_SNP=group_subset[x,1]
          TO=group_subset[x,3]
        }else{
          end_SNP=group_subset[x-1,1]
          TO=group_subset[x-1,3]
        }

        final_table <- rbind.data.frame(final_table,final_table=data.frame("Group"= group_subset[x-1,5],
                                                                           "Start_SNP"=Start_SNP,
                                                                           "End_SNP"=end_SNP,
                                                                           "chrom"=group_subset[x-1,2],
                                                                           "nSNP"=snp_count,
                                                                           "from"=snp_pos1,
                                                                           "to"=TO))

        #reset variable
        snp_count=0
        snp_pos1=group_subset[x,3]
        Start_SNP=group_subset[x,1]
      }

      #upgrade x value
      x <- x+1

    }
  }

  if (nrow(final_table) > 0L)
    rownames(final_table) <- seq_len(nrow(final_table))
  return(final_table)
}

