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
#' runs <- scanRUNS(genotypeFile, method = "sliding", windowSize = 15, threshold = 0.1, minSNP = 15,
#' ROHet = FALSE, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
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

  # Use sample_info from RUNS object (includes 0-run individuals);
  # fall back to unique individuals in the runs data for plain data.frames.
  if (inherits(runs, "RUNS") && !is.null(runs$sample_info)) {
    info_breed <- as.data.frame(runs$sample_info)[, c("group", "id"), drop = FALSE]
  } else {
    info_breed <- unique(.get_runs(runs)[c('group','id')])
  }
  runs <- .get_runs(runs)

  # Early return when no runs are present -- all individuals have Froh = 0
  if (nrow(runs) == 0L) {
    if (genome_wide) {
      message("calculating Froh on all genome")
      info_breed$sum         <- 0
      info_breed$Froh_genome <- 0
    } else {
      message("calculating Froh chromosome by chromosome")
    }
    return(info_breed)
  }

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
#' runs <- scanRUNS(genotypeFile, method = "sliding", windowSize = 15, threshold = 0.1, minSNP = 15,
#' ROHet = FALSE, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
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
#' This function processes the results from \code{scanRUNS} and produces a
#' number of interesting descriptives
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
#' statistics on detected runs. The list contains 9 dataframes:
#' 1) summary_ROH_count_chr: n. of runs per chromosome and breed/group;
#' 2) summary_ROH_percentage_chr: percent distribution of runs per chromosome in each breed/group (sum to 1);
#' 3) summary_ROH_count: n. of runs per size-class (Mb) in each breed/group;
#' 4) summary_ROH_percentage: percent distribution of runs per size-class (Mb) in each breed/group (sum to 1);
#' 5) summary_ROH_mean_chr: average size of runs (Mb) per chromosome and breed/group;
#' 6) summary_ROH_mean_class: average size of runs (Mb) per size-class (Mb) in each breed/group;
#' 7) result_Froh_genome_wide: genome-wide inbreeding (\eqn{F_{ROH}}) for each individual;
#' 8) result_Froh_chromosome_wide: inbreeding (\eqn{F_{ROH}}) per individual and chromosome;
#' 9) result_Froh_class: genome-wide inbreeding (\eqn{F_{ROH}}) per individual and size-class (Mb) of runs.
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
#' runs <- scanRUNS(genotypeFile, method = "sliding", windowSize = 15, threshold = 0.1, minSNP = 15,
#' ROHet = FALSE, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' summaryRuns(runs = runsDF, mapFile = mapFile, genotypeFile = genotypeFile, Class = 2,
#' snpInRuns = FALSE)
#'

summaryRuns <- function(runs, mapFile=NULL, genotypeFile=NULL, Class=2, snpInRuns=FALSE){
  runs_input <- runs   # keep original RUNS object for metadata extraction
  runs <- .get_runs(runs)
  message("Checking files...")
  message(paste("Using class:",Class))

  n_class=Class

  result_Froh_genome_wide <- Froh_inbreeding(runs = runs_input, mapFile = mapFile, genome_wide = TRUE)
  result_Froh_chromosome_wide <- Froh_inbreeding(runs = runs_input, mapFile = mapFile, genome_wide = FALSE)
  result_Froh_class <- Froh_inbreedingClass(runs = runs_input, mapFile = mapFile, Class = n_class)

  if (nrow(runs) == 0L) {
    message("No runs detected -- summary statistics are empty.")
    return(list(
      summary_ROH_count_chr       = data.frame(stringsAsFactors = FALSE),
      summary_ROH_percentage_chr  = data.frame(stringsAsFactors = FALSE),
      summary_ROH_count           = data.frame(stringsAsFactors = FALSE),
      summary_ROH_percentage      = data.frame(stringsAsFactors = FALSE),
      summary_ROH_mean_chr        = data.frame(stringsAsFactors = FALSE),
      summary_ROH_mean_class      = data.frame(stringsAsFactors = FALSE),
      result_Froh_genome_wide     = result_Froh_genome_wide,
      result_Froh_chromosome_wide = result_Froh_chromosome_wide,
      result_Froh_class           = result_Froh_class
    ))
  }

  runs$MB <- runs$lengthBps/1000000
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

    chroms <- sort(unique(runs$CHROMOSOME))
    total  <- length(chroms)
    message(paste('Chromosome founds: ', total))
    chrom_list <- vector("list", total)
    n  <- 0L
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    for (chrom in chroms) {
      runsChrom  <- runs[runs$CHROMOSOME == chrom, ]
      mapKrom    <- mappa[mappa$CHR == chrom, ]
      n          <- n + 1L
      chrom_list[[n]] <- snpInsideRuns(runsChrom, mapKrom, sample_info)
      setTxtProgressBar(pb, n)
    }
    close(pb)
    all_SNPinROH <- do.call(rbind, chrom_list)

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
#' @param nCores number of cores for parallel chromosome processing (default: all
#'   physical cores). On Windows only 1 core is used regardless of this value,
#'   as \code{parallel::mclapply} requires a Unix fork.
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
#' runs <- scanRUNS(genotypeFile, method = "sliding", windowSize = 15, threshold = 0.1, minSNP = 15,
#' ROHet = FALSE, maxOpp = 1, maxMiss = 1, minLengthBps = 100000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' runsDF = readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' tableRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile, threshold = 0.5)
#'

tableRuns <- function(runs=NULL, SnpInRuns=NULL, genotypeFile=NULL, mapFile=NULL,
                      threshold=0.5,
                      nCores=parallel::detectCores(logical=FALSE)) {

  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1)
    stop("Threshold must be between 0 and 1")

  runs_input     <- runs
  if (!is.null(runs)) runs <- .get_runs(runs)
  threshold_used <- threshold * 100
  mappa          <- .get_snp_map(runs_input, mapFile)

  if (!is.null(runs) & is.null(SnpInRuns)) {
    message('I found only Runs data frame. GOOD!')
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

    sample_info <- .get_sample_info(runs_input, genotypeFile)
    chroms      <- sort(unique(runs$CHROMOSOME))

    .one_chrom <- function(chrom) {
      runsC <- runs[runs$CHROMOSOME == chrom, ]
      mapC  <- mappa[mappa$CHR == chrom, ]
      if (nrow(runsC) == 0L || nrow(mapC) == 0L)
        return(NULL)

      snp_pct        <- snpInsideRuns(runsC, mapC, sample_info)
      snp_pct$Number <- match(snp_pct$SNP_NAME, mapC$SNP_NAME)

      group_list <- as.vector(unique(snp_pct$BREED))
      grp_tables <- vector("list", length(group_list))

      for (gi in seq_along(group_list)) {
        grp <- group_list[gi]
        gs  <- snp_pct[snp_pct$BREED == grp &
                       snp_pct$PERCENTAGE > threshold_used, , drop = FALSE]
        if (nrow(gs) == 0L) next

        nums <- gs$Number
        chrs <- gs$CHR
        n    <- nrow(gs)

        is_break  <- c(TRUE, nums[-1L] != nums[-n] + 1L | chrs[-1L] != chrs[-n])
        isl_start <- which(is_break)
        isl_end   <- c(isl_start[-1L] - 1L, n)

        rows <- vector("list", length(isl_start))
        for (j in seq_along(isl_start)) {
          s <- isl_start[j]; e <- isl_end[j]
          rows[[j]] <- data.frame(
            Group     = as.character(gs[e, "BREED"]),
            Start_SNP = gs[s, "SNP_NAME"],
            End_SNP   = gs[e, "SNP_NAME"],
            chrom     = gs[s, "CHR"],
            nSNP      = e - s + 1L,
            from      = gs[s, "POSITION"],
            to        = gs[e, "POSITION"],
            stringsAsFactors = FALSE
          )
        }
        grp_tables[[gi]] <- do.call(rbind, rows)
      }
      do.call(rbind, grp_tables)
    }

    chrom_results <- lapply(chroms, .one_chrom)
    final_table <- do.call(rbind, chrom_results)

  } else if (is.null(runs) & !is.null(SnpInRuns)) {
    message('I found only SNPinRuns data frame. GOOD!')
    all_SNPinROH        <- SnpInRuns
    all_SNPinROH$Number <- seq_len(nrow(all_SNPinROH))

    group_list <- as.vector(unique(all_SNPinROH$BREED))
    grp_tables <- vector("list", length(group_list))

    for (gi in seq_along(group_list)) {
      grp <- group_list[gi]
      gs  <- all_SNPinROH[all_SNPinROH$BREED == grp &
                          all_SNPinROH$PERCENTAGE > threshold_used, , drop = FALSE]
      if (nrow(gs) == 0L) next

      nums <- gs$Number
      chrs <- gs$CHR
      n    <- nrow(gs)

      is_break  <- c(TRUE, nums[-1L] != nums[-n] + 1L | chrs[-1L] != chrs[-n])
      isl_start <- which(is_break)
      isl_end   <- c(isl_start[-1L] - 1L, n)

      rows <- vector("list", length(isl_start))
      for (j in seq_along(isl_start)) {
        s <- isl_start[j]; e <- isl_end[j]
        rows[[j]] <- data.frame(
          Group     = as.character(gs[e, "BREED"]),
          Start_SNP = gs[s, "SNP_NAME"],
          End_SNP   = gs[e, "SNP_NAME"],
          chrom     = gs[s, "CHR"],
          nSNP      = e - s + 1L,
          from      = gs[s, "POSITION"],
          to        = gs[e, "POSITION"],
          stringsAsFactors = FALSE
        )
      }
      grp_tables[[gi]] <- do.call(rbind, rows)
    }
    final_table <- do.call(rbind, grp_tables)

  } else {
    stop('You gave me Runs and SNPinRuns! Please choose one!')
  }

  if (is.null(final_table))
    final_table <- data.frame(Group = character(0), Start_SNP = character(0),
                              End_SNP = character(0), chrom = character(0),
                              nSNP = integer(0), from = integer(0),
                              to = integer(0), stringsAsFactors = FALSE)

  if (nrow(final_table) > 0L)
    rownames(final_table) <- seq_len(nrow(final_table))
  final_table
}

