#####################
## STATISTIC FOR RUNS
#####################

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


#' Function to found max position for each chromosome
#'
#'
#' @param mapFile Plink map file (for SNP position)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @return A data frame with the max position for chromosome
#' @export
#'
#' @examples
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' chromosomeLength(mapFile)
#'

chromosomeLength <- function(mapFile){

  # check for file existance
  if(file.exists(mapFile)){
    # using data.table to read data
    mappa <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }

  name<-c("CHR","SNP_NAME","0","POSITION")
  colnames(mappa)<-name
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
#' Froh = (sum of all ROH for individual) / (Genome length covered by SNP)
#'
#'
#' @param runs R object (dataframe) with results per chromosome
#' @param mapFile Plink map file (for SNP position)
#' @param genome_wide vector of TRUE/FALSE (Analisys genome-wide)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @import reshape2
#'
#' @return A data frame with the max position for chromosome
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
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'
#' Froh_inbreeding(runs, mapFile)
#' Froh_inbreeding(runs, mapFile, genome_wide=FALSE)
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
#' @param runs R object (dataframe) with results per chromosome
#' @param mapFile Plink map file (for SNP position)
#' @param Class group of length (in Mbps) by class (defaul: 0-2, 2-4, 4-8, 8-16, >16)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @return A data frame with the max position for chromosome
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
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'
#' Froh_inbreedingClass(runs, mapFile, Class=2)
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


#' Summary for Runs file
#' Report a list data frame for all analisys
#'
#' @param genotypeFile Plink ped file (for SNP position)
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results per chromosome
#' @param Class group of length (in Mbps) by class (defaul: 0-2, 2-4, 4-8, 8-16, >16)
#' @param snpInRuns function to create a dataframe for SNP inside Runs
#'
#' @details
#' Report a list data frame for all analisys
#'
#' @return FILIPPO
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
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'
#' summaryRuns(runs, mapFile, genotypeFile, Class=2, snpInRuns=FALSE)
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
  summary_ROH_count = ddply(runs,.(CLASS,group),nrow)
  names(summary_ROH_count)[3] <- "nRuns"
  summary_ROH_percentage <- summary_ROH_count[,c(1,2)]
  summary_ROH_percentage$pctRuns <- summary_ROH_count$nRuns/sum(summary_ROH_count$nRuns)

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
    if(file.exists(mapFile)){
      # using data.table to read data
      mappa <- data.table::fread(mapFile, header = F)
    } else {
      stop(paste("file", mapFile, "doesn't exists"))
    }
    names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
    mappa$x <- NULL

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
      snpInRuns <- snpInsideRunsCpp(runsChrom,mapKrom, genotypeFile)
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


#' Summary table  for Runs file
#' Report a list data frame for all analisys
#'
#' @param genotypeFile Plink ped file (for SNP position)
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results per chromosome
#' @param threshold value 0 to 1 (default 0.7)
#' @param SnpInRuns dataframe for SNP inside Runs
#'
#' @details
#' Table
#'
#' @return FILIPPO
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
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
#' colClasses = colClasses)
#'
#' tableRuns(runs=runs, genotypeFile= genotypeFile, mapFile= mapFile, threshold = 0.5)
#'

tableRuns <- function(runs=NULL,SnpInRuns=NULL,genotypeFile, mapFile, threshold = 0.5) {

  #set a threshold
  threshold_used=threshold*100
  message(paste('Threshold used:',threshold_used))

  #read map file
  if(file.exists(mapFile)){
    # using data.table to read data
    mappa <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }
  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL


  if(!is.null(runs) & is.null(SnpInRuns)){
    message('I found only Runs data frame. GOOD!')

    #change colnames in runs file
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

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

    #SNP in ROH
    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom <- mappa[mappa$CHR==chrom,]
      snpInRuns <- snpInsideRunsCpp(runsChrom,mapKrom, genotypeFile)
      all_SNPinROH <- rbind.data.frame(all_SNPinROH,snpInRuns)
      n=n+1
      setTxtProgressBar(pb, n)
    }
    close(pb)
    message("Calculation % SNP in ROH finish") #FILIPPO

  }
  else if (is.null(runs) & !is.null(SnpInRuns)) {
    message('I found only SNPinRuns data frame. GOOD!')
    all_SNPinROH=SnpInRuns

  }
  else{
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
        #print(paste("Group:",grp,'- Chr:',chr_old,'- n SNP in Runs:',snp_count)) #FILIPPO
        #print(paste("x=",x,"diff",diff,group_subset[x-1,3],group_subset[x-1,1],"lunghezza:",length(rownames(group_subset))))
        
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

  rownames(final_table) <- seq(1,length(row.names(final_table)))
  return(final_table)
}

