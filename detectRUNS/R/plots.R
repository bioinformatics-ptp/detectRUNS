####################################
## FUNCTIONS TO MAKE PLOTS FROM RUNS
####################################

# three functions for three different plots:
# i) Runs per animal (animals on the y-axis)
# ii) stacked runs per animal
# iii) n. of times a SNP is in a run in the population

#' Function to plot runs per animal
#'
#' Function to plot runs per animal (see Williams et al. 2016, Animal Genetics)
#' IDs on the y-axis, bps on the x-axis: plots run (TRUE) / no run (FALSE)
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param suppressInds shall we suppress individual IDs on the y-axis? (defaults to FALSE)
#' @param savePlots should plots be saved out to files (one pdf file for all chromosomes) or plotted in the graphical terminal (default)?
#' @param separatePlots should plots for each individual chromosome be saved out to separate files?
#' @param title_prefix title prefix (the base name of graph, if savePlots is TRUE)
#'
#' @return plot of runs by chromosome (pdf files)
#' @export
#'
#' @import utils
#' @importFrom grDevices dev.off pdf
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
#' # plot runs per animal (interactive)
#' plot_Runs(runs, suppressInds=FALSE, savePlots=FALSE, title_prefix="ROHom")
#'

plot_Runs <- function(runs, suppressInds=FALSE, savePlots=FALSE, separatePlots=FALSE, title_prefix=NULL) {

  chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
  list_chr=unique(runs$chrom)
  new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

  plot_list <- list()
  for (chromosome in new_list_chr) {

    #subset by chromosome
    krom <- subset(runs,chrom==chromosome)

    #rearrange subset
    teilsatz <- krom[,c(5,6,2,1)]
    teilsatz <- teilsatz[order(teilsatz$group),]

    #new progressive numerical ID
    newID <- seq(1,length(unique(teilsatz$id)))
    id <- unique(teilsatz$id)
    teilsatz$NEWID=newID[match(teilsatz$id,id)]

    optionen <- ggplot2::scale_y_discrete("IDs",limits=unique(teilsatz$id))
    alfa <- 1
    grosse <- 1

    if (length(id) > 50) {

      optionen <- ggplot2::theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
      alfa <- 0.75
      grosse <- 0.25
    }

    if (suppressInds) optionen <- ggplot2::theme(axis.text.y=element_blank(),
                                                 axis.title.y=element_blank(),axis.ticks.y=element_blank())

    #size in mb
    teilsatz$from <- (teilsatz$from/(10^6))
    teilsatz$to <- (teilsatz$to/(10^6))

    row.names(teilsatz) <- NULL

    teilsatz$id <- as.factor(teilsatz$id)
    teilsatz$id <- factor(teilsatz$id, levels = unique(teilsatz$id[order(teilsatz$NEWID)]))

    p <- ggplot2::ggplot(teilsatz)
    p <- p + ggplot2::geom_segment(data=teilsatz,aes(x = from, y = id, xend = to,
                                                     yend = id,colour=as.factor(group)),alpha=alfa, size=grosse)
    p <- p + ggplot2::xlim(0, max(teilsatz$to)) + ggplot2::ggtitle(paste('Chromosome:',chromosome))
    p <- p + ggplot2::guides(colour=guide_legend(title="Population")) + ggplot2::xlab("Mbps")
    p <- p + optionen

    if(savePlots) {
      plot_list[[chromosome]] <- p
    } else print(p)
  }

  # if (! is.null(title_prefix)) {
  #   titel <- paste(title_prefix, "chromosome", chromosome, sep="_")
  # } else {
  #   titel <- paste("chromosome", chromosome, sep="_")
  # }

  if(savePlots) {

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "all_chromosomes", sep="_")
    } else {
      titel <- "all_chromosomes"
    }
      pdf(paste(titel,".pdf",sep=""))
    
      for(p in plot_list) {
        print(p)
      }

      dev.off()
    
  }

  if(savePlots & separatePlots) {
    for(chromosome in names(plot_list)) {
      if (! is.null(title_prefix)) {
        titel <- paste(title_prefix, "chromosome", chromosome, sep="_")
      } else {
        titel <- paste("chromosome", chromosome, sep="_")
      }
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(plot_list[[chromosome]])
      dev.off()
    }
  }
}


#' Plot stacked runs
#'
#' Function to plot stacked runs along the chromosome (signalling presence of large numbers of runs)
#' Counts on the y-axis, bps on the x-axis: plots run (TRUE) / no run (FALSE)
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
#' @param separatePlots should plots for each individual chromosome be saved out to separate files?
#' @param title_prefix title prefix (the base name of graph, if savePlots is TRUE)
#'
#' @return plot of stacked runs by population and by chromosome (pdf files)
#' @export
#'
#' @import utils
#' @importFrom grDevices dev.off pdf
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
#' # plot runs per animal (interactive)
#' plot_StackedRuns(runs, savePlots=FALSE, title_prefix="ROHom")
#'

plot_StackedRuns <- function(runs, savePlots=FALSE, separatePlots=FALSE, title_prefix=NULL) {

  plot_list <- list()
  #select a POPULATION
  for (rasse in unique(runs$group)){
    print(paste('Current population: ',rasse))
    teilsatz <- subset(runs,runs$group==rasse)

    chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
    list_chr=unique(teilsatz$chrom)
    new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

    #select a chromosome
    for (chromosome in new_list_chr){

      print(paste('CHR: ',chromosome))
      krom <- subset(teilsatz,chrom==chromosome)
      krom <- krom[order(krom$from),]

      #start the order
      yread <- c(); #keeps track of the x space that is used up by segments

      # get x axis limits
      minstart <- min(krom$from);
      maxend <- max(krom$to);

      # initialise yread
      yread[1] <- minstart - 1;
      ypos <- c(); #holds the y pos of the ith segment

      for (r in 1:nrow(krom)){
        read <- krom[r,];
        start <- read$from;
        placed <- FALSE;

        # iterate through yread to find the next availible
        # y pos at this x pos (start)
        y <- 1;
        while(!placed){

          if(yread[y] < start){
            ypos[r] <- y;
            yread[y] <- read$to;
            placed <- TRUE;
          }

          # current y pos is used by another segment, increment
          y <- y + 1;
          # initialize another y pos if we're at the end of the list
          if(y > length(yread)){
            yread[y] <- minstart-1;
          }
        }
      }

      maxy <- length(yread);
      krom$ypos <- ypos;
      utils::head(krom)

      if (! is.null(title_prefix)) {
        titel <- paste(title_prefix, "chr", chromosome, rasse, "stacked", sep="_")
      } else {
        titel <- paste("chr", chromosome, rasse, "stacked", sep="_")
      }

      #PLOT STACKED RUNS
      p <- ggplot2::ggplot()
      p <- p + ggplot2::geom_segment(data=krom, aes(x = from/(10^6), y = ypos, xend = to/(10^6), yend = ypos),
                                     colour="lightcoral", alpha=1, size=0.75)
      p <- p + xlim(0, max(krom$to/(10^6))+10) + ylim(0,length(yread)+1)
      p <- p + ylab('n Runs') + xlab('Chromosome position (Mbps)')
      p <- p + ggplot2::ggtitle(paste("POPULATION: ",rasse,'\nChromosome:',chromosome))

      if(savePlots) {
        plot_list[[paste(rasse,chromosome,sep="_")]] <- p
      } else print(p)

    }
  }
  if(savePlots) {

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "all_chromosomes_stacked", sep="_")
    } else {
      titel <- "all_chromosomes_stacked"
    }
      pdf(paste(titel,".pdf",sep=""))

      for(p in plot_list) {
        print(p)
      }

      dev.off()
    
  }

  if(savePlots & separatePlots) {
    for(chromosome in names(plot_list)) {
      if (! is.null(title_prefix)) {
        titel <- paste(title_prefix, "stacked_chromosome", chromosome, sep="_")
      } else {
        titel <- paste("stacked_chromosome", chromosome, sep="_")
      }
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(plot_list[[chromosome]])
      dev.off()
    }
  }
}


#' Plot N. of times SNP is in runs
#'
#' Function to plot the number of times/percentage a SNP in in a run (population-specific signals)
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param genotypeFile genotype (.ped) file location
#' @param mapFile map file (.map) file location
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
#' @param separatePlots should plots for each individual chromosome be saved out to separate files?
#' @param title_prefix title prefix (the base name of graph, if savePlots is TRUE)
#'
#' @return plot of n. of times a SNP is in a run by chromosome and population (pdf files)
#' @export
#'
#' @importFrom grDevices dev.off pdf
#' @import utils
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' # skipping runs calculation
#' \dontrun{
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' # loading pre-calculated data
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' colClasses <- c(rep("character", 3), rep("numeric", 4)  )
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses)
#'
#' # plot runs per animal (interactive)
#' plot_SnpsInRuns(runs, genotypeFile, mapFile, savePlots=FALSE, title_prefix="ROHom")

plot_SnpsInRuns <- function(runs, genotypeFile, mapFile, savePlots=FALSE, separatePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  if(file.exists(mapFile)){
    # using data.table to read data
    mappa <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }

  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL

  chr_order <- c((0:99),"X","Y","XY","MT","Z","W")
  list_chr=unique(runs$CHROMOSOME)
  new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

  # avoid warnings in testing
  CHR <- NULL
  POSITION <- NULL
  PERCENTAGE <- NULL
  BREED <- NULL

  plot_list <- list()
  for (chromosome in new_list_chr) {

    print(paste("Chromosome is: ",chromosome))
    runsChrom <- runs[runs$CHROMOSOME==chromosome,]
    print(paste("N. of runs:",nrow(runsChrom)))

    mapChrom <- mappa[mappa$CHR==chromosome,]
    print(paste("N.of SNP is",nrow(mapChrom)))

    snpInRuns <- snpInsideRunsCpp(runsChrom, mapChrom, genotypeFile)
    krom <- subset(snpInRuns,CHR==chromosome)

    p <- ggplot(data=krom, aes(x=POSITION/(10^6), y=PERCENTAGE, colour=BREED))
    p <- p + geom_line() +  ggtitle(paste('chr', chromosome, sep=' '))
    p <- p + scale_y_continuous(limits = c(-0, 100)) + xlab("Mbps")
    p <- p + scale_x_continuous(limits = c(-0, max(snpInRuns$POSITION/(10^6))+1))

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "chr", chromosome, "SNP", sep="_")
    } else {
      titel <- paste("chr", chromosome, "SNP", sep="_")
    }

    if(savePlots) {
      plot_list[[chromosome]] <- p
    } else print(p)

  }
  if(savePlots) {

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "all_chromosomes_snpInRun", sep="_")
    } else {
      titel <- "all_chromosomes_snpInRun"
    }
      pdf(paste(titel,".pdf",sep=""))

      for(p in plot_list) {
        print(p)
      }

      dev.off()
    
  }

  if(savePlots & separatePlots) {
    for(chromosome in names(plot_list)) {
      if (! is.null(title_prefix)) {
        titel <- paste(title_prefix, "snpInRun_chromosome", chromosome, sep="_")
      } else {
        titel <- paste("snpInRun_chromosome", chromosome, sep="_")
      }
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(plot_list[[chromosome]])
      dev.off()
    }
  }
}


#' READ ROH OUTPUT FILE FROM PLINK
#' Function to read in the output file from ROH analysis with Plink
#' Relevant columns are selected, converted and renamed
#'
#' @param plinkFile name of output file from Plink ROH analysis #defaults to plink.hom
#'
#' @return data frame formatted to be used with plot and statistics functions (package detectRUNS)
#' @export
#'
#' @examples #not yet
#'
#'

readFromPlink <- function(plinkFile="plink.hom") {

  plinkDatei <- read.table(file=plinkFile, header=TRUE,
                           colClasses = c("character","character","character","character","character",
                                          "character","numeric","numeric","numeric","numeric","character","character","character" ))
  plinkDatei <- plinkDatei[,c("FID","IID","CHR","NSNP","POS1","POS2","KB")]

  #convert kbps to bps
  plinkDatei$KB <- (plinkDatei$KB*1000)

  #rename columns
  names(plinkDatei) <- c("group","id","chrom","nSNP","from","to","lengthBps")

  return(plinkDatei)
}


#' Plot N. of times SNP is in runs - MANHATTAN PLOT
#'
#' Function to plot the number of times/percentage a SNP in in a run (population-specific signals)
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param genotypeFile genotype (.ped) file location
#' @param mapFile map file (.map) file location
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
#' @param title_prefix title prefix (the base name of graph, if savePlots is TRUE)
#' @param main_titel title in plot
#'
#' @return plot of n. of times a SNP is in a run by chromosome and population (pdf files) using manhattan
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @import utils
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
#' # plot runs per animal (interactive)
#' plot_manhattanRuns(runs, genotypeFile, mapFile, savePlots=FALSE, title_prefix="ROHom")
#'

plot_manhattanRuns <- function(runs, genotypeFile, mapFile, savePlots=FALSE, title_prefix=NULL,main_titel=NULL) {

  #change colnames in runs file
  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  # avoid warnings
  BP <- NULL
  P <- NULL
  CHR <- NULL

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
  print("Calculation % SNP in ROH") #FILIPPO
  all_SNPinROH <- data.frame("SNP_NAME"=character(),
                             "CHR"=integer(),
                             "POSITION"=numeric(),
                             "COUNT"=integer(),
                             "BREED"=factor(),
                             "PERCENTAGE"=numeric(),
                             stringsAsFactors=FALSE)

  # create progress bar
  total <- length(unique(runs$CHROMOSOME))
  print(paste('Chromosome founds: ',total)) #FILIPPO
  n=0
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  for (chrom in sort(unique(runs$CHROMOSOME))) {
    runsChrom <- runs[runs$CHROMOSOME==chrom,]
    mapChrom <- mappa[mappa$CHR==chrom,]
    snpInRuns <- snpInsideRunsCpp(runsChrom,mapChrom, genotypeFile)
    all_SNPinROH <- rbind.data.frame(all_SNPinROH,snpInRuns)
    n=n+1
    setTxtProgressBar(pb, n)
  }
  close(pb)
  print("Calculation % SNP in ROH finish") #FILIPPO

  print("Manhattan plot: START") #FILIPPO
  group_list=unique(all_SNPinROH$BREED)
  for (group in group_list){
    print(paste('Processing Groups:',group)) #FILIPPO

    #list of Groups
    subset_group=subset(all_SNPinROH,all_SNPinROH$BREED==group)
    names(subset_group) <- c("SNP","CHR","BP","COUNT", "GROUP","P")
    subset_group <- subset_group[,c(1,2,3,6)]

    #sort a file
    #subset_group=subset_group[order(as.numeric(subset_group$CHR)),]
    subset_group=subset_group[order(subset_group$CHR),]


    row.names(subset_group) <- 1:nrow(subset_group)

    #create a new position
    chrNum <- length(unique(subset_group$CHR))
    chr_order <-c((0:99),"X","Y","XY","MT")
    list_chr=unique(subset_group$CHR)
    new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))
    chroms = new_list_chr

    for (i in chroms){
      index <- which(chroms==i)
      ndx <- which(subset_group[, "CHR"]==i)
      lstMrk <- max(subset_group[ndx, "BP"])
      if (index < chrNum) ndx2 <- which(subset_group[, "CHR"]==chroms[index+1])
      if (index < chrNum) subset_group[ndx2, "BP"] <- as.numeric(subset_group[ndx2, "BP"] + lstMrk)
    }


    #search a crhomosome center
    bpMidVec <- vector(length=chrNum)
    for (i in chroms){
      ndx <- which((subset_group[, 2])==i)
      posSub <- subset_group[ndx, 3]
      bpMidVec[which(chroms==i)] <- ((max(posSub) - min(posSub))/2) + min(posSub)
    }

    #create a title for manhattan plot
    if (! is.null(main_titel)) {
      main_title <- paste(main_titel,group,sep = ' - ')
    }
    else {
      main_title <- paste("Manhattan Plot - % SNP in ROH for",group) #FILIPPO
    }

    #Manhattan plot using ggplot2
    print(paste("Creating Manhattan plot for ",group)) #FILIPPO
    p <- ggplot(subset_group)
    p <- p + geom_point(aes(x=BP, y=P, colour=as.factor(CHR)), alpha=2/3)
    p <- p + scale_color_manual(values=rep(c('red','blue'), round(chrNum/2,0)+1))
    p <- p + scale_size(range = c(0.1, 0.1)) + ylim(0,100)
    p <- p + theme_bw(base_size=11) + theme(legend.position='none')
    p <- p + scale_x_continuous(labels=as.character(chroms), breaks=bpMidVec)
    #p <- p + geom_hline(yintercept=4.08, linetype=1, col='red', lwd=0.5)  #linea significativa ?? #FILIPPO
    roh_plot <- p + ggtitle(main_title) + xlab('CHROMOSOME') + ylab('% SNP in ROH')

    #main_titel=NULL

    #create a title for manhattan plot
    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix,group,sep = ' - ')
    }
    else {
      titel <- paste("Manhattan_Plot_SNP_in_ROH_",group) #FILIPPO
    }

    #Save plot
    if(savePlots) {
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(roh_plot)
      dev.off()
      print(paste('Manhattan plot created for ',group)) #FILIPPO
    }
    else {
      print(roh_plot)
    }
  }

}


#' Plot N. of ROH by sum/mean
#'
#' Function to plot the number of times/percentage a SNP in in a run (population-specific signals)
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param mapFile map file (.map) file location
#' @param method "sum" or "mean" for single individual
#'
#' @return plot of n. of ROH by sum/mean
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @import utils
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
#' plot_SumMeanRuns(runs, mapFile, method='sum')
#' plot_SumMeanRuns(runs, mapFile, method='mean')
#'

plot_SumMeanRuns <- function(runs,mapFile,method=c('sum','mean')){

  # check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  # checking cromsome lengths
  LengthGenome=chromosomeLength(mapFile)

  # avoid warnings
  freq <- NULL

  #start calculation by method
  if (method=="sum") {
    message("Using sum")
    sum_ROH_genome <- ddply(runs,.(id),summarize,sum=sum(lengthBps)/10^6)
    method="Sum"
  } else {
    message("Using mean")
    sum_ROH_genome <- ddply(runs,.(id),summarize,sum=mean(lengthBps)/10^6)
    method="Mean"
  }

  #sum of ROH for Sample
  count_ROH_genome <- count(runs,"id")
  sum_ROH_genome=merge(sum_ROH_genome,count_ROH_genome,by='id')
  sum_ROH_genome=merge(sum_ROH_genome,unique(runs[,c("id","group")]),by='id')
  head(sum_ROH_genome)

  #RESULTS!!!!!
  p <- ggplot(data=sum_ROH_genome, aes(x=sum, y=freq, colour=group)) + geom_point()
  p <- p + xlab(paste(method," of ROH in Mbps" , sep='')) + ylab ("Number of ROH for Individual")
  p

}


#' Violin plot sum/mean ROH
#'
#' Function to plot violin plot
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param method "sum" or "mean" for single individual
#'
#' @return Violin plot of n. of ROH by sum/mean
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @import utils
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
#' plot_ViolinRuns(runs, method="sum")
#' plot_ViolinRuns(runs, method="mean")
#'

plot_ViolinRuns <- function(runs, method=c("sum","mean")) {

  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  #start calculation by method
  if (method=="sum") {
    #use ddply
    mean_roh=ddply(runs,.(id,group),summarize,sum=sum(lengthBps/10^6))
    method="Sum"
  }else{
    mean_roh=ddply(runs,.(id,group),summarize,sum=mean(lengthBps/10^6))
    method="Mean"
  }

  #Violinplot
  p <- ggplot(data=mean_roh, aes(x=group, y=sum, colour=group))
  p <- p + geom_violin (aes(fill=group)) + geom_boxplot(width=0.1)
  p <- p + ylab(paste(method," of ROH in Mbps" , sep=''))
  p

  return(p)

}


#' Plot Inbreeding by Chromosome
#'
#' The function report a plot with the level fo Froh by chromosome by population
#' It's possible choose the long plot or polar plot.
#'
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results per chromosome
#' @param polar dataframe for SNP inside Runs
#'
#' @return plot Inbreeding by chromosome
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
#' plot_InbreedingChr(runs = runs, mapFile = mapFile, polar=TRUE)
#'

plot_InbreedingChr<- function(runs, mapFile , polar=FALSE){

  Chromosome_Inbreeding=Froh_inbreeding(runs = runs,
                                        mapFile = mapFile,
                                        genome_wide = FALSE)

  # avoid warnings
  variable <- NULL
  value <- NULL
  GROUP <- NULL

  #transform data in long format using reshape2
  long_DF=melt(Chromosome_Inbreeding,id.vars = c("IND", "GROUP"))
  compact_DF=dcast(long_DF, GROUP ~variable ,fun.aggregate = mean, na.rm = TRUE)

  #creating list chromosome
  name_val=colnames(compact_DF)
  list_chr=gsub("Chr_","",name_val[2:length(name_val)])

  #final data frame
  final_DF=melt(compact_DF, id.vars = c("GROUP"))

  #ggplot
  p <- ggplot(data=final_DF, aes(x=variable, y=value, colour=GROUP))
  p <- p + geom_line(aes(group=GROUP))+ geom_point()
  p <- p + scale_x_discrete(labels=list_chr)
  p <- p + xlab("Inbreeding by Chromosome") + ylab("Froh")
  p

  if  (polar) {
    p <- p + coord_polar()}

  return(p)
}
