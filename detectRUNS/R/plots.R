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
#' obtained from RUNS.run
#' @param suppressInds shall we suppress individual IDs on the y-axis? (defaults to FALSE)
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
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
#' genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
#' mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotype_path, mapfile_path, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = FALSE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000, minDensity = 1/10)
#'
#' # plot runs per animal (interactive)
#' plotRuns(runs, suppressInds=FALSE, savePlots=FALSE, title_prefix="ROHom")
#'

plotRuns <- function(runs, suppressInds=FALSE, savePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  for (chrom in sort(unique(runs$CHROMOSOME))) {

    #subset by chromosome
    krom <- subset(runs,CHROMOSOME==chrom)

    #rearrange subset
    teilsatz <- krom[,c(5,6,2,1)]
    teilsatz <- teilsatz[order(teilsatz$POPULATION),]

    #new progressive numerical ID
    newID <- seq(1,length(unique(teilsatz$IND)))
    id <- unique(teilsatz$IND)
    teilsatz$NEWID=newID[match(teilsatz$IND,id)]

    optionen <- ggplot2::scale_y_discrete("IDs",limits=unique(teilsatz$IND))
    alfa <- 1
    grosse <- 1

    if (length(id) > 50) {

      optionen <- ggplot2::theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
      alfa <- 0.75
      grosse <- 0.25
    }

    if (suppressInds) optionen <- ggplot2::theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())

    #size in mb
    teilsatz$START <- (teilsatz$START/(10^6))
    teilsatz$END <- (teilsatz$END/(10^6))

    row.names(teilsatz) <- NULL

    teilsatz$IND <- as.factor(teilsatz$IND)
    teilsatz$IND <- factor(teilsatz$IND, levels = unique(teilsatz$IND[order(teilsatz$NEWID)]))

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "chromosome", chrom, sep="_")
    } else {
      titel <- paste("chromosome", chrom, sep="_")
    }

    p <- ggplot2::ggplot(teilsatz)
    p <- p + ggplot2::geom_segment(data=teilsatz,aes(x = START, y = IND, xend = END, yend = IND,colour=as.factor(POPULATION)), alpha=alfa, size=grosse)
    p <- p + ggplot2::xlim(0, max(teilsatz$END)) + ggplot2::ggtitle(paste('Chromosome:',chrom))
    p <- p + ggplot2::guides(colour=guide_legend(title="Population")) + ggplot2::xlab("Mbps")
    p <- p + optionen

    if(savePlots) {
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(p)
      dev.off()
    } else print(p)

  }
}


#' Plot stacked runs
#'
#' Function to plot stacked runs along the chromosome (signalling presence of large numbers of runs)
#' Counts on the y-axis, bps on the x-axis: plots run (TRUE) / no run (FALSE)
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
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
#' genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
#' mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotype_path, mapfile_path, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = FALSE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000, minDensity = 1/10)
#'
#' # plot runs per animal (interactive)
#' plotStackedRuns(runs, savePlots=FALSE, title_prefix="ROHom")
#'

plotStackedRuns <- function(runs, savePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  #select a POPULATION
  for (rasse in unique(runs$POPULATION)){
    print(paste('Current population: ',rasse))
    teilsatz <- subset(runs,runs$POPULATION==rasse)

    #select a chromosome
    for (chrom in sort(unique(teilsatz$CHROMOSOME))){

      print(paste('CHR: ',chrom))
      krom <- subset(teilsatz,CHROMOSOME==chrom)
      krom <- krom[order(krom$START),]

      #start the order
      yread <- c(); #keeps track of the x space that is used up by segments

      # get x axis limits
      minstart <- min(krom$START);
      maxend <- max(krom$END);

      # initialise yread
      yread[1] <- minstart - 1;
      ypos <- c(); #holds the y pos of the ith segment

      for (r in 1:nrow(krom)){
        read <- krom[r,];
        start <- read$START;
        placed <- FALSE;

        # iterate through yread to find the next availible
        # y pos at this x pos (start)
        y <- 1;
        while(!placed){

          if(yread[y] < start){
            ypos[r] <- y;
            yread[y] <- read$END;
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
        titel <- paste(title_prefix, "chr", chrom, rasse, "stacked", sep="_")
      } else {
        titel <- paste("chr", chrom, rasse, "stacked", sep="_")
      }

      #PLOT STACKED RUNS
      p <- ggplot2::ggplot()
      p <- p + ggplot2::geom_segment(data=krom, aes(x = START/(10^6), y = ypos, xend = END/(10^6), yend = ypos), colour="lightcoral", alpha=1, size=0.75)
      p <- p + xlim(0, max(krom$END/(10^6))+10) + ylim(0,length(yread))
      p <- p + ylab('n Runs') + xlab('Chromosome position (Mbps)')
      p <- p + ggplot2::ggtitle(paste("POPULATION: ",rasse,'\nChromosome:',chrom))

      if(savePlots) {
        pdf(paste(titel,".pdf",sep=""),height=8,width=10)
        print(p)
        dev.off()
      } else print(p)

    }
  }

}

#' Plot N. of times SNP is in runs
#'
#' Function to plot the number of times/percentage a SNP in in a run (population-specific signals)
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param genotype_path genotype (.ped) file location
#' @param mapfile_path map file (.map) file location
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
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
#' genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
#' mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotype_path, mapfile_path, windowSize = 20, threshold = 0.1, minSNP = 5,
#' ROHet = FALSE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000, minDensity = 1/10)
#'
#' # plot runs per animal (interactive)
#' plotSnpsInRuns(runs, genotype_path, mapfile_path, savePlots=FALSE, title_prefix="ROHom")
#'


plotSnpsInRuns <- function(runs, genotype_path, mapfile_path, savePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  if(file.exists(mapfile_path)){
    # using data.table to read data
    mappa <- data.table::fread(mapfile_path, header = F)
  } else {
    stop(paste("file", mapfile_path, "doesn't exists"))
  }

  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL

  for (chrom in sort(unique(runs$CHROMOSOME))) {

    print(paste("Chromosome is: ",chrom))
    runsChrom <- runs[runs$CHROMOSOME==chrom,]
    print(paste("N. of runs:",nrow(runsChrom)))

    mapKrom <- mappa[mappa$CHR==chrom,]
    print(paste("N.of SNP is",nrow(mapKrom)))

    snpInRuns <- snp_inside_ROH(runsChrom,mapKrom, genotype_path)
    krom <- subset(snpInRuns,CHR==chrom)

    p <- ggplot(data=krom, aes(x=POSITION/(10^6), y=PERCENTAGE, colour=BREED))
    p <- p + geom_line() +  ggtitle(paste('chr', chrom, sep=' '))
    p <- p + scale_y_continuous(limits = c(-0, 100)) + xlab("Mbps")
    p <- p + scale_x_continuous(limits = c(-0, max(snpInRuns$POSITION/(10^6))+1))

    if (! is.null(title_prefix)) {
      titel <- paste(title_prefix, "chr", chrom, "SNP", sep="_")
    } else {
      titel <- paste("chr", chrom, "SNP", sep="_")
    }

    if(savePlots) {
      pdf(paste(titel,".pdf",sep=""),height=8,width=10)
      print(p)
      dev.off()
    } else print(p)

  }

}
