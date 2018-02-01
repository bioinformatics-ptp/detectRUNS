####################################
## FUNCTIONS TO MAKE PLOTS FROM RUNS
####################################

# three functions for three different plots:
# i) Runs per individual sample (samples on the y-axis)
# ii) stacked runs per sample
# iii) n. of times a SNP is in a run in the population

#' Function to plot runs per individual
#'
#' Function to plot runs per individual (see Williams et al. 2016, Animal Genetics,
#' for an example with animal data)
#' Individual IDs on the y-axis, bps on the x-axis (position along the chromosome)
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param suppressInds shall we suppress individual IDs on the y-axis? (defaults to FALSE)
#' @param savePlots should plots be saved out to files (one pdf file for all chromosomes)
#' or plotted in the graphical terminal (default)?
#' @param separatePlots should plots for each chromosome be saved out to separate files?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
#'
#' @return plot of runs by chromosome
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' # plot runs per animal (interactive)
#' plot_Runs(runs = runsDF, suppressInds = FALSE, savePlots = FALSE, outputName = "ROHom")
#'

plot_Runs <- function(runs, suppressInds=FALSE, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {

  # avoid notes
  chrom <- NULL ; from <- NULL ; to <- NULL ; group <- NULL

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
    p <- p + ggplot2::xlim(0, max(teilsatz$to)) + ggplot2::ggtitle(paste('Chromosome ',chromosome))
    p <- p + ggplot2::guides(colour=guide_legend(title="Population")) + ggplot2::xlab("Mbps")
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + optionen


    if(savePlots & separatePlots) {
      if (! is.null(outputName)) {
        fileNameOutput <- paste(outputName, "Chr", chromosome, '.pdf', sep="_")
      } else { fileNameOutput <- paste("Chr", chromosome, '.pdf',sep="_") }
      ggsave(filename = fileNameOutput , plot = p, device = "pdf")
    } else if (savePlots) { plot_list[[chromosome]] <- p
    } else { print(p)}

  }

  if(savePlots & !separatePlots) {

    if (! is.null(outputName)) {
      fileNameOutput <- paste(outputName, "AllChromosomes.pdf", sep="_")
    } else {
      fileNameOutput <- "Runs_AllChromosome.pdf"
    }
      plot_list_final <- gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1)
      ggsave(filename = fileNameOutput , plot = plot_list_final, device = "pdf")
  }
}


#' Plot stacked runs
#'
#' Function to plot stacked runs along the chromosome (signaling presence of large numbers of runs
#' in specific regions of a chromosome)
#' Counts on the y-axis, bps on the x-axis (position along the chromosome)
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param savePlots should plots be saved out in files (default) or plotted in
#' the graphical terminal?
#' @param separatePlots should plots for chromosomes be saved out to separate files?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' # plot runs per animal (interactive)
#' plot_StackedRuns(runs = runsDF, savePlots = FALSE, outputName = "ROHom")
#'

plot_StackedRuns <- function(runs, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {

  # avoid notes
  chrom <- NULL
  from <- NULL
  to <- NULL

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

      #PLOT STACKED RUNS
      p <- ggplot2::ggplot()
      p <- p + ggplot2::geom_segment(data=krom, aes(x = from/(10^6), y = ypos, xend = to/(10^6), yend = ypos),
                                     colour="lightcoral", alpha=1, size=0.75)
      p <- p + xlim(0, max(krom$to/(10^6))+10) + ylim(0,length(yread)+1)
      p <- p + ylab('n Runs') + xlab('Chromosome position (Mbps)')
      p <- p + ggplot2::ggtitle(paste("Group: ",rasse,'\nChromosome:',chromosome))
      p <- p + theme(plot.title = element_text(hjust = 0.5))


      # Save plots by Chromosome
      if(savePlots & separatePlots) {
        if (! is.null(outputName)) { fileNameOutput <- paste(outputName, "Chr", chromosome, rasse, "Stacked", sep="_")
        } else { fileNameOutput <- paste("Runs_StackedChr", chromosome, rasse,sep="_") }
        ggsave(filename = paste(fileNameOutput,'.pdf',sep='') , plot = p, device = "pdf")
      } else if (savePlots) { plot_list[[chromosome]] <- p
      } else { print(p) }
    }

    # Save plot all Chromosome
    if(savePlots & !separatePlots) {
      if (! is.null(outputName)) { fileNameOutput <- paste(outputName, rasse ,"StackedAllChr.pdf", sep="_")
      } else { fileNameOutput <- paste("Runs_StackedAllChr",rasse,".pdf",sep='_') }
      plot_list_final <- gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1)
      ggsave(filename = fileNameOutput , plot = plot_list_final, device = "pdf")
    }
  }
}


#' Plot the number of times each SNP falls inside runs
#'
#' Function to plot the number of times/percentage each SNP is inside a run
#' (population-specific signals) against the SNP positions in the genome.
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param savePlots should plots be saved out in files (default) or plotted in
#' the graphical terminal?
#' @param separatePlots should plots for each chromosome be saved out to separate files?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
#'
#' @return plot number of times a SNP is in a run by chromosome and population (pdf files)
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' # plot runs per animal (interactive)
#' plot_SnpsInRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile,
#' savePlots = FALSE, outputName = "ROHom")
#'

plot_SnpsInRuns <- function(runs, genotypeFile, mapFile, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {

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
    p <- p + geom_line() +  ggtitle(paste('Chromosome', chromosome, sep=' '))
    p <- p + scale_y_continuous(limits = c(-0, 100)) + xlab("Mbps")
    p <- p + scale_x_continuous(limits = c(-0, max(snpInRuns$POSITION/(10^6))+1))
    p <- p + theme(plot.title = element_text(hjust = 0.5))


    # Save plots by Chromosome
    if(savePlots & separatePlots) {
      if (! is.null(outputName)) { fileNameOutput <- paste(outputName, "Chr", chromosome, "SNPinRuns", sep="_")
      } else { fileNameOutput <- paste("Runs_SNPinRunsChr", chromosome, sep="_") }
      ggsave(filename = paste(fileNameOutput,'.pdf',sep='') , plot = p, device = "pdf")
    } else if (savePlots) { plot_list[[chromosome]] <- p
    } else { print(p) }
  }

  # Save plot all Chromosome
  if(savePlots & !separatePlots) {
    if (! is.null(outputName)) { fileNameOutput <- paste(outputName,"SNPinRunsAllChr.pdf", sep="_")
    } else { fileNameOutput <- paste("SNPinRunsAllChr.pdf",sep='') }
    plot_list_final <- gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1)
    ggsave(filename = fileNameOutput , plot = plot_list_final, device = "pdf")
  }
}


#' Plot the proportion of times SNPs are inside runs - MANHATTAN PLOT
#'
#' Function to plot the proportion of times/percentage each SNP in inside a run
#' (population-specific signals) against SNP position in all chromosomes together
#' Proportions on the y-axis, bps on the x-axis for all analysed chromosomes
#' This is similar to the familiar GWAS Manhattan plot
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param savePlots should plots be saved out in files (default) or plotted in the graphical terminal?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
#' @param plotTitle title in plot (default)
#'
#' @return Manhattan plots of proportion of times SNPs are inside runs,
#' per population (pdf files)
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' # plot runs per animal (interactive)
#' plot_manhattanRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile,
#' savePlots = FALSE, plotTitle = "ROHom")
#'

plot_manhattanRuns <- function(runs, genotypeFile, mapFile, savePlots=FALSE, outputName=NULL, plotTitle=NULL) {

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
    chr_order <-c((0:99),"X","Y","XY","MT","Z","W")
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
    if (! is.null(plotTitle)) {
      main_title <- paste(plotTitle,group,sep = ' - ') # title + group
    } else {
      main_title <- paste("Manhattan Plot - % SNP in Runs for ",group) } # default title

    #Manhattan plot using ggplot2
    print(paste("Creating Manhattan plot for ",group)) #FILIPPO
    p <- ggplot(subset_group)
    p <- p + geom_point(aes(x=BP, y=P, colour=as.factor(CHR)), alpha=2/3)
    p <- p + scale_color_manual(values=rep(c('red','blue'), round(chrNum/2,0)+1))
    p <- p + scale_size(range = c(0.1, 0.1)) + ylim(0,100)
    p <- p + theme_bw(base_size=11) + theme(legend.position='none')
    p <- p + scale_x_continuous(labels=as.character(chroms), breaks=bpMidVec)
    roh_plot <- p + ggtitle(main_title) + xlab('Chromosome') + ylab('% SNP in Runs') + theme(plot.title = element_text(hjust = 0.5))

    #create a title for manhattan plot
    if (! is.null(outputName)) {
      fileNameOutput <- paste(outputName,group,sep = ' - ')
    }
    else {
      fileNameOutput <- paste("Manhattan_Plot_SNP_in_ROH_",group) #FILIPPO
    }

    #Save plot
    if (savePlots){
      ggsave(filename = paste(fileNameOutput,"_",group,".pdf",sep="") , plot = roh_plot, device = "pdf")
    } else { print(roh_plot) }

    print(paste('Manhattan plot created for ',group)) #FILIPPO

  }

}

#' Plot sum of run-lengths (or average run-lengths) against the number of runs per individual
#'
#' Function to plot the sum of run lengths (or the average run length) per individual
#' against the average number of runs per individual. Points can be differentially
#' coloured by group/population. This plot can be useful to identify patterns in
#' the distribution of runs in different groups (e.g. few long runs vs many short runs)
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param mapFile map file (.map) file path
#' @param method "sum" or "mean" of run lengths per individual sample
#' @param savePlots should plots be saved out to files or plotted in the graphical terminal (default)?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)#'
#' @param plotTitle title in plot (default NULL)
#'
#' @return plot of number of runs vs run-length sum/mean per individual sample
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' plot_PatternRuns(runs = runsDF, mapFile = mapFile, method = 'sum')
#' plot_PatternRuns(runs = runsDF, mapFile = mapFile, method = 'mean')
#'

plot_PatternRuns <- function(runs,mapFile,method=c('sum','mean'), outputName = NULL , savePlots = FALSE, plotTitle = NULL){

  # check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  # set title name
  if(!is.null(plotTitle)){
    mainTitle <- paste(plotTitle,method,sep=' - ') # title plot
  }

  # Set output file name
  if(!is.null(outputName) ){
    fileNameOutput <- paste(outputName,'_',method,'.pdf',sep='') # name outputName
  }else{
    fileNameOutput <- paste('RunsPattern_',method,'.pdf',sep='') # name outputName
  }

  # avoid notes
  lengthBps <- NULL
  group <- NULL

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
  if(!is.null(plotTitle)) { p <- p +  ggtitle(mainTitle) + theme(plot.title = element_text(hjust = 0.5)) }

  # Save Plot
  if (savePlots){ ggsave(filename = fileNameOutput , plot = p, device = "pdf") } else { print(p) }

}


#' Violin plot of run length per individual (either sum or mean)
#'
#' Function to produce violin plots of the distribution of runs lengths per group
#' The sum of run lengths, or its average, per individual sample is used to
#' characterize the distribution of runs
#'
#' @param runs a data.frame with runs per individual (group, id, chrom, nSNP, start, end, length)
#' @param method "sum" or "mean" of run lengths per individual samples
#' @param savePlots should plots be saved out to files or plotted in the graphical terminal (default)?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
#' @param plotTitle title in plot (default NULL)
#'
#' @return Violin plot of the distribution of runs-lengths (sum or mean)
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
#' runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
#'
#' plot_ViolinRuns(runs = runsDF, method = "sum" , savePlots = FALSE)
#' plot_ViolinRuns(runs = runsDF, method = "mean" , savePlots = FALSE)
#'

plot_ViolinRuns <- function(runs, method=c("sum","mean"), outputName = NULL, plotTitle = NULL , savePlots = FALSE) {

  # Check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  # Set plot title
  if(!is.null(plotTitle)){
    mainTitle <- paste(plotTitle,method,sep=' - ') # title plot
  }

  # Set output file name
  if(!is.null(outputName) ){
    fileNameOutput <- paste(outputName,'_',method,'_ViolinPlot.pdf',sep='') # name outputName
  }else{
    fileNameOutput <- paste('ViolinPlot_',method,'.pdf',sep='') # name outputName
  }

  # Avoid notes
  lengthBps <- NULL
  group <- NULL

  # Start calculation by method
  if (method=="sum") {
    mean_roh=ddply(runs,.(id,group),summarize,sum=sum(lengthBps/10^6))
    method="Sum"
  }else{
    mean_roh=ddply(runs,.(id,group),summarize,sum=mean(lengthBps/10^6))
    method="Mean"
  }

  # Violin Plot
  p <- ggplot(data=mean_roh, aes(x=group, y=sum, colour=group))
  p <- p + geom_violin (aes(fill=group)) + geom_boxplot(width=0.1)
  p <- p + ylab(paste(method," of ROH in Mbps" , sep=''))
  if(!is.null(plotTitle)) { p <- p +  ggtitle(mainTitle) + theme(plot.title = element_text(hjust = 0.5)) }
  if (savePlots){ ggsave(filename = fileNameOutput , plot = p, device = "pdf") } else { print(p) }

}


#' Plot Froh-based inbreeding coefficients by group
#'
#' The function plots the distribution of inbreeding/consanguinity coefficients
#' per chromosome and/or group. Three types of plots can be produces: barplots, boxplots,
#' violin plots. With \code{style="All"} all three plots are produced.
#'
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results on detected runs
#' @param groupSplit plots split by group (defaults to TRUE)
#' @param style type of plot: ChrBarPlot, ChrBoxPlot, FrohBoxPlot, All (all plots)
#' @param savePlots should plots be saved out to files or plotted in the graphical terminal (default)?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)
#' @param plotTitle title in plot (default NULL)
#'
#' @return plots of the distribution of inbreeding by chromosome and group
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
#' plot_InbreedingChr(runs = runsDF, mapFile = mapFile, style='All')
#'

plot_InbreedingChr<- function(runs, mapFile , groupSplit=TRUE, style=c("ChrBarPlot","ChrBoxPlot","FrohBoxPlot","All"),
                              outputName = NULL, plotTitle = NULL , savePlots = FALSE){

  # check method
  method <- match.arg(style)

  Chromosome_Inbreeding=Froh_inbreeding(runs = runs,
                                        mapFile = mapFile,
                                        genome_wide = FALSE)
  Genome_Inbreeding=Froh_inbreeding(runs = runs,
                                    mapFile = mapFile,
                                    genome_wide = TRUE)

  # Set plot title
  if(!is.null(plotTitle)){
    mainTitle1 <- paste(plotTitle,sep='') # title ChrBarPlot
    mainTitle2 <- paste(plotTitle,sep='') # title ChrBoxPlot
    mainTitle3 <- paste(plotTitle,sep='') # title FrohBoxPlot
  }

  # Set output file name
  if(!is.null(outputName) ){
    fileNameOutput1 <- paste(outputName,'_BarPlot.pdf',sep='') # title ChrBarPlot
    fileNameOutput2 <- paste(outputName,'_BoxPlot.pdf',sep='') # title ChrBoxPlot
    fileNameOutput3 <- paste(outputName,'_Froh.pdf',sep='') # title FrohBoxPlot
  }else{
    fileNameOutput1 <- paste('ChrBarPlot','.pdf',sep='') # title ChrBarPlot
    fileNameOutput2 <- paste('ChrBoxPlot','.pdf',sep='') # title ChrBoxPlot
    fileNameOutput3 <- paste('BoxPlot_Froh','.pdf',sep='') # title FrohBoxPlot
  }

  # avoid warnings
  variable <- NULL ; value <- NULL ; group <- NULL ; Froh_genome <- NULL

  #transform data in long format using reshape2
  long_DF=melt(Chromosome_Inbreeding,id.vars = c("id", "group"))
  compact_DF=dcast(long_DF, group ~ variable ,fun.aggregate = mean, na.rm = TRUE)

  #creating list chromosome
  name_val=colnames(compact_DF)
  list_chr=gsub("Chr_","",name_val[2:length(name_val)])

  #final data frame
  final_DF=melt(compact_DF, id.vars = c("group"))

  ########
  # Plot BarPlot, BoxPlot, Froh BoxPlot
  # BarPlot by Chromosome style = ChrBarPlot
  head(final_DF)
  if (style == "ChrBarPlot" | style == "All") {
    g1 <- ggplot(data=final_DF, aes(x=variable, y=value, fill=group))
    g1 <- g1 +  geom_bar(stat="identity", position=position_dodge())
    g1 <- g1 +  scale_x_discrete(labels=list_chr)
    g1 <- g1 +  xlab("Inbreeding by Chromosome") + ylab("Froh")
    if(!is.null(plotTitle)) { g1 <- g1 +  ggtitle(mainTitle1) + theme(plot.title = element_text(hjust = 0.5)) }
    if (groupSplit) { g1 <- g1 + facet_grid(group ~. ) + guides(fill=FALSE) }     # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput1 , plot = g1, device = "pdf") } else { print(g1) }
  }

  # BoxPlot by Chromosome by group - style = ChrBoxPlot
  head(long_DF)
  if (style == "ChrBoxPlot" | style == "All") {
    g2 <- ggplot(data=long_DF, aes(x=variable, y=value, fill=group))
    g2 <- g2 + geom_boxplot()
    g2 <- g2 + scale_x_discrete(labels=list_chr)
    g2 <- g2 + xlab("Inbreeding by Chromosome") + ylab("Froh")
    if(!is.null(plotTitle)) { g2 <- g2 +  ggtitle(mainTitle2) + theme(plot.title = element_text(hjust = 0.5)) }
    if (groupSplit) { g2 <- g2 + facet_grid(group ~. ) + guides(fill=FALSE) }    # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput2 , plot = g2, device = "pdf") } else { print(g2) }
  }

  # BoxPlot Froh by group - style = FrohBoxPlot
  head(Genome_Inbreeding)
  if (style == "FrohBoxPlot" | style == "All") {
    g3 <- ggplot(data=Genome_Inbreeding, aes(x=group, y=Froh_genome, colour=group))
    g3 <- g3 + geom_violin(aes(fill=group))
    g3 <- g3 + geom_boxplot(width=0.1)
    g3 <- g3  + ylab("Froh") # +xlab("group")
    if(!is.null(plotTitle)) { g3 <- g3 +  ggtitle(mainTitle3) + theme(plot.title = element_text(hjust = 0.5)) }
    if (savePlots){ ggsave(filename = fileNameOutput3 , plot = g3, device = "pdf") } else { print(g3) }
  }

}


#' Plot Distribution of runs
#'
#' This function the distribution of runs per group. The average run length per size-class,
#' the average run length per chromosome (and group), the percent distribution of runs
#' per size-class and group, and the proportion of runs per chromosome are plotted.
#' With \code{style="All"} all three plots are produced.
#'
#' @param mapFile Plink map file (for SNP position)
#' @param runs R object (dataframe) with results on detected runs
#' @param groupSplit plots split by group (defaults to TRUE)
#' @param style type of plot: MeanClass, MeanChr, RunsPCT, RunsPCT_Chr, All (all plots)
#' @param savePlots should plots be saved out to files or plotted in the graphical terminal (default)?
#' @param outputName title prefix (the base name of graph, if savePlots is TRUE)#'
#' @param plotTitle title in plot (default NULL)
#' @param Class group of length (in Mbps) by class (default: 0-2, 2-4, 4-8, 8-16, >16)
#'
#' @return plot Distribution Runs
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
#' plot_InbreedingChr(runs = runsDF, mapFile = mapFile, style='All')
#'

plot_DistributionRuns <- function(runs, mapFile , groupSplit=TRUE, style=c("MeanClass","MeanChr","RunsPCT","RunsPCT_Chr","All") ,
                              savePlots=FALSE, outputName=NULL, plotTitle=NULL, Class=2){
  # check method
  method <- match.arg(style)

  # avoid warnings
  group =NULL ; CLASS=NULL ; MB=NULL ; chrom=NULL ; value=NULL

  # Set plot title
  if(!is.null(plotTitle)){
    mainTitle1 <- paste(plotTitle,'Mean Length (Mb) by Class',sep=' - ')      # title MeanClass
    mainTitle2 <- paste(plotTitle,'Mean Length (Mb) by Chromosome',sep=' - ') # title MeanChr
    mainTitle3 <- paste(plotTitle,'Percentage Runs by Class',sep=' - ')       # title RunsPCT
    mainTitle4 <- paste(plotTitle,'Percentage Runs by Chromosome',sep=' - ')  # title RunsPCT_Chr
  }else{
    mainTitle1 <- paste('Mean Length (Mb) by Class',sep=' - ')        # title MeanClass
    mainTitle2 <- paste('Mean Length (Mb) by Chromosome',sep=' - ')   # title MeanChr
    mainTitle3 <- paste('Percentage Runs by Class',sep=' - ')         # title RunsPCT
    mainTitle4 <- paste('Percentage Runs by Chromosome',sep=' - ')    # title RunsPCT_Chr
  }

  # Set output file name
  if(!is.null(outputName) ){
    fileNameOutput1 <- paste(outputName,'_MeanClass.pdf',sep='')      # title MeanClass
    fileNameOutput2 <- paste(outputName,'_MeanChr.pdf',sep='')        # title MeanChr
    fileNameOutput3 <- paste(outputName,'_RunsPCT.pdf',sep='')        # title RunsPCT
    fileNameOutput4 <- paste(outputName,'_RunsPCT_Chr.pdf',sep='')    # title RunsPCT_Chr
  }else{
    fileNameOutput1 <- paste('MeanClass','.pdf',sep='')     # title MeanClass
    fileNameOutput2 <- paste('MeanChr','.pdf',sep='')       # title MeanChr
    fileNameOutput3 <- paste('RunsPercentage','.pdf',sep='')       # title RunsPCT
    fileNameOutput4 <- paste('RunsPercentage_Chr','.pdf',sep='')   # title RunsPCT_Chr
  }


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

  head(runs)

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


  #RESULTS!!!!!
  summary_ROH_count_chr =  ddply(runs,.(chrom,group),nrow)
  summary_ROH_count_chr1=dcast(summary_ROH_count_chr, chrom ~ group , value.var = "V1")
  rownames(summary_ROH_count_chr1)=summary_ROH_count_chr1$chrom
  summary_ROH_count_chr1$chrom=NULL
  summary_ROH_count_chr=summary_ROH_count_chr1
  summary_ROH_percentage_chr= as.data.frame(t(as.data.frame( t(summary_ROH_count_chr)/colSums(summary_ROH_count_chr,na.rm=TRUE))))
  summary_ROH_percentage_chr$chrom=row.names(summary_ROH_percentage_chr)


  ########
  # Plot MeanClass, MeanChr, RunsPCT, RunsPCT_Chr
  # Runs mean by class
  if (style == "MeanClass" | style == "All") {
    long_DF=melt(summary_ROH_mean_class,id.vars = c("CLASS"))
    colnames(long_DF)[colnames(long_DF)=='variable'] <- 'group'
    g1 <- ggplot(data=long_DF, aes(x=CLASS, y=value, fill=group))
    g1 <- g1 + geom_bar(stat="identity", position=position_dodge())
    g1 <- g1 + xlab("Class Length Category") + ylab("Mean (Mb)") + scale_x_discrete(limits=unique(long_DF$chrom))
    g1 <- g1 +  ggtitle(mainTitle1) + theme(plot.title = element_text(hjust = 0.5))
    if (groupSplit) { g1 <- g1 + facet_grid(group ~. ) + guides(fill=FALSE) }    # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput1 , plot = g1, device = "pdf") } else { print(g1) }
  }

  # Runs Mean by Chromosome
  if (style == "MeanChr" | style == "All") {
    summary_ROH_mean_chr=reorderDF(summary_ROH_mean_chr)
    long_DF=melt(summary_ROH_mean_chr,id.vars = c("chrom"))
    colnames(long_DF)[colnames(long_DF)=='variable'] <- 'group'
    g2 <- ggplot(data=long_DF, aes(x=chrom, y=value, fill=group))
    g2 <- g2 + geom_bar(stat="identity", position=position_dodge())
    g2 <- g2 + xlab("Chromosome") + ylab("Mean (Mb)") + scale_x_discrete(limits=unique(long_DF$chrom))
    g2 <- g2 +  ggtitle(mainTitle2) + theme(plot.title = element_text(hjust = 0.5))
    if (groupSplit) { g2 <- g2 + facet_grid(group ~. ) + guides(fill=FALSE) }    # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput2 , plot = g2, device = "pdf") } else { print(g2) }
  }

  # Runs percentage by Class
  if (style == "RunsPCT" | style == "All") {
    long_DF=melt(summary_ROH_percentage,id.vars = c("CLASS"))
    colnames(long_DF)[colnames(long_DF)=='variable'] <- 'group'
    g3 <- ggplot(data=long_DF, aes(x=CLASS, y=value, fill=group))
    g3 <- g3 + geom_bar(stat="identity", position=position_dodge()) + scale_x_discrete(limits=unique(long_DF$CLASS))
    g3 <- g3 + xlab("Class Length Category") + ylab("Frequency")
    g3 <- g3 +  ggtitle(mainTitle3) + theme(plot.title = element_text(hjust = 0.5))
    if (groupSplit) { g3 <- g3 + facet_grid(group ~. ) + guides(fill=FALSE) }    # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput3 , plot = g3, device = "pdf") } else { print(g3) }
  }

  # Runs percentage by Chromosome
  if (style == "RunsPCT_Chr" | style == "All") {
    summary_ROH_percentage_chr = reorderDF(summary_ROH_percentage_chr)
    long_DF=melt(summary_ROH_percentage_chr,id.vars = c("chrom"))
    colnames(long_DF)[colnames(long_DF)=='variable'] <- 'group'
    g4 <- ggplot(data=long_DF, aes(x=chrom, y=value, fill=group))
    g4 <- g4 + geom_bar(stat="identity", position=position_dodge()) + scale_x_discrete(limits=unique(long_DF$chrom))
    g4 <- g4 + xlab("Chromosome") + ylab("Frequency")
    g4 <- g4 +  ggtitle(mainTitle4) + theme(plot.title = element_text(hjust = 0.5))
    if (groupSplit) { g4 <- g4 + facet_grid(group ~. ) + guides(fill=FALSE) }    # if you want split or not!
    if (savePlots){ ggsave(filename = fileNameOutput4 , plot = g4, device = "pdf") } else { print(g4) }
  }
}

