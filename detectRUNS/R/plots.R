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
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#'
#' # plot runs per animal (interactive)
#' plot_Runs(runs, suppressInds=FALSE, savePlots=FALSE, title_prefix="ROHom")
#'

plot_Runs <- function(runs, suppressInds=FALSE, savePlots=FALSE, title_prefix=NULL) {
  # suppress notes
  IND <- NULL
  LENGTH <- NULL
  CHROMOSOME <- NULL
  START <- NULL
  END <- NULL
  POPULATION <- NULL

  # then names runs
  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  chr_order <-c((0:99),"X","Y","XY","MT")
  list_chr=unique(runs$CHROMOSOME)
  new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

  for (chrom in new_list_chr) {

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

    if (suppressInds) optionen <- ggplot2::theme(axis.text.y=element_blank(),
                                                 axis.title.y=element_blank(),axis.ticks.y=element_blank())

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
    p <- p + ggplot2::geom_segment(data=teilsatz,aes(x = START, y = IND, xend = END,
                                                     yend = IND,colour=as.factor(POPULATION)),alpha=alfa, size=grosse)
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
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#'
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#'
#' # plot runs per animal (interactive)
#' plot_StackedRuns(runs, savePlots=FALSE, title_prefix="ROHom")
#'

plot_StackedRuns <- function(runs, savePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  # avoid warnings in testing
  CHROMOSOME <- NULL
  START <- NULL
  END <- NULL

  #select a POPULATION
  for (rasse in unique(runs$POPULATION)){
    print(paste('Current population: ',rasse))
    teilsatz <- subset(runs,runs$POPULATION==rasse)

    chr_order <-c((0:99),"X","Y","XY","MT")
    list_chr=unique(teilsatz$CHROMOSOME)
    new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

    #select a chromosome
    for (chrom in new_list_chr){

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
      p <- p + ggplot2::geom_segment(data=krom, aes(x = START/(10^6), y = ypos, xend = END/(10^6), yend = ypos),
                                     colour="lightcoral", alpha=1, size=0.75)
      p <- p + xlim(0, max(krom$END/(10^6))+10) + ylim(0,length(yread)+1)
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
#' @param genotypeFile genotype (.ped) file location
#' @param mapFile map file (.map) file location
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
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#' # calculating runs of Homozygosity
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#'
#' # plot runs per animal (interactive)
#' plot_SnpsInRuns(runs, genotypeFile, mapFile,
#' savePlots=FALSE, title_prefix="ROHom")

plot_SnpsInRuns <- function(runs, genotypeFile, mapFile, savePlots=FALSE, title_prefix=NULL) {

  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  if(file.exists(mapFile)){
    # using data.table to read data
    mappa <- data.table::fread(mapFile, header = F)
  } else {
    stop(paste("file", mapFile, "doesn't exists"))
  }

  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL

  chr_order <-c((0:99),"X","Y","XY","MT")
  list_chr=unique(runs$CHROMOSOME)
  new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

  # avoid warnings in testing
  CHR <- NULL
  POSITION <- NULL
  PERCENTAGE <- NULL
  BREED <- NULL

  for (chrom in new_list_chr) {

    print(paste("Chromosome is: ",chrom))
    runsChrom <- runs[runs$CHROMOSOME==chrom,]
    print(paste("N. of runs:",nrow(runsChrom)))

    mapKrom <- mappa[mappa$CHR==chrom,]
    print(paste("N.of SNP is",nrow(mapKrom)))

    snpInRuns <- snpInsideRuns(runsChrom,mapKrom, genotypeFile)
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

  plinkDatei <- read.table(file=plinkFile, header=TRUE)
  plinkDatei <- plinkDatei[,c("FID","IID","CHR","NSNP","POS1","POS2","KB")]

  #convert kbps to bps
  plinkDatei$KB <- (plinkDatei$KB*1000)

  #rename columns
  names(plinkDatei) <- c("breed","id","chrom","nSNP","von","bis","lengthBps")

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
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
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
    mapKrom <- mappa[mappa$CHR==chrom,]
    snpInRuns <- snpInsideRuns(runsChrom,mapKrom, genotypeFile)
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
      if (index < chrNum) subset_group[ndx2, "BP"] <- subset_group[ndx2, "BP"] + lstMrk
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
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
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

  names(runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  # avoid warnings
  IND <- NULL
  LENGTH <- NULL
  freq <- NULL
  GROUP <- NULL

  #start calculation by method
  if (method=="sum") {
    message("Using sum")
    sum_ROH_genome <- ddply(runs,.(IND),summarize,sum=sum(LENGTH)/1000000)
    method="Sum"
  } else {
    message("Using mean")
    sum_ROH_genome <- ddply(runs,.(IND),summarize,sum=mean(LENGTH)/1000000)
    method="Mean"
  }

  #sum of ROH for Sample
  count_ROH_genome <- count(runs,"IND")
  sum_ROH_genome=merge(sum_ROH_genome,count_ROH_genome,by='IND')
  sum_ROH_genome=merge(sum_ROH_genome,unique(runs[,c("IND","GROUP")]),by='IND')
  head(sum_ROH_genome)

  #RESULTS!!!!!
  p <- ggplot(data=sum_ROH_genome, aes(x=sum, y=freq, colour=GROUP)) + geom_point()
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
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
#'
#' plot_ViolinRuns(runs, method="sum")
#' plot_ViolinRuns(runs, method="mean")
#'

plot_ViolinRuns <- function(runs, method=c("sum","mean")){

  names(runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  # avoid warnings
  IND <- NULL
  GROUP <- NULL
  LENGTH <- NULL

  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  #start calculation by method
  if (method=="sum") {
    #use ddply
    mean_roh=ddply(runs,.(IND,GROUP),summarize,sum=sum(LENGTH/1000000))
    method="Sum"
  }else{
    mean_roh=ddply(runs,.(IND,GROUP),summarize,sum=mean(LENGTH/1000000))
    method="Mean"
  }

  #Violinplot
  p <- ggplot(data=mean_roh, aes(x=GROUP, y=sum, colour=GROUP))
  p <- p + geom_violin (aes(fill=GROUP)) + geom_boxplot(width=0.1)
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
#' @param runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
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
#' runs <- RUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,  minSNP = 15,
#' ROHet = FALSE,  maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)
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
  p <- p + scale_x_discrete(labels=list_chr)  #+coord_polar()
  p <- p + xlab("Inbreeding by Chromosome") + ylab("Froh")
  p

  if  (polar) {
    p <- p + coord_polar()}

  return(p)
}
