
# Valuating detectRUNS performance

# clean up
rm(list = ls())

# source library
source('helper.R')

# importing libraries
library(detectRUNS)
library(microbenchmark)
library(ggplot2)
library(data.table)

# define a parameters list
parameters <- list(
  windowSize=20,
  threshold=0.1,
  minSNP=5,
  ROHet=TRUE,
  maxOppositeGenotype=1,
  maxMiss=1,
  maxGap=10^6,
  minLengthBps=1000,
  minDensity=1/10,
  maxOppRun=NULL,
  maxMissRun=NULL
)

# how many times perform test
times <- 10

# how many points in X axis
x_points <- 10

# get genotype data
genotypeFile  <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
genotype <- read.table(genotypeFile, sep = " ", header = FALSE, stringsAsFactors = FALSE)

# get only ped data
ped <- genotype[ , -c(1:6)]

# convert into raw data
raw <- t(apply(ped, 1, pedConvertCpp))

# clean unuseful data
rm(list=c("genotype", "ped"))

# read first two columns form genotype
animals <- readPOPCpp(genotypeFile = genotypeFile)

# get only one individual. Get index
# idx <- 1
idx <- which(animals$ID=="H70")

# get an animal
animal <- animals[idx, ]

# remap animal correctly
animal <- list(FID=animal$POP, IID=animal$ID)

# remove unuseful rows
x <- raw[idx, ]

# get map data
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
mapfile <- fread(mapFile, header = F)

# setting colnames
colnames(mapfile) <- c("Chrom","SNP","cM","bps")

# calculate sequence. 11 elements, then remove the first
steps <- ceiling(seq(1, length(x), length.out = (x_points+1) ))[-1]

# calculating runs of Homozygosity
runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,
                        minSNP = 15, ROHet = FALSE,  maxOppositeGenotype = 1,
                        maxMiss = 1,  minLengthBps = 100000,  minDensity = 1/10000)

# fix column names
names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

# a dataframe in which i will store everything
tests <- data.frame(fun=character(), step=integer(), time=integer(), language=character())

# iterate over times steps
for (i in steps) {
  # get a subset
  subset_map <- mapfile[1:i, ]
  subset_genotype <- x[1:i]

  # calculate gaps (only one chromosome)
  gaps <- diff(subset_map$bps)

  ##############################################################################
  # Test Windows

  # debug
  message(paste("Test sliding window: step", i))

  # calculate sliding window
  y <- slidingWindow(subset_genotype, gaps, parameters$windowSize, step=1, parameters$maxGap,
                     parameters$ROHet, parameters$maxOppositeGenotype, parameters$maxMiss)

  test_sliding <- microbenchmark(
    slidingWindow(subset_genotype, gaps, parameters$windowSize, step=1, parameters$maxGap,
                  parameters$ROHet, parameters$maxOppositeGenotype, parameters$maxMiss),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("sliding", times)
  test_step <- rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_sliding$time, language=test_language)
  tests <- rbind(tests, tmp)

  # check cpp slidingWindow
  test_slidingCpp <- microbenchmark(
    slidingWindowCpp(subset_genotype, gaps, parameters$windowSize, step=1, parameters$maxGap,
                     parameters$ROHet, parameters$maxOppositeGenotype, parameters$maxMiss),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("sliding", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingCpp$time, language=test_language)
  tests <- rbind(tests, tmp)

  ##############################################################################
  # Test snpInRun

  # debug
  message(paste("Test snpInRun: step", i))

  # vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
  snpRun <- snpInRun(y$windowStatus, parameters$windowSize, parameters$threshold)

  test_snpInRun <- microbenchmark(
    snpInRun(y$windowStatus, parameters$windowSize, parameters$threshold),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("snpInRun", times)
  test_step = rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_snpInRun$time, language=test_language)
  tests <- rbind(tests, tmp)

  # check cpp snpInRun
  test_snpInRunCpp <- microbenchmark(
    snpInRunCpp(y$windowStatus, parameters$windowSize, parameters$threshold),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("snpInRun", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_snpInRunCpp$time, language=test_language)
  tests <- rbind(tests, tmp)

  ##############################################################################
  # Test slidingRuns

  # debug
  message(paste("Test slidingRuns: step", i))

  test_slidingRuns <- microbenchmark(
    slidingRuns(subset_genotype, animal, subset_map, gaps, parameters, cpp=FALSE),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("slidingRuns", times)
  test_step = rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingRuns$time, language=test_language)
  tests <- rbind(tests, tmp)

  # check cpp slidingRuns
  test_slidingRunsCpp <- microbenchmark(
    slidingRuns(subset_genotype, animal, subset_map, gaps, parameters, cpp=TRUE),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("slidingRuns", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingRunsCpp$time, language=test_language)
  tests <- rbind(tests, tmp)

  ##############################################################################
  # Test consecutiveRuns

  # debug
  message(paste("Test consecutiveRuns: step", i))

  test_consecutiveRuns <- microbenchmark(
    consecutiveRuns(subset_genotype, animal, subset_map, parameters$ROHet, parameters$minSNP,
                    parameters$maxOppositeGenotype, parameters$maxMiss, parameters$minLengthBps,
                    parameters$maxGap),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("consecutiveRuns", times)
  test_step = rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingRuns$time, language=test_language)
  tests <- rbind(tests, tmp)

  # check cpp consecutiveRuns
  test_consecutiveRunsCpp <- microbenchmark(
    consecutiveRunsCpp(subset_genotype, animal, subset_map, parameters$ROHet, parameters$minSNP,
                       parameters$maxOppositeGenotype, parameters$maxMiss, parameters$minLengthBps,
                       parameters$maxGap),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("consecutiveRuns", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_consecutiveRunsCpp$time, language=test_language)
  tests <- rbind(tests, tmp)

  ##############################################################################
  # Test snpInsideRuns

  # debug
  message(paste("Test snpInsideRuns: step", i))

  # get temporary variables
  mappa <- subset_map
  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL

  # snpInsideRuns needs to be launched against a single chromosome. For testing
  # purpose, we will consider mappa and runs as a unique chromosome

  test_snpInsideRuns <- microbenchmark(
    snpInsideRuns(runs, mappa, genotypeFile),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("snpInsideRuns", times)
  test_step = rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_snpInsideRuns$time, language=test_language)
  tests <- rbind(tests, tmp)

  # check cpp snpInsideRuns
  test_snpInsideRunsCpp <- microbenchmark(
    snpInsideRunsCpp(runs, mappa, genotypeFile),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("snpInsideRuns", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_snpInsideRunsCpp$time, language=test_language)
  tests <- rbind(tests, tmp)
}

# as described by http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
testsc <- summarySE(tests, measurevar="time", groupvars=c("fun","step", "language"))

plotGraph <- function(mydata) {
  # Standard error of the mean
  graph <- ggplot(testsc, aes(x=step, y=time, colour=language)) +
    geom_errorbar(aes(ymin=time-se, ymax=time+se), width=.1) +
    geom_line() +
    geom_point()

  formatter1e6 <- function(x){
    x/10e6
  }

  # scaling data to millisecond
  graph <- graph +
    scale_y_continuous(labels = formatter1e6) +
    ylab("time (ms)") +
    xlab("N° of SNPs")

  # plot subgraps
  graph <- graph + facet_wrap(~fun, ncol = 2, scales = "free_y")

  # change labels
  graph <- graph +
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.5))) +
    theme(plot.title = element_text(size = rel(2))) +
    theme_bw()

  # return graph object
  return(graph)
}

# get a graph object
graph <- plotGraph(testsc)

# write graph in a file
png("test_performance.png", width = 1024, height = 768)
print(graph)
dev.off()

