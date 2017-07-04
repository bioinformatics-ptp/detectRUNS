
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
# genotype_path  <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
genotype <- read.table(genotype_path, sep = " ", header = FALSE, stringsAsFactors = FALSE)

# get only ped data
ped <- genotype[ , -c(1:6)]

# convert into raw data
raw <- t(apply(ped, 1, pedConvertCpp))

# clean unuseful data
rm(list=c("genotype", "ped"))

# read first two columns form genotype
animals <- readPOPCpp(genotype_path = genotype_path)

# get only one individual. Get index
# idx <- 1
idx <- which(animals$ID=="Chill_12")

# get an animal
animal <- animals[idx, ]

# remap animal correctly
animal <- list(FID=animal$POP, IID=animal$ID)

# remove unuseful rows
x <- raw[idx, ]

# get map data
# mapfile_path <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
mapfile <- fread(mapfile_path, header = F)

# setting colnames
colnames(mapfile) <- c("Chrom","SNP","cM","bps")

# calculate sequence. 11 elements, then remove the first
steps <- ceiling(seq(1, length(x), length.out = (x_points+1) ))[-1]

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

  # # check cpp consecutiveRuns
  # test_consecutiveRuns <- microbenchmark(
  #   consecutiveRuns(subset_genotype, animal, subset_map, parameters$ROHet, parameters$minSNP,
  #                   parameters$maxOppositeGenotype, parameters$maxMiss, parameters$minLengthBps,
  #                   parameters$maxGap, cpp=TRUE),
  #   unit = 'ms',
  #   times = times
  # )
  #
  # test_fun <- rep("consecutiveRuns", times)
  # test_step <- rep(i, times)
  # test_language <- rep("Cpp", times)
  # tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingRunsCpp$time, language=test_language)
  # tests <- rbind(tests, tmp)

}

# write runs and return TRUE/FALSE if RUNS are written out or not
# this will save always in the same output file
# is_run <- writeRUN(as.character(x$IID),dRUN,ROHet,as.character(x$FID))

# call the whole function
# test_RUN <- microbenchmark(
#   RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
#           ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000,
#           minDensity = 1/10
#   ),
#   unit = 'ms',
#   times = 10
# )

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

