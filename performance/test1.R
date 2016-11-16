
# Valuating detectRUNS performance

# clean up
rm(list = ls())

# source library
source('helper.R')

# importing libraries
library(detectRUNS)
library(microbenchmark)
library(ggplot2)
library(bigmemory)
library(data.table)

# parameters
windowSize <- 20
threshold <- 0.1
minSNP <- 5
ROHet <- TRUE
maxOppositeGenotype <- 1
maxMiss <- 1
maxGap <- 10^6
minLengthBps <- 1000
minDensity <- 1/10

# how many times perform test
times <- 100

# how many points in X axis
x_points <- 10

# get genotype data
genotype_path  <- system.file("extdata", "subsetChillingham.raw", package = "detectRUNS")
genotype <- read.big.matrix(genotype_path, sep = " ", header = T, type = "integer")

# load IID and breed
colClasses <- c(
  rep("character", 2),
  rep("NULL", ncol(genotype)-2)
)

animals <- fread(genotype_path, sep = " ", header = T, colClasses = colClasses)

# get only one individual. Get index
idx <- which(animals$IID=="Chill_12")

# remove unuseful columns
x <- genotype[idx, -c(1:6)]

# get map data
mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
mapfile <- fread(mapfile_path, header = F)

# calculate sequence. 11 elements, then remove the first
steps <- ceiling(seq(1, ncol(x)-2, length.out = (x_points+1) ))[-1]

# a dataframe in which i will store everything
tests <- data.frame(fun=character(), step=integer(), time=integer(), language=character())

# iterate over times steps
for (i in steps ) {
  # get a subset
  subset_map <- mapfile$bps[1:i]
  subset_genotype <- x[1:i]

  # calculate gaps (only one chromosome)
  gaps <- diff(subset_map)

  # value diff
  test_diff <- microbenchmark(
    diff(subset_map),
    unit = "ms", #microseconds
    times = times
  )

  test_fun <- rep("diff", times)
  test_step <- rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_diff$time, language=test_language)
  tests <- rbind(tests, tmp)

  # calculate sliding window
  y <- slidingWindow(subset_genotype, gaps, windowSize, step=1, ROHet=ROHet, maxOppositeGenotype, maxMiss, maxGap)

  test_sliding <- microbenchmark(
    slidingWindow(subset_genotype, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss),
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
    slidingWindowCpp(subset_genotype, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("sliding", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_slidingCpp$time, language=test_language)
  tests <- rbind(tests, tmp)


  # vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
  snpRun <- snpInRun(y,windowSize,threshold)

  test_snpInRun <- microbenchmark(
    snpInRun(y,windowSize,threshold),
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
    snpInRunCpp(y,windowSize,threshold),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("snpInRun", times)
  test_step <- rep(i, times)
  test_language <- rep("Cpp", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_snpInRunCpp$time, language=test_language)
  tests <- rbind(tests, tmp)

  # a data.frame with RUNS per animal
  dRUN <- createRUNdf(snpRun, chillingham_map, minSNP, minLengthBps, minDensity)

  test_createRUNdf <- microbenchmark(
    createRUNdf(snpRun, chillingham_map, minSNP, minLengthBps, minDensity),
    unit = 'ms',
    times = times
  )

  test_fun <- rep("createRUNdf", times)
  test_step <- rep(i, times)
  test_language <- rep("R", times)
  tmp <- data.frame(fun=test_fun, step=test_step, time=test_createRUNdf$time, language=test_language)
  tests <- rbind(tests, tmp)
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

