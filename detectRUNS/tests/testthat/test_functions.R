library(detectRUNS)
context("Testing functions")

test_that("Testing snpInRun", {
  # importing data
  data("chillingam")

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

  # get genotype data
  genotype <- chillingham_genotype

  # get map data
  mapFile <- chillingham_map

  # get animals
  animals <- genotype[ ,c(1,2)]

  #remove unnecessary fields from the .raw file
  genotype <- genotype[ ,-c(1:6)]

  # require "plyr"
  n_of_individuals <- vector(length = nrow(genotype))

  # define an internal function
  is_run <- function(x, animal) {
    gaps <- diff(mapFile$bps)
    y <- slidingWindow(x, gaps, windowSize, step=1, ROHet=ROHet, maxOppositeGenotype, maxMiss, maxGap);

    # calculate snpRun (R mode)
    snpRun <- snpInRun(y, windowSize, threshold)

    # calculate snpRun (cpp)
    snpRunCpp <- snpInRunCpp(y, windowSize, threshold)

    # test every record
    expect_identical(snpRun, snpRunCpp)
  }

  for (i in 1:nrow(genotype)) {
    n_of_individuals[i] <- is_run(
      as.vector(genotype[i, ]),
      animal = animals[i, ]
    )
  }

})
