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
  genotype <- chillingham_genotype[,-c(3,4,5,6)]

  # get only one individual
  x <- genotype[genotype$IID=="Chill_12", ]

  # calc gaps
  gaps <- diff(chillingham_map$bps)

  # calc sliding windows
  y <- slidingWindow(x, gaps, windowSize, step=1, ROHet=ROHet, maxOppositeGenotype, maxMiss, maxGap)

  # calculate snpRun (R mode)
  snpRun <- snpInRun(y, windowSize, threshold)

  # calculate snpRun (cpp)
  snpRunCpp <- snpInRunCpp(y, windowSize, threshold)

  expect_identical(snpRun, snpRunCpp)

})
