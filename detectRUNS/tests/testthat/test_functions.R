library(testthat)
library(detectRUNS)
context("Testing functions")

test_that("Testing snpInRun", {
  # importing data
  data("chillingam")

  # parameters
  windowSize <- 10
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

  # remove unnecessary fields from the .raw file
  genotype <- genotype[ ,-c(1:6)]

  # require "plyr"
  n_of_individuals <- vector(length = nrow(genotype))

  # calculate gaps
  gaps <- diff(mapFile$bps)

  # define an internal function
  is_run <- function(x) {
    y <- slidingWindow(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # calculate snpRun (R mode)
    snpRun <- snpInRun(y, windowSize, threshold)

    # calculate snpRun (cpp)
    snpRunCpp <- snpInRunCpp(y, windowSize, threshold)

    # test every record
    expect_identical(snpRun, snpRunCpp)
  }

  for (i in 1:nrow(genotype)) {
    n_of_individuals[i] <- is_run( as.integer(genotype[i, ]) )
  }

})

test_that("Testing genoConvert", {
  # create a genotype of 0/1/2
  geno012 <- c(1, 2, 0, 1, 2, NA, 0, NA)
  geno01 <- c(1, 0, 0, 1, 0, NA, 0, NA)

  # test R genoConvert
  test <- genoConvert(geno012)
  expect_identical(test, geno01)

  # testing Cpp genoConvert
  test <- genoConvertCpp(geno012)
  expect_equal(test, geno01)
})

test_that("Testing slidingWindow", {
  # importing data
  data("chillingam")

  # parameters
  windowSize <- 10
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

  # remove unnecessary fields from the .raw file
  genotype <- genotype[ ,-c(1:6)]

  # define an internal function
  is_run <- function(x) {
    gaps <- diff(mapFile$bps)

    # call R function
    y <- slidingWindow(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # call cppFunction
    test <- slidingWindowCpp(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # testing function
    expect_identical(test, y)
  }

  for (i in 1:nrow(genotype)) {
    is_run( as.integer(genotype[i, ]) )
  }

})

test_that("Testing homoZygotTest", {
  # setting parameters
  maxHom <- 1
  maxMiss <- 1
  maxGap <- 10^6

  x <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  gaps <- c(3721, 3871, 7059, 4486, 7545, 4796, 3043, 9736, 3495, 5051,
            9607, 6555, 11934, 6410, 3415, 1302, 3110, 6609, 3292)
  test <- homoZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_true(test)

  # check Cpp function
  test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_true(test)

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)
  test <- homoZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # revert, and change maxGap
  x[1:2] <- c(0, 0)
  gaps[10] <- maxGap + gaps[10]
  test <- homoZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
            615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)
  test <- homoZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)
})

test_that("Testing heteroZygotTest", {
  # setting parameters
  maxHom <- 1
  maxMiss <- 1
  maxGap <- 10^6

  x <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
         1, 1, 1, 1, 0, 0, 1, 0, 0, 0)
  gaps <- c(4374, 8744, 5123, 14229, 5344, 690, 8566, 5853, 2369, 3638,
            4848, 600, 2333, 976, 2466, 2269, 5411, 6021, 4367)
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
            615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_true(test)

  # check Cpp function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_true(test)

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # revert, and change maxGap
  x[1:2] <- c(0, 1)
  gaps[10] <- maxGap + gaps[10]
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)

  # check Cpp function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
  expect_false(test)
})

