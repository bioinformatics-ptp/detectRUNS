library(testthat)
library(detectRUNS)
context("Testing functions")

# get file paths
genotype_path  <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
raw_path <- system.file("extdata", "subsetChillingham.raw", package = "detectRUNS")

# importing data once
chillingham_genotype <- read.table(genotype_path, sep = " ", header = FALSE, stringsAsFactors = FALSE)
chillingham_map <- read.delim(mapfile_path, header = FALSE)
chillingham_raw <- read.table(raw_path, sep=" ", header = TRUE)

test_that("Testing snpInRun", {
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

  # get genotype data (using raw for tests)
  genotype <- chillingham_raw

  # remove unnecessary fields from the .ped file
  genotype <- genotype[ ,-c(1:6)]

  # get map data
  mapFile <- chillingham_map

  # setting colnames
  names(mapFile) <- c("Chrom","SNP","cM","bps")

  # require "plyr"
  n_of_individuals <- vector(length = nrow(genotype))

  # calculate gaps
  gaps <- diff(mapFile$bps)

  # define an internal function
  is_run <- function(x, i) {
    y <- slidingWindow(x, gaps, windowSize, step=1, maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # calculate snpRun (R mode)
    snpRun <- snpInRun(y$windowStatus, windowSize, threshold)

    # calculate snpRun (cpp)
    snpRunCpp <- snpInRunCpp(y$windowStatus, windowSize, threshold)

    # test every record
    info = paste("step", i)
    expect_identical(snpRun, snpRunCpp, info=info)
  }

  # using raw data for testing functions
  for (i in 1:nrow(genotype)) {
    n_of_individuals[i] <- is_run( as.integer(genotype[i, ]), i )
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

test_that("Testing pedConvert with correct values", {
  # create a PED like genotype
  ped <- c("B", "B", "A", "A", "B", "A", "1", "1", "2", "1", "2", "2", "A", "C", "G", "G", "0", "0", "5", "5", "N", "N", "-", "-")
  geno01 <- c(0, 0, 1, 0, 1, 0, 1, 0, NA, NA, NA, NA)

  # testing Cpp pedConvertCpp
  test <- pedConvertCpp(ped)
  expect_equal(test, geno01)
})

test_that("Testing pedConvert with odd input", {
  # create a PED like genotype
  ped <- c("B", "B", "A")

  # testing Cpp pedConvertCpp
  expect_error(pedConvertCpp(ped), "Need .ped input with 2 allel")
})

test_that("Testing pedConvert with a missing value in a pair", {
  # create ped like genotype
  ped <- c("B", "B", "A", "A", "5", "A", "1", "1", "2", "1", "2", "2", "A", "C", "G", "G", "0", "0", "5", "5", "N", "N")

  # testing pedConvertCpp
  expect_error(test <- pedConvertCpp(ped), "Found only one allele missing in a pair")
})

test_that("Testing data conversion", {
  ped <- chillingham_genotype[ , -c(1:6)]
  raw <-  chillingham_raw[ , -c(1:6)]

  t1 <- apply(ped, 1, pedConvertCpp)
  t2 <- apply(raw, 1, genoConvertCpp)

  # testing conversion
  expect_identical(t1, t2)
})

test_that("Testing slidingWindow", {
  # setting parameters
  windowSize <- 3
  step <- 1
  maxGap <- 1000
  maxOppositeGenotype <- 1
  maxMiss <- 1

  # setting values for RoHet
  data <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA)
  gaps <- rep(100, length(data)-1)

  # update a gap: setting value higher than threshold
  gaps[6] <- gaps[6] + maxGap

  # setting excpected values
  windowStatus <- c(F, F, T, T, F, F, T, T, F, F, T, T, T)
  oppositeAndMissingGenotypes <- c("0", "0", "0", "9", "9", "0", "9")
  names(oppositeAndMissingGenotypes) <- c(1, 2, 3, 10, 11, 13, 15)
  expected <- list(windowStatus=windowStatus,
                   oppositeAndMissingGenotypes=oppositeAndMissingGenotypes)

  # calling function
  test <- slidingWindow(data, gaps, windowSize, step=step, maxGap=maxGap,
                        ROHet=TRUE, maxOppositeGenotype, maxMiss)

  # testing values
  expect_equal(expected, test, info="testing ROHet")

  # setting data for RoHom
  data <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, NA, NA, 0, 1, 0, NA)

  # calling function
  test <- slidingWindow(data, gaps, windowSize, step=step, maxGap=maxGap,
                        ROHet=FALSE, maxOppositeGenotype, maxMiss)

  # testing values
  expect_equal(expected, test, info="testing ROHom")

})

test_that("Testing slidingWindowCpp", {
  # setting parameters
  windowSize <- 3
  step <- 1
  maxGap <- 1000
  maxOppositeGenotype <- 1
  maxMiss <- 1

  # setting values for RoHet
  data <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA)
  gaps <- rep(100, length(data)-1)

  # update a gap: setting value higher than threshold
  gaps[6] <- gaps[6] + maxGap

  # setting excpected values
  windowStatus <- c(F, F, T, T, F, F, T, T, F, F, T, T, T)
  oppositeAndMissingGenotypes <- c("0", "0", "0", "9", "9", "0", "9")
  names(oppositeAndMissingGenotypes) <- c(1, 2, 3, 10, 11, 13, 15)
  expected <- list(windowStatus=windowStatus,
                   oppositeAndMissingGenotypes=oppositeAndMissingGenotypes)

  # calling function
  test <- slidingWindowCpp(data, gaps, windowSize, step=step, maxGap=maxGap,
                           ROHet=TRUE, maxOppositeGenotype, maxMiss)

  # testing values
  expect_equal(expected, test, info="testing ROHet")

  # setting data for RoHom
  data <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, NA, NA, 0, 1, 0, NA)

  # calling function
  test <- slidingWindowCpp(data, gaps, windowSize, step=step, maxGap=maxGap,
                           ROHet=FALSE, maxOppositeGenotype, maxMiss)

  # testing values
  expect_equal(expected, test, info="testing ROHom")

})

test_that("Testing snpInRun", {
  # setting variables
  windowSize <- 3
  threshold <- 0.5

  # defining slidingWindow result
  windowStatus <- c(F, F, T, T, F, F, T, T, F, F, T, T, T)
  oppositeAndMissingGenotypes <- c("0", "0", "0", "9", "9", "0", "9")
  names(oppositeAndMissingGenotypes) <- c(1, 2, 3, 10, 11, 13, 15)
  res <- list(windowStatus=windowStatus,
              oppositeAndMissingGenotypes=oppositeAndMissingGenotypes)

  # setting expected value
  expected <- c(F, F, F, T, T, F, F, T, T, F, F, T, T, T, T)

  # calling function
  test <- snpInRun(res$windowStatus, windowSize, threshold)

  # testing values
  expect_equal(expected, test)
})

test_that("Testing snpInRunCpp", {
  # setting variables
  windowSize <- 3
  threshold <- 0.5

  # defining slidingWindow result
  windowStatus <- c(F, F, T, T, F, F, T, T, F, F, T, T, T)
  oppositeAndMissingGenotypes <- c("0", "0", "0", "9", "9", "0", "9")
  names(oppositeAndMissingGenotypes) <- c(1, 2, 3, 10, 11, 13, 15)
  res <- list(windowStatus=windowStatus,
              oppositeAndMissingGenotypes=oppositeAndMissingGenotypes)

  # setting expected value
  expected <- c(F, F, F, T, T, F, F, T, T, F, F, T, T, T, T)

  # calling function
  test <- snpInRunCpp(res$windowStatus, windowSize, threshold)

  # testing values
  expect_equal(expected, test)
})


test_that("Testing homoZygotTest", {
  # setting parameters
  maxHet <- 1
  maxMiss <- 1
  maxGap <- 10^6
  i <- 175
  x <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  gaps <- c(3721, 3871, 7059, 4486, 7545, 4796, 3043, 9736, 3495, 5051,
            9607, 6555, 11934, 6410, 3415, 1302, 3110, 6609, 3292)
  windowSize <- length(x)

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(181)
  expected <- list("windowStatus"=TRUE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTest(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a homozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # defining expected value
  oppositeAndMissingSNP <- c(9, 9, 0)
  names(oppositeAndMissingSNP) <- c(175, 176, 181)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTest(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 0)
  gaps[10] <- maxGap + gaps[10]

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(181)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTest(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with big gap")

  # test with a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
            615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)

  # defining expected value
  oppositeAndMissingSNP <- rep(0, 19)
  names(oppositeAndMissingSNP) <- seq(i+1, i+19)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTest(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with heterozygot window")

})

test_that("Testing homoZygotTestCpp", {
  # setting parameters
  maxHet <- 1
  maxMiss <- 1
  maxGap <- 10^6
  i <- 175
  x <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  gaps <- c(3721, 3871, 7059, 4486, 7545, 4796, 3043, 9736, 3495, 5051,
            9607, 6555, 11934, 6410, 3415, 1302, 3110, 6609, 3292)
  windowSize <- length(x)

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(181)
  expected <- list("windowStatus"=TRUE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a homozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # defining expected value
  oppositeAndMissingSNP <- c(9, 9, 0)
  names(oppositeAndMissingSNP) <- c(175, 176, 181)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 0)
  gaps[10] <- maxGap + gaps[10]

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(181)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with big gap")

  # test with a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
            615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)

  # defining expected value
  oppositeAndMissingSNP <- rep(0, 19)
  names(oppositeAndMissingSNP) <- seq(i+1, i+19)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with heterozygot window")

})

test_that("Testing heteroZygotTest", {
  # setting parameters
  maxHom <- 1
  maxMiss <- 1
  maxGap <- 10^6
  i <- 150
  x <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
         1, 1, 1, 1, 0, 0, 1, 0, 0, 0)
  gaps <- c(4374, 8744, 5123, 14229, 5344, 690, 8566, 5853, 2369, 3638,
            4848, 600, 2333, 976, 2466, 2269, 5411, 6021, 4367)

  windowSize <- length(x)

  # defining expected value
  oppositeAndMissingSNP <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  names(oppositeAndMissingSNP) <- c(150, 151, 152, 153, 154, 155, 156, 157,
                                    164, 165, 167, 168, 169)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a homozygot window")

  # test a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(150)
  expected <- list("windowStatus"=TRUE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a heterozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # defining expected value
  oppositeAndMissingSNP <- c(9, 9)
  names(oppositeAndMissingSNP) <- c(150, 151)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 1)
  gaps[10] <- maxGap + gaps[10]

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(150)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTest(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with big gap")

})

test_that("Testing heteroZygotTestCpp", {
  # setting parameters
  maxHom <- 1
  maxMiss <- 1
  maxGap <- 10^6
  i <- 150
  x <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
         1, 1, 1, 1, 0, 0, 1, 0, 0, 0)
  gaps <- c(4374, 8744, 5123, 14229, 5344, 690, 8566, 5853, 2369, 3638,
            4848, 600, 2333, 976, 2466, 2269, 5411, 6021, 4367)

  windowSize <- length(x)

  # defining expected value
  oppositeAndMissingSNP <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  names(oppositeAndMissingSNP) <- c(150, 151, 152, 153, 154, 155, 156, 157,
                                    164, 165, 167, 168, 169)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a homozygot window")

  # test a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(150)
  expected <- list("windowStatus"=TRUE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test a heterozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # defining expected value
  oppositeAndMissingSNP <- c(9, 9)
  names(oppositeAndMissingSNP) <- c(150, 151)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 1)
  gaps[10] <- maxGap + gaps[10]

  # defining expected value
  oppositeAndMissingSNP <- c(0)
  names(oppositeAndMissingSNP) <- c(150)
  expected <- list("windowStatus"=FALSE,
                   "oppositeAndMissingSNP"=oppositeAndMissingSNP)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap, i, windowSize)

  # check for identity
  expect_identical(test, expected, info = "test with big gap")
})

test_that("Testing loading pop from ped", {
  # read data. I need to know only number of columns
  conn  <- file(genotype_path, open = "r")
  oneLine <- readLines(conn, n = 1, warn = FALSE)
  genotype.sample <- (strsplit(oneLine, " "))
  genotype.sample <- as.character(genotype.sample[[1]])
  close(conn)

  # read animals properly
  colClasses <- c(
    rep("character", 2),
    rep("NULL", length(genotype.sample)-2)
  )

  # loading reference_pops
  reference_pops <- read.table(genotype_path, sep = " ", header = FALSE, colClasses = colClasses)
  names(reference_pops) <- c("POP","ID")

  # load pops with a Cpp function
  test_pops <- readPOPCpp(genotype_path)

  # testing
  expect_equal(test_pops, reference_pops)

})
