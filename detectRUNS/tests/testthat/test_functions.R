library(testthat)
library(detectRUNS)
context("Testing functions")

# get file paths
genotype_path  <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")
raw_path <- system.file("extdata", "subsetChillingham.raw", package = "detectRUNS")

# inporting data once
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
  genotype_raw <- chillingham_raw

  # remove unnecessary fields from the .ped file
  genotype <- genotype[ ,-c(1:6)]
  genotype_raw <- genotype_raw[ ,-c(1:6)]

  # converting ped file into raw matrix
  genotype <- apply(genotype, 1, pedConvertCpp)
  genotype <- t(genotype)

  # get map data
  mapFile <- chillingham_map

  # setting colnames
  names(mapFile) <- c("Chrom","SNP","cM","bps")

  # calculating gaps
  gaps <- diff(mapFile$bps)

  # define an internal function
  is_run <- function(x, x_raw) {
    # call R function
    y <- slidingWindow(x_raw, gaps, windowSize, step=1, maxGap=maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # call cppFunction
    test <- slidingWindowCpp(x, gaps, windowSize, step=1, maxGap=maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # testing function
    expect_identical(test, y)

    # call R function for RoHom
    y <- slidingWindow(x_raw, gaps, windowSize, step=1, maxGap=maxGap, ROHet=FALSE, maxOppositeGenotype, maxMiss);

    # call cppFunction
    test <- slidingWindowCpp(x, gaps, windowSize, step=1, maxGap=maxGap, ROHet=FALSE, maxOppositeGenotype, maxMiss);

    # testing function
    expect_identical(test, y)

    # call R function and change steps
    y <- slidingWindow(x_raw, gaps, windowSize, step=5, maxGap=maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # call cppFunction
    test <- slidingWindowCpp(x, gaps, windowSize, step=5, maxGap=maxGap, ROHet=ROHet, maxOppositeGenotype, maxMiss);

    # testing function
    expect_identical(test, y)
  }

  # test genotype with ped file for Cpp function and raw file for R functions
  for (i in 1:nrow(genotype)) {
    is_run( genotype[i, ], as.integer(genotype_raw[i, ]) )
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
