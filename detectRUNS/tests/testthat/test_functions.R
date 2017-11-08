library(testthat)
library(detectRUNS)
context("Testing functions")

# get file paths
genotypeFile  <- "test.ped"
mapFile <- "test.map"
rawFile <- "test.raw"

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
  # read genotype
  genotype <- read.table(genotypeFile, sep = " ", header = FALSE,
                         stringsAsFactors = FALSE, colClasses = "character")
  genotype <- genotype[ , -c(1:6)]

  # read raw data
  raw <- read.table(rawFile, sep=" ", header = TRUE)
  raw <-  raw[ , -c(1:6)]

  t1 <- apply(genotype, 1, pedConvertCpp)
  t2 <- apply(raw, 1, genoConvertCpp)

  # testing conversion
  expect_identical(t1, t2)
})

test_that("Testing slidingWindow", {
  # setting parameters
  windowSize <- 3
  step <- 1
  maxGap <- 1000
  maxOppWindow <- 1
  maxMissWindow <- 1

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
                        ROHet=TRUE, maxOppWindow, maxMissWindow)

  # testing values
  expect_equal(expected, test, info="testing ROHet")

  # setting data for RoHom
  data <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, NA, NA, 0, 1, 0, NA)

  # calling function
  test <- slidingWindow(data, gaps, windowSize, step=step, maxGap=maxGap,
                        ROHet=FALSE, maxOppWindow, maxMissWindow)

  # testing values
  expect_equal(expected, test, info="testing ROHom")

})

test_that("Testing findOppositeAndMissing", {
  # define values
  data <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA)

  # testing ROHet
  reference <- c("0", "0", "0", "9", "9", "0", "9")
  names(reference) <- c("1", "2", "3", "10", "11", "13", "15")

  test <- findOppositeAndMissing(data, ROHet=TRUE)
  expect_equal(reference, test, info="testing ROHet")

  # testing ROHom
  reference <- c("0", "0", "0", "0", "0", "0", "9", "9", "0", "0", "9")
  names(reference) <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "15")

  test <- findOppositeAndMissing(data, ROHet=FALSE)
  expect_equal(reference, test, info="testing ROHom")

})

test_that("Testing slidingWindowCpp", {
  # setting parameters
  windowSize <- 3
  step <- 1
  maxGap <- 1000
  maxOppWindow <- 1
  maxMissWindow <- 1

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
                           ROHet=TRUE, maxOppWindow, maxMissWindow)

  # testing values
  expect_equal(expected, test, info="testing ROHet")

  # setting data for RoHom
  data <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, NA, NA, 0, 1, 0, NA)

  # calling function
  test <- slidingWindowCpp(data, gaps, windowSize, step=step, maxGap=maxGap,
                           ROHet=FALSE, maxOppWindow, maxMissWindow)

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

test_that("Testing createRUNdf", {
  ####################
  # setting variables

  data <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0,
            1, 0, 0, 1, NA, 1, 0, 0, 1, 0, 1, NA)
  minSNP <- 2
  minLengthBps <- 100
  minDensity <- 0.1
  maxGap <- 1000
  windowSize <- 3
  step <- 1
  maxOppWindow <- 1
  maxMissWindow <- 1
  threshold <- 0.5

  # defining gaps
  gaps <- rep(100, length(data)-1)

  # update a gap: setting value higher than threshold
  gaps[5] <- gaps[5] + maxGap

  # update gap to increase run size
  gaps[7] <- gaps[7] + 100

  # defining slidingWindow result
  res <- slidingWindow(data, gaps, windowSize, step=step, maxGap=maxGap,
                       ROHet=TRUE, maxOppWindow, maxMissWindow)

  # define snpInRun results
  snpRun <- snpInRun(res$windowStatus, windowSize, threshold)

  # snpRun is: c(F, F, T, T, F, F, T, T, T, F, F, T, T,
  #              F, F, F, T, T, T, F, F, F, T, T, T)

  # defining a fake mapfile
  mapFile <- data.frame(Chrom=rep(1, length(data)),
                     SNP=seq(1, length(data)),
                     cM=rep(0, length(data)),
                     bps=rep(0, length(data)))

  # updating positions
  for (i in 1:length(gaps)) {
    mapFile[i+1, "bps"] <- mapFile[i, "bps"] + gaps[i]
  }

  mapFile[1, "bps"] <- 1

  # defining expected dataframe with all datas
  from=as.numeric(c(200, 1600, 2200, 2700, 3300))
  to=as.numeric(c(300, 1900, 2300, 2900, 3500))
  nSNP=as.integer(c(2, 3, 2, 3, 3))
  chrom=as.character(rep(1, 5))
  lengthBps=as.numeric(c(100, 300, 100, 200, 200))

  reference <- data.frame(from=from,
                          to=to,
                          nSNP=nSNP,
                          chrom=chrom,
                          lengthBps=lengthBps,
                          stringsAsFactors = FALSE)

  # reference with nOpp and nMiss columns
  #   from   to nSNP chrom lengthBps nOpp nMiss
  # 1  200  300    2     1       100    0     0
  # 2 1600 1900    3     1       300    0     0
  # 3 2200 2300    2     1       100    1     0
  # 4 2700 2900    3     1       200    0     1
  # 5 3300 3500    3     1       200    1     1

  #####################################
  # calling function without filtering

  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes)

  # testing values
  expect_equal(reference, test, info = "whithout filtering")

  ##############################
  # testing for minlength = 300

  expected <- reference[2, ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      300,
                      minDensity,
                      res$oppositeAndMissingGenotypes)

  # testing values
  expect_equal(expected, test, info = "filtering by minlength")

  ##################################
  # testing with no missing in RUNs

  expected <- reference[1:3, ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxMissRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noMissing")

  ###################################
  # testing with no opposite in RUNs

  expected <- reference[c(1, 2, 4), ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxOppRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noOpposite")

  #################################################
  # testing with no opposite and no missing in RUNs

  expected <- reference[c(1:2), ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxOppRun=0,
                      maxMissRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noOpposite and noMissing")

  ###################################
  # testing with no opposite in RUNs and no missing and with minSNP

  expected <- reference[c(0), ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapFile,
                      10, # minSNP
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxOppRun=0,
                      maxMissRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noOpposite and noMissing and minSNP = 10")

  ##############################
  # testing on different chroms

  mapFile[9:nrow(mapFile), "Chrom"] <- 2
  mapFile[9:nrow(mapFile), "bps"] <- mapFile[9:nrow(mapFile), "bps"] - 1800

  # re calculate gaps
  gaps <- diff(mapFile$bps)

  # defining expected dataframe with all datas (there is a min length > 100)
  from=as.numeric(c(200, 1600, 400, 900, 1500))
  to=as.numeric(c(300, 1800, 500, 1100, 1700))
  nSNP=as.integer(c(2, 2, 2, 3, 3))
  chrom=as.character(c(rep(1, 2), rep(2, 3)))
  lengthBps=as.numeric(c(100, 200, 100, 200, 200))

  reference <- data.frame(from=from,
                          to=to,
                          nSNP=nSNP,
                          chrom=chrom,
                          lengthBps=lengthBps,
                          stringsAsFactors = FALSE)

  # reference with nOpp and nMiss columns
  # from   to nSNP chrom lengthBps
  # 1  200  300    2     1       100
  # 2  400  500    2     2       100
  # 3  900 1100    3     2       200
  # 4 1500 1700    3     2       200

  # calling function
  test <- createRUNdf(snpRun,
                      mapFile,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes)

  # testing values
  expect_equal(reference, test, info = "testing on different chroms")

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

  # homo array
  x <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  gaps <- c(3721, 3871, 7059, 4486, 7545, 4796, 3043, 9736, 3495, 5051,
            9607, 6555, 11934, 6410, 3415, 1302, 3110, 6609, 3292)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap)

  # check for identity
  expect_true(test, info = "test a homozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 0)
  gaps[10] <- maxGap + gaps[10]

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test with big gap")

  # test with a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
            615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)

  # calling function
  test <- homoZygotTestCpp(x, gaps, maxHet, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test with heterozygot window")

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

  # homo array
  x <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
         1, 1, 1, 1, 0, 0, 1, 0, 0, 0)
  gaps <- c(4374, 8744, 5123, 14229, 5344, 690, 8566, 5853, 2369, 3638,
            4848, 600, 2333, 976, 2466, 2269, 5411, 6021, 4367)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test a homozygot window")

  # test a hetero array
  x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)

  # check for identity
  expect_true(test, info = "test a heterozygot window")

  # insert two missing values (> naxMIss)
  x[1:2] <- c(NA, NA)

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test with missing values")

  # revert, and change maxGap
  x[1:2] <- c(0, 1)
  gaps[10] <- maxGap + gaps[10]

  # calling function
  test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)

  # check for identity
  expect_false(test, info = "test with big gap")
})

test_that("Testing loading pop from ped", {
  # read data. I need to know only number of columns
  conn  <- file(genotypeFile, open = "r")
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
  reference_pops <- read.table(genotypeFile, sep = " ", header = FALSE, colClasses = colClasses)
  names(reference_pops) <- c("POP","ID")

  # load pops with a Cpp function
  test_pops <- readPOPCpp(genotypeFile)

  # testing
  expect_equal(test_pops, reference_pops)

})

test_that("Testing consecutiveRuns", {
  # setting parameters
  maxGap <- 1000
  maxOppRun <- 1
  maxMissRun <- 1
  minLengthBps <- 100
  ROHet <- TRUE
  minSNP <- 1

  # setting values for RoHet
  indGeno <- c(0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA, 0, 1, 1)
  ind <- list(FID="foo", IID="bar")
  gaps <- rep(100, length(indGeno)-1)

  # update a gap: setting value higher than threshold
  gaps[4] <- gaps[4] + maxGap

  # defining a fake mapfile
  mapFile <- data.frame(Chrom=c(rep(1, 9), rep(2, 7)),
                        SNP=seq(1, length(indGeno)),
                        cM=rep(0, length(indGeno)),
                        bps=rep(0, length(indGeno)))

  # updating positions
  for (i in 1:length(gaps)) {
    mapFile[i+1, "bps"] <- mapFile[i, "bps"] + gaps[i]
  }

  mapFile[1, "bps"] <- 1

  # defining expected dataframe with all datas
  from <- as.numeric(c(100, 1400, 1900, 2400))
  to <- as.numeric(c(300, 1700, 2200, 2500))
  nSNP <- as.integer(c(3, 4, 4, 2))
  chrom <- as.character(c(rep(1, 2), rep(2, 2)))
  lengthBps <- as.numeric(c(200, 300, 300, 100))
  group <- as.character(rep(ind$FID), 4)
  id <- as.character(rep(ind$IID), 4)

  reference <- data.frame(group=group,
                          id=id,
                          chrom=chrom,
                          nSNP=nSNP,
                          from=from,
                          to=to,
                          lengthBps=lengthBps,
                          stringsAsFactors = FALSE)

  # reference is
  #   group  id chrom nSNP from   to lengthBps
  # 1   foo bar     1    3  100  300       200
  # 2   foo bar     1    4 1400 1700       300
  # 3   foo bar     2    4 1900 2200       300
  # 4   foo bar     2    2 2400 2500       100

  # calling consecutiveRuns
  test <- consecutiveRuns(indGeno, ind, mapFile, ROHet=ROHet, minSNP=minSNP,
                          maxOppositeGenotype=maxOppRun, maxMiss=maxMissRun,
                          minLengthBps=minLengthBps,maxGap=maxGap)

  # testing functions
  expect_equal(test, reference)

})

test_that("Testing consecutiveRunsCpp", {
  # setting parameters
  maxGap <- 1000
  maxOppRun <- 1
  maxMissRun <- 1
  minLengthBps <- 100
  ROHet <- TRUE
  minSNP <- 1

  # setting values for RoHet
  indGeno <- c(0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA, 0, 1, 1)
  ind <- list(FID="foo", IID="bar")
  gaps <- rep(100, length(indGeno)-1)

  # update a gap: setting value higher than threshold
  gaps[4] <- gaps[4] + maxGap

  # defining a fake mapfile
  mapFile <- data.frame(Chrom=c(rep(1, 9), rep(2, 7)),
                        SNP=seq(1, length(indGeno)),
                        cM=rep(0, length(indGeno)),
                        bps=rep(0, length(indGeno)))

  # updating positions
  for (i in 1:length(gaps)) {
    mapFile[i+1, "bps"] <- mapFile[i, "bps"] + gaps[i]
  }

  mapFile[1, "bps"] <- 1

  # defining expected dataframe with all datas
  from <- as.numeric(c(100, 1400, 1900, 2400))
  to <- as.numeric(c(300, 1700, 2200, 2500))
  nSNP <- as.integer(c(3, 4, 4, 2))
  chrom <- as.character(c(rep(1, 2), rep(2, 2)))
  lengthBps <- as.numeric(c(200, 300, 300, 100))
  group <- as.character(rep(ind$FID), 4)
  id <- as.character(rep(ind$IID), 4)

  reference <- data.frame(group=group,
                          id=id,
                          chrom=chrom,
                          nSNP=nSNP,
                          from=from,
                          to=to,
                          lengthBps=lengthBps,
                          stringsAsFactors = FALSE)

  # reference is
  #   group  id chrom nSNP from   to lengthBps
  # 1   foo bar     1    3  100  300       200
  # 2   foo bar     1    4 1400 1700       300
  # 3   foo bar     2    4 1900 2200       300
  # 4   foo bar     2    2 2400 2500       100

  # calling consecutiveRuns
  test <- consecutiveRunsCpp(indGeno, ind, mapFile, ROHet=ROHet, minSNP=minSNP,
                             maxOppositeGenotype=maxOppRun, maxMiss=maxMissRun,
                             minLengthBps=minLengthBps,maxGap=maxGap)

  # testing functions
  expect_equal(test, reference)

})

test_that("Testing snpInsideRuns", {
  # defining mapChrom
  mappa <- data.table::fread(mapFile, header = FALSE)
  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL
  chrom <- "24"

  # subsetting mapChrom
  set.seed(42)
  mapChrom <- mappa[mappa$CHR==chrom, ]
  mapChrom <- mapChrom[sort(sample(nrow(mapChrom), 10)), ]

  # loading pre-calculated data
  runsFile <- "test.ROHet.sliding.csv"
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE,
                    colClasses = colClasses)

  # fix column names and define runsChrom
  names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")
  runsChrom <- runs[runs$CHROMOSOME==chrom, ]

  # get snps inside runs
  reference <- snpInsideRuns(runsChrom, mapChrom, genotypeFile)
  test <- snpInsideRunsCpp(runsChrom, mapChrom, genotypeFile)

  # testing functions
  expect_equivalent(test, reference)

})
