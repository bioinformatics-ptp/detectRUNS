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
  maxOppositeGenotype <- 1
  maxMiss <- 1
  threshold <- 0.5

  # defining gaps
  gaps <- rep(100, length(data)-1)

  # update a gap: setting value higher than threshold
  gaps[5] <- gaps[5] + maxGap

  # update gap to increase run size
  gaps[7] <- gaps[7] + 100

  # defining slidingWindow result
  res <- slidingWindow(data, gaps, windowSize, step=step, maxGap=maxGap,
                       ROHet=TRUE, maxOppositeGenotype, maxMiss)

  # define snpInRun results
  snpRun <- snpInRun(res$windowStatus, windowSize, threshold)

  # snpRun is: c(F, F, T, T, F, F, T, T, T, F, F, T, T,
  #              F, F, F, T, T, T, F, F, F, T, T, T)

  # defining a fake mapfile
  mapa <- data.frame(Chrom=rep(1, length(data)),
                     SNP=seq(1, length(data)),
                     cM=rep(0, length(data)),
                     bps=rep(0, length(data)))

  # updating positions
  for (i in 1:length(gaps)) {
    mapa[i+1, "bps"] <- mapa[i, "bps"] + gaps[i]
  }

  mapa[1, "bps"] <- 1

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
                      mapa,
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
                      mapa,
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
                      mapa,
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
                      mapa,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxOppRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noOpposite")

  # testing with no opposite and no missing in RUNs
  expected <- reference[c(1:2), ]
  row.names(expected) <- NULL

  test <- createRUNdf(snpRun,
                      mapa,
                      minSNP,
                      minLengthBps,
                      minDensity,
                      res$oppositeAndMissingGenotypes,
                      maxOppRun=0,
                      maxMissRun=0)

  # testing values
  expect_equal(expected, test, info = "filtering by noOpposite and noMissing")

  ##############################
  # testing on different chroms

  mapa[9:nrow(mapa), "Chrom"] <- 2

  # defining expected dataframe with all datas
  from=as.numeric(c(200, 1600, 1800, 2200, 2700, 3300))
  to=as.numeric(c(300, 1600, 1900, 2300, 2900, 3500))
  nSNP=as.integer(c(2, 1, 2, 2, 3, 3))
  chrom=as.character(c(rep(1, 2), rep(2, 4)))
  lengthBps=as.numeric(c(100, 0, 200, 100, 200, 200))

  reference <- data.frame(from=from,
                          to=to,
                          nSNP=nSNP,
                          chrom=chrom,
                          lengthBps=lengthBps,
                          stringsAsFactors = FALSE)

  # reference with nOpp and nMiss columns
  #   from   to nSNP chrom lengthBps nOpp nMiss
  # 1  200  300    2     1       100    0     0
  # 2 1600 1600    1     1         0    0     0
  # 3 1800 1900    2     2       200    0     0
  # 4 2200 2300    2     2       100    1     0
  # 5 2700 2900    3     2       200    0     1
  # 6 3300 3500    3     2       200    1     1

  # calling function
  test <- createRUNdf(snpRun,
                      mapa,
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
