library(testthat)
library(detectRUNS)
context("Testing RUNS")

# get file paths: reference file need to be changed or removed
genotypeFile  <- "test.ped"
mapFile <- "test.map"

test_that("detected ROHet are identical", {
  # testing slinding windows
  test_sliding <- slidingRUNS.run(genotypeFile, mapFile, windowSize=15, threshold=0.1, minSNP=15,
                         ROHet=FALSE, maxOppWindow=1, maxMissWindow=1,  minLengthBps=100000,
                         minDensity=1/10000, maxOppRun=NULL, maxMissRun=NULL)

  # reading rohet reference: this need to be updated
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  reference_rohet <- read.csv2("test.ROHet.sliding.csv", header = T, stringsAsFactors = FALSE, colClasses = colClasses)

  # compare rohet table
  expect_equal(test_sliding, reference_rohet, info = "testing sliding window approach")

  # testing slinding windows
  test_consecutive <- consecutiveRUNS.run(genotypeFile, mapFile, minSNP=15,
                          ROHet=FALSE,  minLengthBps=100000,
                          maxOppRun=1, maxMissRun=1)

  # reading rohet reference: this need to be updated
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  reference_rohet <- read.csv2("test.ROHet.consecutive.csv", header = T, stringsAsFactors = FALSE, colClasses = colClasses)

  # compare rohet table
  expect_equal(test_consecutive, reference_rohet, info = "testing consecutive approach")
})

test_that("Marker differ in size", {
  # read mapFile
  mapFile <- read.delim(mapFile, header = F)
  names(mapFile) <- c("Chrom","SNP","cM","bps")

  # subset mapfile
  mapFile <- mapFile[100, ]

  # write a fake mapfile in a temporary dir
  fake_mapfile = tempfile(fileext = ".map")
  write.table(mapFile, fake_mapfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # test function
  expect_error(slidingRUNS.run(genotypeFile, fake_mapfile), "Number of markers differ")

  # clean up
  file.remove(fake_mapfile)
})

test_that("No file path throws error", {
  # test for errors
  expect_error(slidingRUNS.run("fake_genotype", mapFile), "doesn't exists")
  expect_error(slidingRUNS.run(genotypeFile, "fake_map"), "doesn't exists")
})
