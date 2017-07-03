library(testthat)
library(detectRUNS)
context("Testing RUNS")

# TODO:
# Data for test should be small (don't test on real data)

# get file paths: reference file need to be changed or removed
genotype_path  <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")

test_that("detected ROHet are identical", {
  # testing slinding windows
  test_rohet <- RUNS.run(genotype_path, mapfile_path, windowSize=20, threshold=0.1, minSNP=5,
                         ROHet=TRUE, maxOppositeGenotype=1, maxMiss=1,  minLengthBps=1000,
                         minDensity=1/10, maxOppRun=NULL, maxMissRun=NULL, method='slidingWindow')

  # reading rohet reference: this need to be updated
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  reference_path <- system.file("extdata", "detected.ROHet.csv", package = "detectRUNS")
  reference_rohet <- read.csv2(reference_path, header = T, stringsAsFactors = FALSE, colClasses = colClasses)

  # compare rohet table
  expect_equal(test_rohet, reference_rohet)
})

test_that("Marker differ in size", {
  # read mapFile
  mapFile <- read.delim(mapfile_path, header = F)
  names(mapFile) <- c("Chrom","SNP","cM","bps")

  # subset mapfile
  mapFile <- mapFile[100, ]

  # write a fake mapfile (TODO: in temporary dir?)
  fake_mapfile = "fake_mapfile.map"
  write.table(mapFile, fake_mapfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # test function
  expect_error(RUNS.run(genotype_path, fake_mapfile), "Number of markers differ")

  # clean up
  file.remove(fake_mapfile)
})

test_that("No file path throws error", {
  # test for errors
  expect_error(RUNS.run("fake_genotype", mapfile_path), "doesn't exists")
  expect_error(RUNS.run(genotype_path, "fake_map"), "doesn't exists")
})
