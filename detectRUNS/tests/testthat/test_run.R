library(testthat)
library(detectRUNS)
context("Testing RUNS")

test_that("detected ROHet are identical", {
  # loading dataset
  data("chillingam")

  # calling function
  test_rohet <- RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
                            ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000,
                            minDensity = 1/10)

  # reading rohet reference
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  reference_path <- system.file("extdata", "detected.ROHet.csv", package = "detectRUNS")
  reference_rohet <- read.csv2(reference_path, header = T, stringsAsFactors = FALSE, colClasses = colClasses)

  # compare rohet table
  expect_equal(test_rohet, reference_rohet)
})

test_that("Dealing with pair alleles (ped)", {
  # loading dataset
  data("chillingam")

  # get a ped file
  genotype_path  <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")

  # check warning message
  expect_warning(RUNS.run(genotype_path, chillingham_map), "genotype with 2 alles for marker detected")
})

test_that("Marker differ in size", {
  # loading dataset
  data("chillingam")

  # subset mapfile
  mapFile <- chillingham_map[100, ]

  expect_error(RUNS.run(chillingham_genotype, mapFile), "Number of markers differ")
})

test_that("No file nor dataframe throws error", {
  # loading dataset
  data("chillingam")

  # test for errors
  expect_error(RUNS.run("fake_genotype", chillingham_map), "doesn't exists")
  expect_error(RUNS.run(chillingham_genotype, "fake_map"), "doesn't exists")
})

test_that("Use file path instead of dataframe", {
  # get file paths
  genotype_path  <- system.file("extdata", "subsetChillingham.raw", package = "detectRUNS")
  mapfile_path <- system.file("extdata", "subsetChillingham.map", package = "detectRUNS")

  test_rohet <- RUNS.run(genotype_path, mapfile_path, windowSize = 20, threshold = 0.1, minSNP = 5,
                            ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000,
                            minDensity = 1/10)

  # reading rohet reference
  colClasses <- c(rep("character", 3), rep("numeric", 4)  )
  reference_path <- system.file("extdata", "detected.ROHet.csv", package = "detectRUNS")
  reference_rohet <- read.csv2(reference_path, header = T, stringsAsFactors = FALSE, colClasses = colClasses)

  # compare rohet table
  expect_equal(test_rohet, reference_rohet)

})
