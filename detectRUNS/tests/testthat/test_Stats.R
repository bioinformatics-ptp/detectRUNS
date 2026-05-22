library(testthat)
library(detectRUNS)
context("Testing Stats")

# get file paths
genotypeFile <- "test.ped"
mapFile <- "test.map"
rawFile <- "test.raw"
runsFile <- "test.ROHet.consecutive.csv"

test_that("Test tableRuns", {
  colClasses <- c(
    "character", "character", "character", "character", "numeric",
    "integer", "integer"
  )
  reference <- read.csv2("test.tableRuns.csv", colClasses = colClasses)

  runsDF <- readExternalRuns(inputFile = runsFile, program = "detectRUNS")
  test <- tableRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile, threshold = 0.5)
  test$Group <- as.character(test$Group)

  expect_equal(reference, test)

  expect_error(
    tableRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile, threshold = 50),
    "Threshold must be between 0 and 1"
  )
})
