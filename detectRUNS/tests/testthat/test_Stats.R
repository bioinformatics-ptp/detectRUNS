library(testthat)
library(detectRUNS)
context("Testing Stats")

# get file paths
genotypeFile  <- "test.ped"
mapFile <- "test.map"
rawFile <- "test.raw"
runsFile <- "test.ROHet.consecutive.csv"

test_that("Test tableRuns", {
  colClasses <- c("factor", "character", "character", "character", "numeric",
                  "integer", "integer", "numeric")
  reference <- read.csv2("test.tableRuns.csv", colClasses = colClasses)
  levels(reference$Group) <- c("Jacobs", "Navajo-Churro")
  
  runsDF <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')
  test <- tableRuns(runs = runsDF, genotypeFile = genotypeFile, mapFile = mapFile, threshold = 0.5)
  
  expect_equal(reference, test)
})
