library(testthat)
library(detectRUNS)
context("Testing plots")

runsFile <- "test.ROHet.sliding.csv"
genotypeFile  <- "test.ped"
mapFile <- "test.map"

test_that("Test plot_manhattanRuns", {
  # loading data from CSV
  runs <- readExternalRuns(inputFile = runsFile, program = 'detectRUNS')

  # plotting data
  plot_manhattanRuns(runs, genotypeFile, mapFile, savePlots = F)
})