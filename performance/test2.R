
# clean up
rm(list = ls())

# importing libraries
library(detectRUNS)
library(microbenchmark)

genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package = "detectRUNS")

# read runs
runs <- readExternalRuns(runsFile, program = "detectRUNS")

# how many times perform test
times <- 100

test_tableRuns <- microbenchmark(
  R = tableRuns(runs, mapFile = mapFile, genotypeFile = genotypeFile),
  Cpp = tableRunsCpp(runs, mapFile = mapFile, genotypeFile = genotypeFile),
  unit = 'ms',
  times = times
)

boxplot(test_tableRuns)
