library(detectRUNS)
context("Testing RUNS")

test_that("n_individuals is equal to 3", {
  # loading dataset
  data("chillingam")

  # calling function
  n_individuals <- RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
                            ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000,
                            minDensity = 1/10)

  # testing function
  expect_equal(n_individuals, 3)

  # reading rohet, reference and test
  reference_path <- system.file("extdata", "detected.ROHet.csv", package = "detectRUNS")
  reference_rohet <- read.csv2(reference_path, header = T)
  test_rohet <- read.csv2("detected.ROHet.csv", header = T)

  # compare rohet table
  expect_identical(test_rohet, reference_rohet)

  # cleanup
  file.remove("detected.ROHet.csv")
  file.remove("genotype.raw")
  file.remove("plink.map")
})
