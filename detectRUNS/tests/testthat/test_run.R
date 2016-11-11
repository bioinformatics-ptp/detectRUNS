library(detectRUNS)
context("Testing RUNS")

test_that("detected ROHet are identical", {
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

  #  sort dataframe
  # reference_rohet <- reference_rohet[order(reference_rohet$breed, reference_rohet$id), ]

  # change row_names
  # row.names(reference_rohet) <- 1:nrow(reference_rohet)

  # reading test
  test_rohet <- read.csv2("detected.ROHet.csv", header = T)

  # sort and add row names
  # test_rohet <- test_rohet[order(test_rohet$breed, test_rohet$id), ]
  # row.names(test_rohet) <- 1:nrow(test_rohet)

  # compare rohet table
  expect_identical(test_rohet, reference_rohet)

  # cleanup
  file.remove("detected.ROHet.csv")
  file.remove("genotype.raw")
  file.remove("plink.map")
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

  n_individuals <- RUNS.run(chillingham_genotype, chillingham_map, windowSize = 20, threshold = 0.1, minSNP = 5,
                            ROHet = TRUE, maxOppositeGenotype = 1, maxMiss = 1,  minLengthBps = 1000,
                            minDensity = 1/10)

  # testing function
  expect_equal(n_individuals, 3)

  reference_path <- system.file("extdata", "detected.ROHet.csv", package = "detectRUNS")
  reference_rohet <- read.csv2(reference_path, header = T)

  # reading test
  test_rohet <- read.csv2("detected.ROHet.csv", header = T)

  # compare rohet table
  expect_identical(test_rohet, reference_rohet)

  # cleanup
  file.remove("detected.ROHet.csv")
  file.remove("genotype.raw")
  file.remove("plink.map")

})
