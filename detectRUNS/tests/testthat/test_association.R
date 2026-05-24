library(testthat)
library(detectRUNS)
context("Testing runsAssociation")

# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------
#
# 10 animals (ind1..ind10), 3 ROH regions:
#   region A  chr1_1000_5000  : present in ind1..ind6  (6/10 = 60%)
#   region B  chr1_6000_9000  : present in ind1..ind2  (2/10 = 20%)
#   region C  chr2_1000_3000  : present in ind1        (1/10 = 10%)
#
# Phenotype:
#   ind1..ind6  -> 2.0   (carriers of region A)
#   ind7..ind10 -> 1.0   (non-carriers)
#
# With minFreq = 0.30: region A (60%) and region B (20%) pass;
#                      region C (10%) is filtered out.
# With minFreq = 0.50: only region A (60%) passes.
#
# For region A: perfect separation -> beta = 1.0, p << 0.05
# For region B: 2 carriers with pheno 2.0, 8 non-carriers with mean 1.5
#               beta = 2.0 - 1.5 = 0.5, p will be >0.05 (only 2 carriers)

make_runs <- function() {
  data.frame(
    group     = "PopA",
    id        = c(paste0("ind", 1:6),   # region A carriers
                  paste0("ind", 1:2),   # region B carriers
                  "ind1"),              # region C carrier
    chrom     = c(rep(1, 6), rep(1, 2), 2),
    nSNP      = 20L,
    from      = c(rep(1000L, 6), rep(6000L, 2), 1000L),
    to        = c(rep(5000L, 6), rep(9000L, 2), 3000L),
    lengthBps = c(rep(4000L, 6), rep(3000L, 2), 2000L),
    stringsAsFactors = FALSE
  )
}

make_pheno <- function() {
  data.frame(
    id    = paste0("ind", 1:10),
    pheno = c(rep(2.0, 6), rep(1.0, 4)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------

test_that("returns data.frame with correct columns", {
  res <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.50)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("region", "n_animals", "beta", "pvalue",
                      "pvalue_bonferroni", "pvalue_fdr"))
})

test_that("region A: known beta = 1 and p-value significant", {
  res <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.50)
  row_A <- res[res$region == "chr1_1000_5000", ]
  expect_equal(nrow(row_A), 1L)
  expect_equal(row_A$n_animals, 6L)
  expect_equal(row_A$beta, 1.0, tolerance = 1e-10)
  expect_lt(row_A$pvalue, 0.05)
  # with 1 region tested, Bonferroni = BH = raw p-value
  expect_equal(row_A$pvalue_bonferroni, row_A$pvalue, tolerance = 1e-10)
  expect_equal(row_A$pvalue_fdr,        row_A$pvalue, tolerance = 1e-10)
})

test_that("minFreq filters correctly", {
  # minFreq = 0.50: only region A (60%) survives
  res_strict <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.50)
  expect_true(all(res_strict$region == "chr1_1000_5000"))

  # minFreq = 0.10: all three regions survive
  res_loose <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.10)
  expect_equal(nrow(res_loose), 3L)
  expect_setequal(res_loose$region,
                  c("chr1_1000_5000", "chr1_6000_9000", "chr2_1000_3000"))
})

test_that("accepts a plain data.frame and ROH S3 object identically", {
  plain_runs <- make_runs()
  res_plain  <- runsAssociation(plain_runs, make_pheno(), minFreq = 0.50)

  # Wrap in a minimal ROH S3 object (same structure as package internals)
  roh_obj <- structure(
    list(runs = plain_runs),
    class = "RUNS"
  )
  res_s3 <- runsAssociation(roh_obj, make_pheno(), minFreq = 0.50)

  expect_equal(res_plain, res_s3)
})

test_that("custom phenoCol name works", {
  pheno2 <- make_pheno()
  names(pheno2)[names(pheno2) == "pheno"] <- "milk_yield"
  res <- runsAssociation(make_runs(), pheno2, minFreq = 0.50, phenoCol = "milk_yield")
  expect_equal(nrow(res), 1L)
  expect_equal(res$beta, 1.0, tolerance = 1e-10)
})

test_that("no ROH passes minFreq returns empty data.frame", {
  res <- suppressMessages(
    runsAssociation(make_runs(), make_pheno(), minFreq = 0.99)
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
  expect_named(res, c("region", "n_animals", "beta", "pvalue",
                      "pvalue_bonferroni", "pvalue_fdr"))
})

test_that("region with no variation returns NA for beta and pvalue", {
  # Give every animal in pheno the same ROH (ind1..ind10 all carriers)
  # -> after merge, all individuals are carriers -> no variation -> NA
  runs_all <- data.frame(
    group     = "PopA",
    id        = paste0("ind", 1:10),
    chrom     = 1,
    nSNP      = 20L,
    from      = 1000L,
    to        = 5000L,
    lengthBps = 4000L,
    stringsAsFactors = FALSE
  )
  res <- runsAssociation(runs_all, make_pheno(), minFreq = 0.10)
  expect_equal(nrow(res), 1L)
  expect_true(is.na(res$beta))
  expect_true(is.na(res$pvalue))
})

test_that("error on duplicate 'id' values in pheno", {
  bad_pheno <- make_pheno()
  bad_pheno <- rbind(bad_pheno, bad_pheno[1, ])
  expect_error(runsAssociation(make_runs(), bad_pheno),
               "duplicate 'id' values")
})

test_that("NA phenotype values are excluded from regression but not from frequency", {
  pheno_na <- make_pheno()
  pheno_na$pheno[7:8] <- NA   # two non-carriers get NA phenotype

  res <- runsAssociation(make_runs(), pheno_na, minFreq = 0.50)
  row_A <- res[res$region == "chr1_1000_5000", ]

  # n_animals counts only phenotyped carriers (ind1..ind6, all have pheno)
  expect_equal(row_A$n_animals, 6L)
  # beta should still be 1: carriers mean=2, non-carriers with pheno (ind7..ind10 - NA) mean=1
  expect_equal(row_A$beta, 1.0, tolerance = 1e-10)
  expect_false(is.na(row_A$pvalue))
})

test_that("error on missing 'id' column in pheno", {
  bad_pheno <- data.frame(animal = "ind1", pheno = 1.0)
  expect_error(runsAssociation(make_runs(), bad_pheno),
               "must have a column named 'id'")
})

test_that("error on missing phenotype column", {
  expect_error(
    runsAssociation(make_runs(), make_pheno(), phenoCol = "nonexistent"),
    "Column 'nonexistent' not found"
  )
})

test_that("error on invalid minFreq", {
  expect_error(runsAssociation(make_runs(), make_pheno(), minFreq = 0),
               "'minFreq' must be a single numeric value in \\(0, 1\\]")
  expect_error(runsAssociation(make_runs(), make_pheno(), minFreq = 1.5),
               "'minFreq' must be a single numeric value in \\(0, 1\\]")
  expect_error(runsAssociation(make_runs(), make_pheno(), minFreq = "0.1"),
               "'minFreq' must be a single numeric value in \\(0, 1\\]")
})

test_that("error on missing required columns in runs", {
  bad_runs <- make_runs()
  bad_runs$from <- NULL
  expect_error(runsAssociation(bad_runs, make_pheno()),
               "missing required columns: from")
})

test_that("no overlap between pheno and runs returns empty data.frame", {
  # pheno individual is not in runs -> frequency = 0 for all regions -> empty result
  pheno_other <- data.frame(id = "unknown_animal", pheno = 1.0,
                            stringsAsFactors = FALSE)
  res <- suppressMessages(runsAssociation(make_runs(), pheno_other))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
  expect_named(res, c("region", "n_animals", "beta", "pvalue",
                      "pvalue_bonferroni", "pvalue_fdr"))
})

test_that("Bonferroni and FDR adjusted p-values are correct", {
  res <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.10)
  n   <- nrow(res)

  # Bonferroni: raw p * number of tests, capped at 1
  expected_bonf <- pmin(res$pvalue * n, 1)
  expect_equal(res$pvalue_bonferroni, expected_bonf, tolerance = 1e-10)

  # BH-FDR: must equal p.adjust(pvalue, "BH")
  expected_fdr <- p.adjust(res$pvalue, method = "BH")
  expect_equal(res$pvalue_fdr, expected_fdr, tolerance = 1e-10)

  # adjusted p-values must be >= raw p-values (correction never makes things more significant)
  expect_true(all(res$pvalue_bonferroni >= res$pvalue, na.rm = TRUE))
  expect_true(all(res$pvalue_fdr        >= res$pvalue, na.rm = TRUE))
})

test_that("n_animals counts correctly for each region", {
  res <- runsAssociation(make_runs(), make_pheno(), minFreq = 0.10)
  n <- setNames(res$n_animals, res$region)
  expect_equal(n[["chr1_1000_5000"]], 6L)
  expect_equal(n[["chr1_6000_9000"]], 2L)
  expect_equal(n[["chr2_1000_3000"]], 1L)
})
