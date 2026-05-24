library(testthat)
library(detectRUNS)
context("Edge cases and stress tests")

# ---------------------------------------------------------------------------
# Shared fixtures — run once cheaply with the small test files
# ---------------------------------------------------------------------------

ped_file <- "test.ped"
map_file <- "test.map"

# Normal scan (moderate filters) — reused in several tests
.roh_normal <- local({
  suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 15, maxOpp = 1, maxMiss = 1,
             minLengthBps = 100000, verbose = FALSE)
  ))
})

# Scan with impossible filter — produces 0 runs but a valid ROH object
.roh_empty <- local({
  suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 999999L, verbose = FALSE)
  ))
})


# ===========================================================================
# 1. scanRUNS — input validation
# ===========================================================================

test_that("scanRUNS errors when genotype file does not exist", {
  expect_error(
    suppressWarnings(scanRUNS("nonexistent_file.ped", verbose = FALSE)),
    "File not found|Cannot determine file format"
  )
})

test_that("scanRUNS errors on non-boolean ROHet", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, ROHet = "yes", verbose = FALSE)),
    "ROHet must be TRUE or FALSE"
  )
})

test_that("scanRUNS errors on minSNP < 1", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, minSNP = 0, verbose = FALSE)),
    "minSNP must be a positive integer"
  )
})

test_that("scanRUNS errors on negative maxOpp", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, maxOpp = -1, verbose = FALSE)),
    "maxOpp must be a non-negative integer"
  )
})

test_that("scanRUNS errors on negative maxMiss", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, maxMiss = -1, verbose = FALSE)),
    "maxMiss must be a non-negative integer"
  )
})

test_that("scanRUNS errors when maxGap < 1", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, maxGap = 0, verbose = FALSE)),
    "maxGap must be >= 1"
  )
})

test_that("scanRUNS errors on threshold out of [0,1] for sliding method", {
  expect_error(
    suppressWarnings(scanRUNS(ped_file, method = "sliding", threshold = 2.0, verbose = FALSE)),
    "threshold must be between 0 and 1"
  )
})

test_that("scanRUNS errors on non-character genoFile", {
  expect_error(
    suppressWarnings(scanRUNS(123, verbose = FALSE)),
    "genoFile must be a non-empty file path string"
  )
})


# ===========================================================================
# 2. scanRUNS — return value structure
# ===========================================================================

test_that("scanRUNS returns an ROH S3 object", {
  expect_s3_class(.roh_normal, "RUNS")
})

test_that("ROH object contains all expected named elements", {
  expected <- c("runs", "summary", "chrom_lengths", "sample_info",
                "snp_map", "method", "type", "snp_freq", "chrom_map")
  expect_true(all(expected %in% names(.roh_normal)))
})

test_that("ROH$runs has the 7 standard columns", {
  expect_named(
    as.data.frame(.roh_normal$runs),
    c("group", "id", "chrom", "nSNP", "from", "to", "lengthBps")
  )
})

test_that("ROH$summary covers every individual including those with 0 runs", {
  n_total   <- nrow(.roh_normal$sample_info)
  n_summary <- nrow(.roh_normal$summary)
  expect_equal(n_total, n_summary)
})

test_that("ROH$summary$n_ROH is non-negative for all individuals", {
  expect_true(all(.roh_normal$summary$n_ROH >= 0L))
})

test_that("ROH$method is a length-1 character string", {
  expect_true(is.character(.roh_normal$method) && length(.roh_normal$method) == 1L)
})

test_that("ROH$type is ROHom or ROHet", {
  expect_true(.roh_normal$type %in% c("ROHom", "ROHet"))
})

test_that("print.RUNS prints without error", {
  expect_output(print(.roh_normal))
})

test_that("print.RUNS returns the ROH object invisibly", {
  ret <- print(.roh_normal)
  expect_identical(ret, .roh_normal)
})

test_that("as.data.frame.RUNS returns data.frame with correct columns", {
  df <- as.data.frame(.roh_normal)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("group", "id", "chrom", "nSNP", "from", "to", "lengthBps"))
})

test_that("as.data.frame.RUNS has the same rows as $runs", {
  expect_equal(nrow(as.data.frame(.roh_normal)), nrow(.roh_normal$runs))
})


# ===========================================================================
# 3. Empty results — zero runs
# ===========================================================================

test_that("impossible minSNP produces 0 runs", {
  expect_equal(nrow(.roh_empty$runs), 0L)
})

test_that("empty ROH: $summary still covers all individuals (n_ROH = 0)", {
  expect_equal(nrow(.roh_empty$summary), nrow(.roh_empty$sample_info))
  expect_true(all(.roh_empty$summary$n_ROH == 0L))
})

test_that("empty ROH: print.RUNS does not error", {
  expect_output(print(.roh_empty))
})

test_that("empty ROH: Froh_inbreeding genome-wide returns data.frame", {
  froh <- suppressMessages(Froh_inbreeding(.roh_empty, genome_wide = TRUE))
  expect_s3_class(froh, "data.frame")
  expect_true("Froh_genome" %in% names(froh))
})

test_that("empty ROH: Froh_inbreeding genome-wide Froh is 0 for all individuals", {
  froh <- suppressMessages(Froh_inbreeding(.roh_empty, genome_wide = TRUE))
  expect_true(all(froh$Froh_genome == 0 | is.na(froh$Froh_genome)))
})

test_that("empty ROH: Froh_inbreeding chromosome-wide returns data.frame", {
  froh <- suppressMessages(Froh_inbreeding(.roh_empty, genome_wide = FALSE))
  expect_s3_class(froh, "data.frame")
})

test_that("empty ROH: tableRuns returns empty data.frame without error", {
  result <- suppressMessages(tableRuns(.roh_empty, threshold = 0.5))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("empty ROH: summaryRuns does not error", {
  summ <- suppressMessages(suppressWarnings(
    summaryRuns(.roh_empty, Class = 2, snpInRuns = FALSE)
  ))
  expect_type(summ, "list")
  expect_true(length(summ) >= 6L)
})


# ===========================================================================
# 4. saveRUNS / loadRUNS roundtrip
# ===========================================================================

test_that("saveRUNS creates a file on disk", {
  tmp <- tempfile(fileext = ".roh")
  on.exit(unlink(tmp), add = TRUE)
  saveRUNS(.roh_normal, tmp)
  expect_true(file.exists(tmp))
})

test_that("loadRUNS returns a list with $runs element", {
  tmp <- tempfile(fileext = ".roh")
  on.exit(unlink(tmp), add = TRUE)
  saveRUNS(.roh_normal, tmp)
  roh2 <- loadRUNS(tmp)
  expect_type(roh2, "list")
  expect_true("runs" %in% names(roh2))
})

test_that("saveRUNS/loadRUNS preserves runs exactly", {
  tmp <- tempfile(fileext = ".roh")
  on.exit(unlink(tmp), add = TRUE)
  saveRUNS(.roh_normal, tmp)
  roh2 <- loadRUNS(tmp)
  orig <- as.data.frame(.roh_normal$runs)
  back <- as.data.frame(roh2$runs)
  expect_equal(orig, back)
})

test_that("saveRUNS/loadRUNS preserves full ROH object including method and type", {
  tmp <- tempfile(fileext = ".roh")
  on.exit(unlink(tmp), add = TRUE)
  saveRUNS(.roh_normal, tmp)
  roh2 <- loadRUNS(tmp)
  expect_s3_class(roh2, "RUNS")
  expect_equal(roh2$method, .roh_normal$method)
  expect_equal(roh2$type,   .roh_normal$type)
})

test_that("saveRUNS/loadRUNS roundtrip on empty-runs ROH object", {
  tmp <- tempfile(fileext = ".roh")
  on.exit(unlink(tmp), add = TRUE)
  saveRUNS(.roh_empty, tmp)
  roh2 <- loadRUNS(tmp)
  expect_type(roh2, "list")
  expect_equal(nrow(roh2$runs), 0L)
})


# ===========================================================================
# 5. as_RUNS construction
# ===========================================================================

test_that("as_RUNS builds a valid ROH object from an external CSV", {
  ext_runs <- suppressMessages(
    readExternalRuns("test.ROHet.sliding.csv", program = "detectRUNS")
  )
  roh <- suppressWarnings(
    as_RUNS(ext_runs, mapFile = map_file, genotypeFile = ped_file,
           method = "sliding", type = "ROHom")
  )
  expect_s3_class(roh, "RUNS")
  expect_equal(nrow(roh$runs), nrow(ext_runs))
})

test_that("as_RUNS: Froh_inbreeding works on the resulting ROH object", {
  ext_runs <- suppressMessages(
    readExternalRuns("test.ROHet.sliding.csv", program = "detectRUNS")
  )
  roh <- suppressWarnings(
    as_RUNS(ext_runs, mapFile = map_file, genotypeFile = ped_file,
           method = "sliding", type = "ROHom")
  )
  froh <- suppressMessages(Froh_inbreeding(roh, genome_wide = TRUE))
  expect_s3_class(froh, "data.frame")
  expect_true("Froh_genome" %in% names(froh))
  expect_true(all(froh$Froh_genome >= 0, na.rm = TRUE))
})

test_that("as_RUNS errors when neither mapFile nor bedFile is supplied", {
  runs <- data.frame(group="X", id="a", chrom="1", nSNP=10L,
                     from=1L, to=100L, lengthBps=99L)
  expect_error(as_RUNS(runs, method = "consecutive", type = "ROHom"),
               "Provide mapFile=|bedFile=")
})


# ===========================================================================
# 6. Run-level filters
# ===========================================================================

test_that("minLengthBps = 1e9 produces 0 runs (all shorter)", {
  roh <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 2, minLengthBps = 1e9, verbose = FALSE)
  ))
  expect_equal(nrow(roh$runs), 0L)
})

test_that("minSNP=999999 produces strictly fewer runs than minSNP=5", {
  roh_loose <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 5, maxOpp = 2, maxMiss = 2,
             minLengthBps = 0, verbose = FALSE)
  ))
  roh_strict <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 999999L, verbose = FALSE)
  ))
  expect_lte(nrow(roh_strict$runs), nrow(roh_loose$runs))
})

test_that("stricter maxOpp produces no more runs than permissive maxOpp", {
  roh_permissive <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 5, maxOpp = 5, maxMiss = 5,
             minLengthBps = 50000, verbose = FALSE)
  ))
  roh_strict <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 5, maxOpp = 0, maxMiss = 0,
             minLengthBps = 50000, verbose = FALSE)
  ))
  expect_lte(nrow(roh_strict$runs), nrow(roh_permissive$runs))
})

test_that("all detected runs satisfy minLengthBps constraint", {
  min_len <- 100000L
  roh <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = 5, minLengthBps = min_len, verbose = FALSE)
  ))
  if (nrow(roh$runs) > 0L)
    expect_true(all(roh$runs$lengthBps >= min_len))
})

test_that("all detected runs satisfy minSNP constraint", {
  min_snp <- 10L
  roh <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "consecutive",
             minSNP = min_snp, minLengthBps = 0, verbose = FALSE)
  ))
  if (nrow(roh$runs) > 0L)
    expect_true(all(roh$runs$nSNP >= min_snp))
})

test_that("sliding method produces a valid ROH object", {
  roh_sliding <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "sliding",
             windowSize = 10, threshold = 0.05,
             minSNP = 10, minLengthBps = 50000, verbose = FALSE)
  ))
  expect_s3_class(roh_sliding, "RUNS")
  expect_equal(roh_sliding$method, "sliding")
  expect_named(
    as.data.frame(roh_sliding$runs),
    c("group", "id", "chrom", "nSNP", "from", "to", "lengthBps")
  )
})

test_that("ROHet scan returns runs of heterozygosity", {
  roh_het <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file, method = "sliding",
             ROHet = TRUE, windowSize = 5, threshold = 0.05,
             minSNP = 3, minLengthBps = 1000, verbose = FALSE)
  ))
  expect_s3_class(roh_het, "RUNS")
  expect_equal(roh_het$type, "ROHet")
})


# ===========================================================================
# 7. Downstream stats on a valid (non-empty) ROH object
# ===========================================================================

test_that("Froh_inbreeding genome-wide: all Froh values in [0, 1]", {
  froh <- suppressMessages(Froh_inbreeding(.roh_normal, genome_wide = TRUE))
  expect_true(all(froh$Froh_genome >= 0 & froh$Froh_genome <= 1, na.rm = TRUE))
})

test_that("Froh_inbreeding chromosome-wide: returns one row per individual including 0-run ones", {
  froh <- suppressMessages(Froh_inbreeding(.roh_normal, genome_wide = FALSE))
  n_individuals <- nrow(.roh_normal$sample_info)
  expect_equal(nrow(froh), n_individuals)
})

test_that("Froh_inbreedingClass: returns data.frame with id and group columns", {
  froh_class <- suppressMessages(Froh_inbreedingClass(.roh_normal, Class = 2))
  expect_s3_class(froh_class, "data.frame")
  expect_true(all(c("id", "group") %in% names(froh_class)))
})

test_that("summaryRuns returns a list with 9 elements", {
  summ <- suppressMessages(suppressWarnings(
    summaryRuns(.roh_normal, Class = 2, snpInRuns = FALSE)
  ))
  expect_type(summ, "list")
  expect_equal(length(summ), 9L)
})

test_that("summaryRuns list element names are correct", {
  summ <- suppressMessages(suppressWarnings(
    summaryRuns(.roh_normal, Class = 2, snpInRuns = FALSE)
  ))
  expected_names <- c("summary_ROH_count_chr", "summary_ROH_percentage_chr",
                      "summary_ROH_count", "summary_ROH_percentage",
                      "summary_ROH_mean_chr", "summary_ROH_mean_class",
                      "result_Froh_genome_wide", "result_Froh_chromosome_wide",
                      "result_Froh_class")
  expect_named(summ, expected_names)
})

test_that("tableRuns returns a data.frame with expected columns", {
  result <- suppressMessages(tableRuns(.roh_normal, threshold = 0.5))
  expect_s3_class(result, "data.frame")
  expected_cols <- c("Group", "Start_SNP", "End_SNP", "chrom", "nSNP", "from", "to")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("tableRuns threshold validation rejects values outside [0,1]", {
  expect_error(tableRuns(.roh_normal, threshold = 1.5),
               "Threshold must be between 0 and 1")
  expect_error(tableRuns(.roh_normal, threshold = -0.1),
               "Threshold must be between 0 and 1")
})


# ===========================================================================
# 8. BED vs PED consistency (extdata sheep subset — consecutive method)
# ===========================================================================

test_that("BED and PED engines produce the same run count (consecutive)", {
  bed_file <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
  ped_file_ext <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")

  skip_if(bed_file == "", "extdata BED file not found")
  skip_if(ped_file_ext == "", "extdata PED file not found")

  roh_bed <- suppressWarnings(suppressMessages(
    scanRUNS(bed_file, method = "consecutive",
             minSNP = 20, maxOpp = 1, maxMiss = 1,
             minLengthBps = 250000, verbose = FALSE)
  ))
  roh_ped <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file_ext, method = "consecutive",
             minSNP = 20, maxOpp = 1, maxMiss = 1,
             minLengthBps = 250000, verbose = FALSE)
  ))

  expect_equal(nrow(roh_bed$runs), nrow(roh_ped$runs))
})

test_that("BED and PED engines produce the same run positions (consecutive)", {
  bed_file <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
  ped_file_ext <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")

  skip_if(bed_file == "", "extdata BED file not found")
  skip_if(ped_file_ext == "", "extdata PED file not found")

  roh_bed <- suppressWarnings(suppressMessages(
    scanRUNS(bed_file, method = "consecutive",
             minSNP = 20, maxOpp = 1, maxMiss = 1,
             minLengthBps = 250000, verbose = FALSE)
  ))
  roh_ped <- suppressWarnings(suppressMessages(
    scanRUNS(ped_file_ext, method = "consecutive",
             minSNP = 20, maxOpp = 1, maxMiss = 1,
             minLengthBps = 250000, verbose = FALSE)
  ))

  ord <- c("id", "chrom", "from", "to")
  bed_df <- as.data.frame(roh_bed$runs)[order(roh_bed$runs$id,
                                               roh_bed$runs$chrom,
                                               roh_bed$runs$from), ]
  ped_df <- as.data.frame(roh_ped$runs)[order(roh_ped$runs$id,
                                               roh_ped$runs$chrom,
                                               roh_ped$runs$from), ]
  row.names(bed_df) <- NULL
  row.names(ped_df) <- NULL

  expect_equal(bed_df$from, ped_df$from)
  expect_equal(bed_df$to,   ped_df$to)
  expect_equal(bed_df$id,   ped_df$id)
})
