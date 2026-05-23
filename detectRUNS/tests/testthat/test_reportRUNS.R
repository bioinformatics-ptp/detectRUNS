###############################################################################
## Tests for reportRUNS() and scan_params / meta in ROH objects
###############################################################################

library(testthat)
library(detectRUNS)

bed_file <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
ped_file <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
map_file <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")

skip_if_no_data <- function() {
  if (!nzchar(bed_file)) skip("extdata not available")
}

# ---------------------------------------------------------------------------
# scan_params and meta stored in ROH object (BED)
# ---------------------------------------------------------------------------
test_that("scanRUNS BED stores scan_params", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  expect_false(is.null(roh$scan_params))
  sp <- roh$scan_params
  expect_equal(sp$input_format, "bed")
  expect_equal(sp$minSNP, 10L)
  expect_equal(sp$minLengthBps, 50000L)
  expect_false(is.na(sp$genoFile))
  expect_true(is.na(sp$mapFile))
})

test_that("scanRUNS BED stores meta", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  expect_false(is.null(roh$meta))
  m <- roh$meta
  expect_true(nzchar(m$pkg_version))
  expect_true(nzchar(m$r_version))
  expect_true(nzchar(m$timestamp))
  expect_true(nzchar(m$platform))
})

# ---------------------------------------------------------------------------
# scan_params and meta stored in ROH object (PED)
# ---------------------------------------------------------------------------
test_that("scanRUNS PED stores scan_params", {
  skip_if_no_data()
  roh <- scanRUNS(ped_file, mapFile = map_file, method = "sliding",
                  minSNP = 10, minLengthBps = 50000)
  expect_false(is.null(roh$scan_params))
  sp <- roh$scan_params
  expect_equal(sp$input_format, "ped")
  expect_false(is.na(sp$mapFile))
  expect_false(is.na(sp$genoFile))   # ped path stored in genoFile
  expect_true(is.na(sp$bimFile))
  expect_true(is.na(sp$famFile))
})

test_that("scanRUNS PED stores meta", {
  skip_if_no_data()
  roh <- scanRUNS(ped_file, mapFile = map_file, method = "sliding",
                  minSNP = 10, minLengthBps = 50000)
  expect_false(is.null(roh$meta))
  expect_true(nzchar(roh$meta$timestamp))
})

# ---------------------------------------------------------------------------
# print.ROH shows Params and Scanned lines
# ---------------------------------------------------------------------------
test_that("print.ROH shows Params and Scanned lines", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  out <- capture.output(print(roh))
  expect_true(any(grepl("Params", out)))
  expect_true(any(grepl("Scanned", out)))
  expect_true(any(grepl("minSNP", out)))
})

# ---------------------------------------------------------------------------
# reportRUNS basic smoke test
# ---------------------------------------------------------------------------
test_that("reportRUNS produces markdown file and returns list", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  tmp <- tempdir()
  out <- reportRUNS(roh, output_dir = tmp, prefix = "test_md",
                    format = "markdown", include_plots = FALSE,
                    snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
  expect_true(file.exists(out$report_file))
  expect_true(grepl("\\.md$", out$report_file))
  expect_true(is.list(out))
  expect_true(!is.null(out$summary))
  file.remove(out$report_file)
})

test_that("reportRUNS produces html file", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  tmp <- tempdir()
  out <- reportRUNS(roh, output_dir = tmp, prefix = "test_html",
                    format = "html", include_plots = FALSE,
                    snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
  expect_true(file.exists(out$report_file))
  expect_true(grepl("\\.html$", out$report_file))
  html <- readLines(out$report_file)
  expect_true(any(grepl("<html", html, fixed = TRUE)))
  file.remove(out$report_file)
})

test_that("reportRUNS errors when output_dir does not exist", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  expect_error(
    reportRUNS(roh, output_dir = "/nonexistent_xyz_dir", verbose = FALSE),
    "does not exist"
  )
})

test_that("reportRUNS errors when file exists and overwrite=FALSE", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  tmp <- tempdir()
  pfx <- paste0("test_overwrite_", as.integer(Sys.time()))
  out <- reportRUNS(roh, output_dir = tmp, prefix = pfx,
                    format = "markdown", include_plots = FALSE,
                    snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
  expect_error(
    reportRUNS(roh, output_dir = tmp, prefix = pfx,
               format = "markdown", include_plots = FALSE,
               snp_table = FALSE, froh_class = FALSE,
               overwrite = FALSE, verbose = FALSE),
    "already exists"
  )
  file.remove(out$report_file)
})

test_that("reportRUNS overwrite=TRUE replaces existing file", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  tmp <- tempdir()
  pfx <- paste0("test_overwrite2_", as.integer(Sys.time()))
  out1 <- reportRUNS(roh, output_dir = tmp, prefix = pfx,
                     format = "markdown", include_plots = FALSE,
                     snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
  expect_no_error(
    reportRUNS(roh, output_dir = tmp, prefix = pfx,
               format = "markdown", include_plots = FALSE,
               snp_table = FALSE, froh_class = FALSE,
               overwrite = TRUE, verbose = FALSE)
  )
  file.remove(out1$report_file)
})

test_that("reportRUNS errors when 'runs' is not an ROH object", {
  expect_error(reportRUNS(data.frame()), "ROH object")
})

# ---------------------------------------------------------------------------
# reportRUNS with plots
# ---------------------------------------------------------------------------
test_that("reportRUNS generates plot files when include_plots=TRUE", {
  skip_if_no_data()
  roh <- scanRUNS(bed_file, method = "sliding", minSNP = 10, minLengthBps = 50000)
  tmp <- tempdir()
  out <- reportRUNS(roh, output_dir = tmp, prefix = "test_plots",
                    format = "markdown", include_plots = TRUE,
                    snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
  expect_true(length(out$plots) > 0L)
  for (p in out$plots) expect_true(file.exists(p))
  file.remove(out$report_file)
  unlink(out$plot_dir, recursive = TRUE)
})
