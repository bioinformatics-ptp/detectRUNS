## Tests for reportRUNS() and scan_params / meta stored in ROH objects
##
## All base tests use test.bed / test.ped (local, fast, CRAN-safe).
## Tests that involve full report rendering are also run on CRAN but are
## cheap (test.bed is 20 ind × 563 SNPs; report uses include_plots=FALSE).

library(testthat)
library(detectRUNS)

context("reportRUNS and scan metadata")

# ---------------------------------------------------------------------------
# Shared fixtures — local files only (fast, CRAN-safe)
# ---------------------------------------------------------------------------
.roh_bed <- local({
    suppressWarnings(suppressMessages(
        scanRUNS("test.bed", method = "sliding",
                 windowSize = 5, threshold = 0.05,
                 minSNP = 5, minLengthBps = 10000, verbose = FALSE)
    ))
})

.roh_ped <- local({
    suppressWarnings(suppressMessages(
        scanRUNS("test.ped", mapFile = "test.map",
                 method = "sliding",
                 windowSize = 5, threshold = 0.05,
                 minSNP = 5, minLengthBps = 10000, verbose = FALSE)
    ))
})


# ---------------------------------------------------------------------------
# scan_params stored in ROH object
# ---------------------------------------------------------------------------
test_that("scanRUNS BED stores scan_params", {
    expect_false(is.null(.roh_bed$scan_params))
    sp <- .roh_bed$scan_params
    expect_equal(sp$input_format, "bed")
    expect_equal(sp$minSNP, 5L)
    expect_equal(sp$minLengthBps, 10000L)
    expect_false(is.na(sp$genoFile))
    expect_true(is.na(sp$mapFile))
})

test_that("scanRUNS PED stores scan_params", {
    expect_false(is.null(.roh_ped$scan_params))
    sp <- .roh_ped$scan_params
    expect_equal(sp$input_format, "ped")
    expect_false(is.na(sp$mapFile))
    expect_false(is.na(sp$genoFile))
    expect_true(is.na(sp$bimFile))
    expect_true(is.na(sp$famFile))
})

test_that("scanRUNS BED stores meta with version, timestamp, platform", {
    expect_false(is.null(.roh_bed$meta))
    m <- .roh_bed$meta
    expect_true(nzchar(m$pkg_version))
    expect_true(nzchar(m$r_version))
    expect_true(nzchar(m$timestamp))
    expect_true(nzchar(m$platform))
})

test_that("scanRUNS PED stores meta", {
    expect_false(is.null(.roh_ped$meta))
    expect_true(nzchar(.roh_ped$meta$timestamp))
})

# ---------------------------------------------------------------------------
# print.RUNS shows Params and Scanned lines
# ---------------------------------------------------------------------------
test_that("print.RUNS shows Params and Scanned lines", {
    out <- capture.output(print(.roh_bed))
    expect_true(any(grepl("Params",  out)))
    expect_true(any(grepl("Scanned", out)))
    expect_true(any(grepl("minSNP",  out)))
})

# ---------------------------------------------------------------------------
# reportRUNS — writes to tempdir, CRAN-safe
# ---------------------------------------------------------------------------
test_that("reportRUNS produces a markdown file", {
    out <- reportRUNS(.roh_bed, output_dir = tempdir(), prefix = "cran_test_md",
                      format = "markdown", include_plots = FALSE,
                      snp_table = FALSE, froh_class = FALSE,
                      overwrite = TRUE, verbose = FALSE)
    expect_true(file.exists(out$report_file))
    expect_true(grepl("\\.md$", out$report_file))
    expect_type(out, "list")
    expect_false(is.null(out$summary))
    unlink(out$report_file)
})

test_that("reportRUNS produces an HTML file", {
    out <- reportRUNS(.roh_bed, output_dir = tempdir(), prefix = "cran_test_html",
                      format = "html", include_plots = FALSE,
                      snp_table = FALSE, froh_class = FALSE,
                      overwrite = TRUE, verbose = FALSE)
    expect_true(file.exists(out$report_file))
    expect_true(grepl("\\.html$", out$report_file))
    html <- readLines(out$report_file, warn = FALSE)
    expect_true(any(grepl("<html", html, fixed = TRUE)))
    unlink(out$report_file)
})

test_that("reportRUNS errors when output_dir does not exist", {
    expect_error(
        reportRUNS(.roh_bed, output_dir = "/nonexistent_xyz_dir", verbose = FALSE),
        "does not exist"
    )
})

test_that("reportRUNS errors when file exists and overwrite=FALSE", {
    pfx <- paste0("cran_ow_", as.integer(Sys.time()))
    out <- reportRUNS(.roh_bed, output_dir = tempdir(), prefix = pfx,
                      format = "markdown", include_plots = FALSE,
                      snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
    expect_error(
        reportRUNS(.roh_bed, output_dir = tempdir(), prefix = pfx,
                   format = "markdown", include_plots = FALSE,
                   snp_table = FALSE, froh_class = FALSE,
                   overwrite = FALSE, verbose = FALSE),
        "already exists"
    )
    unlink(out$report_file)
})

test_that("reportRUNS overwrite=TRUE replaces existing file", {
    pfx <- paste0("cran_ow2_", as.integer(Sys.time()))
    out1 <- reportRUNS(.roh_bed, output_dir = tempdir(), prefix = pfx,
                       format = "markdown", include_plots = FALSE,
                       snp_table = FALSE, froh_class = FALSE, verbose = FALSE)
    expect_no_error(
        reportRUNS(.roh_bed, output_dir = tempdir(), prefix = pfx,
                   format = "markdown", include_plots = FALSE,
                   snp_table = FALSE, froh_class = FALSE,
                   overwrite = TRUE, verbose = FALSE)
    )
    unlink(out1$report_file)
})

test_that("reportRUNS errors when 'runs' is not a RUNS object", {
    expect_error(reportRUNS(data.frame()), "RUNS object")
})

# ---------------------------------------------------------------------------
# reportRUNS with plots — more expensive, skip on CRAN
# ---------------------------------------------------------------------------
test_that("reportRUNS generates plot files when include_plots=TRUE", {
    skip_on_cran()
    out <- reportRUNS(.roh_bed, output_dir = tempdir(), prefix = "cran_plots",
                      format = "markdown", include_plots = TRUE,
                      snp_table = FALSE, froh_class = FALSE,
                      overwrite = TRUE, verbose = FALSE)
    expect_true(length(out$plots) > 0L)
    for (p in out$plots) expect_true(file.exists(p))
    unlink(out$report_file)
    unlink(out$plot_dir, recursive = TRUE)
})
