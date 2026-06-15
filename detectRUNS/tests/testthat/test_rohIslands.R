## Tests for runsIslands() — permutation-based ROH island detection
##
## Primary fixtures use test.bed (20 ind, 563 SNPs) — fast, CRAN-safe.
## n_perm = 3 throughout except for one skip_on_cran() integration test.

library(testthat)
library(detectRUNS)

context("runsIslands")

# ===========================================================================
# Shared fixtures — test.bed only (fast, local, CRAN-safe)
# ===========================================================================

# BED scan (consecutive) — provides snp_freq required by runsIslands
.roh_bed <- local({
    suppressWarnings(suppressMessages(
        scanRUNS("test.bed", method = "consecutive",
                 minSNP = 5, maxOpp = 1, maxMiss = 1,
                 minLengthBps = 10000, verbose = FALSE)
    ))
})

# BED scan (sliding) — for method-coverage tests
.roh_sliding <- local({
    suppressWarnings(suppressMessages(
        scanRUNS("test.bed", method = "sliding",
                 windowSize = 5, threshold = 0.05,
                 minSNP = 5, maxOpp = 1, maxMiss = 1,
                 minLengthBps = 10000, verbose = FALSE)
    ))
})

# PED scan — no snp_freq (BED engine needed for that)
.roh_ped <- local({
    suppressWarnings(suppressMessages(
        scanRUNS("test.ped", method = "consecutive",
                 minSNP = 5, maxOpp = 1, maxMiss = 1,
                 minLengthBps = 10000, verbose = FALSE)
    ))
})

# Pre-computed reference — n_perm=3, seed=1
.islands_ref <- local({
    suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 1L, verbose = FALSE)
    )
})


# ===========================================================================
# 1. Input validation
# ===========================================================================

test_that("runsIslands rejects non-ROH input", {
    expect_error(runsIslands(list(runs = data.frame())), "must be a RUNS object")
    expect_error(runsIslands(data.frame()),              "must be a RUNS object")
    expect_error(runsIslands("not an ROH"),              "must be a RUNS object")
})

test_that("runsIslands rejects PED-sourced ROH (no snp_freq)", {
    expect_error(runsIslands(.roh_ped), "snp_freq.*NULL|BED-format")
})

test_that("runsIslands rejects ROH with missing scan_params", {
    roh_no_sp <- .roh_bed
    roh_no_sp$scan_params <- NULL
    expect_error(runsIslands(roh_no_sp), "scan_params.*NULL|re-run scanRUNS")
})

test_that("runsIslands rejects non-existent bed_path", {
    expect_error(runsIslands(.roh_bed, bed_path = "no_such_file.bed"),
                 "not found|does not exist")
})

test_that("runsIslands rejects invalid n_perm", {
    expect_error(runsIslands(.roh_bed, n_perm = 0,   verbose = FALSE), "n_perm")
    expect_error(runsIslands(.roh_bed, n_perm = -1,  verbose = FALSE), "n_perm")
    expect_error(runsIslands(.roh_bed, n_perm = "a", verbose = FALSE), "n_perm")
})

test_that("runsIslands rejects percentile outside (0, 1)", {
    expect_error(runsIslands(.roh_bed, percentile = 0,    verbose = FALSE), "percentile")
    expect_error(runsIslands(.roh_bed, percentile = 1,    verbose = FALSE), "percentile")
    expect_error(runsIslands(.roh_bed, percentile = -0.1, verbose = FALSE), "percentile")
    expect_error(runsIslands(.roh_bed, percentile = 1.5,  verbose = FALSE), "percentile")
})

test_that("runsIslands rejects non-scalar seed", {
    expect_error(runsIslands(.roh_bed, seed = c(1, 2), verbose = FALSE), "seed")
})


# ===========================================================================
# 2. Return structure
# ===========================================================================

test_that("runsIslands returns RunsIslands S3 object", {
    expect_s3_class(.islands_ref, "RunsIslands")
})

test_that("RunsIslands has all required fields", {
    expect_true(all(c("islands", "snp_table", "thresholds",
                      "n_samples", "n_perm", "percentile") %in% names(.islands_ref)))
})

test_that("$islands is a data.table", {
    expect_s3_class(.islands_ref$islands, "data.table")
})

test_that("$snp_table is a data.table", {
    expect_s3_class(.islands_ref$snp_table, "data.table")
})

test_that("$islands has required columns", {
    expected <- c("SNP_NAME", "CHR", "POSITION",
                  "snp_freq", "threshold", "is_island", "pct_animals")
    expect_true(all(expected %in% names(.islands_ref$islands)))
})

test_that("$snp_table has required columns", {
    expected <- c("SNP_NAME", "CHR", "POSITION",
                  "snp_freq", "threshold", "is_island", "pct_animals")
    expect_true(all(expected %in% names(.islands_ref$snp_table)))
})

test_that("$thresholds is a named numeric vector", {
    expect_type(.islands_ref$thresholds, "double")
    expect_false(is.null(names(.islands_ref$thresholds)))
    expect_true(length(.islands_ref$thresholds) > 0L)
})

test_that("$n_perm and $percentile reflect inputs", {
    expect_equal(.islands_ref$n_perm, 3L)
    expect_equal(.islands_ref$percentile, 0.99)
})

test_that("$n_samples matches ROH sample count", {
    expect_equal(.islands_ref$n_samples, nrow(.roh_bed$sample_info))
})

test_that("snp_table row count matches SNPs in BIM", {
    expect_equal(nrow(.islands_ref$snp_table), nrow(.roh_bed$snp_map))
})

test_that("islands is exactly the is_island == TRUE rows of snp_table", {
    expect_equal(nrow(.islands_ref$islands),
                 sum(.islands_ref$snp_table$is_island))
})


# ===========================================================================
# 3. Content correctness
# ===========================================================================

test_that("all island SNPs have snp_freq > threshold", {
    if (nrow(.islands_ref$islands) == 0L) skip("no islands with n_perm=3")
    isl <- .islands_ref$islands
    expect_true(all(isl$snp_freq > isl$threshold))
})

test_that("non-island SNPs have snp_freq <= threshold", {
    non_isl <- .islands_ref$snp_table[.islands_ref$snp_table$is_island == FALSE, ]
    if (nrow(non_isl) == 0L) skip("all SNPs are islands")
    expect_true(all(non_isl$snp_freq <= non_isl$threshold))
})

test_that("thresholds are non-negative", {
    expect_true(all(.islands_ref$thresholds >= 0))
})

test_that("pct_animals is in [0, 100]", {
    pct <- .islands_ref$snp_table$pct_animals
    expect_true(all(pct >= 0 & pct <= 100))
})

test_that("snp_freq values in snp_table match roh$snp_freq", {
    expect_equal(.islands_ref$snp_table$snp_freq,
                 as.integer(.roh_bed$snp_freq))
})

test_that("chromosome names in thresholds match CHR column of snp_table", {
    chrs_table  <- sort(unique(as.character(.islands_ref$snp_table$CHR)))
    chrs_thresh <- sort(names(.islands_ref$thresholds))
    expect_equal(chrs_table, chrs_thresh)
})

test_that("per-SNP threshold matches its chromosome entry in $thresholds", {
    df <- as.data.frame(.islands_ref$snp_table)
    for (chr in names(.islands_ref$thresholds)) {
        rows <- df[as.character(df$CHR) == chr, ]
        if (nrow(rows) == 0L) next
        expect_true(all(rows$threshold == .islands_ref$thresholds[[chr]]))
    }
})


# ===========================================================================
# 4. Determinism and seed
# ===========================================================================

test_that("same seed produces identical results", {
    isl2 <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 1L, verbose = FALSE)
    )
    expect_equal(.islands_ref$thresholds, isl2$thresholds)
    expect_equal(.islands_ref$islands,    isl2$islands)
})

test_that("seed = 0 runs without error (random-seed path)", {
    expect_no_error(suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 0L, verbose = FALSE)
    ))
})

test_that("n_perm = 1 completes without error", {
    isl1 <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 1L, seed = 42L, verbose = FALSE)
    )
    expect_s3_class(isl1, "RunsIslands")
    expect_true(all(isl1$thresholds >= 0))
})


# ===========================================================================
# 5. Monotonicity
# ===========================================================================

test_that("lower percentile gives >= island count than higher percentile", {
    isl_hi <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 7L, percentile = 0.99,
                    verbose = FALSE)
    )
    isl_lo <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 7L, percentile = 0.50,
                    verbose = FALSE)
    )
    expect_gte(nrow(isl_lo$islands), nrow(isl_hi$islands))
})

test_that("lower percentile gives <= thresholds than higher percentile", {
    isl_hi <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 7L, percentile = 0.99,
                    verbose = FALSE)
    )
    isl_lo <- suppressWarnings(
        runsIslands(.roh_bed, n_perm = 3L, seed = 7L, percentile = 0.10,
                    verbose = FALSE)
    )
    for (nm in names(isl_hi$thresholds)) {
        expect_lte(isl_lo$thresholds[[nm]], isl_hi$thresholds[[nm]])
    }
})


# ===========================================================================
# 6. Method coverage (sliding window)
# ===========================================================================

test_that("runsIslands works with sliding window method", {
    isl_sw <- suppressWarnings(
        runsIslands(.roh_sliding, n_perm = 3L, seed = 3L, verbose = FALSE)
    )
    expect_s3_class(isl_sw, "RunsIslands")
    expect_true(all(c("islands", "snp_table", "thresholds") %in% names(isl_sw)))
    expect_true(all(isl_sw$thresholds >= 0))
})

test_that("sliding and consecutive give same SNP count in snp_table", {
    isl_sw <- suppressWarnings(
        runsIslands(.roh_sliding, n_perm = 3L, seed = 3L, verbose = FALSE)
    )
    expect_equal(nrow(isl_sw$snp_table), nrow(.islands_ref$snp_table))
})


# ===========================================================================
# 7. Explicit bed_path override
# ===========================================================================

test_that("runsIslands works when bed_path passed explicitly", {
    roh_no_bp <- .roh_bed
    roh_no_bp$bed_path <- NULL
    isl <- suppressWarnings(
        runsIslands(roh_no_bp, bed_path = "test.bed",
                    n_perm = 3L, seed = 1L, verbose = FALSE)
    )
    expect_equal(isl$thresholds, .islands_ref$thresholds)
})


# ===========================================================================
# 8. print.RunsIslands
# ===========================================================================

test_that("print.RunsIslands produces output without error", {
    expect_output(print(.islands_ref), "RunsIslands")
    expect_output(print(.islands_ref), "n_perm")
    expect_output(print(.islands_ref), "percentile")
})

test_that("print.RunsIslands invisibly returns x", {
    out <- withVisible(print(.islands_ref))
    expect_false(out$visible)
    expect_identical(out$value, .islands_ref)
})


# ===========================================================================
# 9. summary.RunsIslands
# ===========================================================================

test_that("summary.RunsIslands returns a data.table", {
    expect_s3_class(summary(.islands_ref), "data.table")
})

test_that("summary.RunsIslands has required columns", {
    out <- summary(.islands_ref)
    expect_true(all(c("CHR", "start_bp", "end_bp", "n_snps",
                      "peak_pct", "width_mb") %in% names(out)))
})

test_that("summary.RunsIslands start_bp <= end_bp", {
    out <- summary(.islands_ref)
    if (nrow(out) == 0L) skip("no islands")
    expect_true(all(out$start_bp <= out$end_bp))
})

test_that("summary.RunsIslands n_snps >= 1", {
    out <- summary(.islands_ref)
    if (nrow(out) == 0L) skip("no islands")
    expect_true(all(out$n_snps >= 1L))
})

test_that("summary.RunsIslands total SNP count matches $islands row count", {
    out <- summary(.islands_ref)
    if (nrow(out) == 0L) skip("no islands")
    expect_equal(sum(out$n_snps), nrow(.islands_ref$islands))
})

test_that("summary.RunsIslands returns empty data.table when no islands", {
    roh_no_isl <- .islands_ref
    snp_copy <- as.data.frame(.islands_ref$snp_table)
    snp_copy$is_island <- FALSE
    roh_no_isl$snp_table <- data.table::as.data.table(snp_copy)
    out <- summary(roh_no_isl)
    expect_s3_class(out, "data.table")
    expect_equal(nrow(out), 0L)
})


# ===========================================================================
# 10. plot.RunsIslands
# ===========================================================================

.with_null_dev <- function(expr) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
    force(expr)
}

test_that("plot.RunsIslands returns a ggplot object", {
    .with_null_dev({
        p <- plot(.islands_ref)
        expect_s3_class(p, "gg")
    })
})

test_that("plot.RunsIslands returns invisibly", {
    .with_null_dev({
        out <- withVisible(plot(.islands_ref))
        expect_false(out$visible)
    })
})

test_that("plot.RunsIslands accepts custom colours without error", {
    .with_null_dev({
        expect_no_error(plot(.islands_ref,
                             col_island    = "darkgreen",
                             col_snp       = c("black", "grey40"),
                             col_threshold = "orange"))
    })
})

test_that("plot.RunsIslands works when no islands are present", {
    roh_no_isl <- .islands_ref
    snp_copy <- as.data.frame(.islands_ref$snp_table)
    snp_copy$is_island <- FALSE
    roh_no_isl$snp_table <- data.table::as.data.table(snp_copy)
    roh_no_isl$islands   <- roh_no_isl$snp_table[integer(0), ]
    .with_null_dev(expect_no_error(plot(roh_no_isl)))
})


# ===========================================================================
# 11. Integration test — larger dataset (skip on CRAN)
# ===========================================================================

test_that("runsIslands with extdata Kijas and n_perm=50 produces stable islands", {
    skip_on_cran()
    bed <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
    skip_if(bed == "", "extdata not found")
    roh <- suppressWarnings(suppressMessages(
        scanRUNS(bed, method = "consecutive",
                 minSNP = 15, maxOpp = 1, maxMiss = 1,
                 minLengthBps = 100000, verbose = FALSE)
    ))
    isl <- suppressWarnings(
        runsIslands(roh, n_perm = 50L, seed = 42L, verbose = FALSE)
    )
    expect_s3_class(isl, "RunsIslands")
    expect_true(all(isl$thresholds >= 0))
    if (nrow(isl$islands) > 0L) {
        expect_true(all(isl$islands$snp_freq > isl$islands$threshold))
    }
})
