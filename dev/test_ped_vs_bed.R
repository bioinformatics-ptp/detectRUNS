###############################################################################
## PED (R algorithm) vs BED (C++ engine) — exhaustive parameter comparison
##
## Dataset: package test data (20 individuals, 563 SNPs, chr 24)
##   detectRUNS/tests/testthat/test.ped / test.map / test.bed
##
## For each scenario, scanRUNS() is run twice:
##   - PED path  → R-based algorithm (funktionen.R)
##   - BED path  → C++ engine (scan_roh.cpp)
##
## A scenario PASSES when both paths return identical runs
## (same rows after sorting on id + chrom + from + to).
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L) {
    project_root <- dirname(dirname(normalizePath(script_path)))
    setwd(project_root)
}

PED <- "detectRUNS/tests/testthat/test.ped"
MAP <- "detectRUNS/tests/testthat/test.map"
BED <- "detectRUNS/tests/testthat/test.bed"

stopifnot(file.exists(PED), file.exists(MAP), file.exists(BED))
cat("Test data: 20 individuals | 563 SNPs | chr 24\n\n")

# ---------------------------------------------------------------------------
# Parameter scenarios
# Covers: ROHom/ROHet, sliding/consecutive,
#         loose → strict on every axis, maxOppRun/maxMissRun, minDensity
# ---------------------------------------------------------------------------
SCENARIOS <- list(

    # ---- sliding ROHom ----
    list(label="sliding_ROHom_default",
         method="sliding", ROHet=FALSE,
         minSNP=3, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=1000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_strict",
         method="sliding", ROHet=FALSE,
         minSNP=20, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=250000, maxGap=1e6,
         minDensity=1/1000, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_veryStrict",
         method="sliding", ROHet=FALSE,
         minSNP=30, maxOpp=0, maxMiss=0, windowSize=20,
         threshold=0.05, minLengthBps=500000, maxGap=5e5,
         minDensity=1/1000, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_lenient",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.10, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_highThreshold",
         method="sliding", ROHet=FALSE,
         minSNP=10, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.50, minLengthBps=100000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_smallWindow",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=1, maxMiss=1, windowSize=5,
         threshold=0.05, minLengthBps=50000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_largeWindow",
         method="sliding", ROHet=FALSE,
         minSNP=10, maxOpp=2, maxMiss=2, windowSize=30,
         threshold=0.05, minLengthBps=100000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_maxOppRun1",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=1, maxMissRun=NULL),

    list(label="sliding_ROHom_maxMissRun2",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=2),

    list(label="sliding_ROHom_maxOppRun0_maxMissRun0",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=0, maxMissRun=0),

    list(label="sliding_ROHom_density_1per50",
         method="sliding", ROHet=FALSE,
         minSNP=5, maxOpp=2, maxMiss=2, windowSize=15,
         threshold=0.05, minLengthBps=100000, maxGap=1.5e6,
         minDensity=1/50, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHom_tinyLength",
         method="sliding", ROHet=FALSE,
         minSNP=2, maxOpp=1, maxMiss=1, windowSize=5,
         threshold=0.05, minLengthBps=100, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    # maxGap=0 is rejected by the PED path ("maxGap must be >= 1"); skip.

    # ---- sliding ROHet ----
    list(label="sliding_ROHet_default",
         method="sliding", ROHet=TRUE,
         minSNP=3, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=1000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHet_strict",
         method="sliding", ROHet=TRUE,
         minSNP=15, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=100000, maxGap=1e6,
         minDensity=1/1000, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHet_lenient",
         method="sliding", ROHet=TRUE,
         minSNP=3, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.10, minLengthBps=10000, maxGap=2e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="sliding_ROHet_maxOppRun1",
         method="sliding", ROHet=TRUE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=1, maxMissRun=NULL),

    # ---- consecutive ROHom ----
    list(label="consecutive_ROHom_default",
         method="consecutive", ROHet=FALSE,
         minSNP=3, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=1000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHom_strict",
         method="consecutive", ROHet=FALSE,
         minSNP=20, maxOpp=0, maxMiss=0, windowSize=15,
         threshold=0.05, minLengthBps=250000, maxGap=1e6,
         minDensity=1/1000, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHom_lenient",
         method="consecutive", ROHet=FALSE,
         minSNP=5, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=50000, maxGap=2e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHom_zeroOppMiss",
         method="consecutive", ROHet=FALSE,
         minSNP=5, maxOpp=0, maxMiss=0, windowSize=15,
         threshold=0.05, minLengthBps=50000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHom_largeGap",
         method="consecutive", ROHet=FALSE,
         minSNP=5, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=50000, maxGap=5e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    # ---- consecutive ROHet ----
    list(label="consecutive_ROHet_default",
         method="consecutive", ROHet=TRUE,
         minSNP=3, maxOpp=1, maxMiss=1, windowSize=15,
         threshold=0.05, minLengthBps=1000, maxGap=1e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHet_strict",
         method="consecutive", ROHet=TRUE,
         minSNP=15, maxOpp=0, maxMiss=0, windowSize=15,
         threshold=0.05, minLengthBps=100000, maxGap=1e6,
         minDensity=1/1000, maxOppRun=NULL, maxMissRun=NULL),

    list(label="consecutive_ROHet_lenient",
         method="consecutive", ROHet=TRUE,
         minSNP=3, maxOpp=3, maxMiss=3, windowSize=10,
         threshold=0.05, minLengthBps=10000, maxGap=2e6,
         minDensity=0, maxOppRun=NULL, maxMissRun=NULL)
)

cat(sprintf("Running %d scenarios...\n\n", length(SCENARIOS)))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
.sort_runs <- function(df) {
    df <- as.data.frame(df)
    df$chrom <- as.character(df$chrom)
    df[order(df$id, df$chrom, df$from, df$to), ]
}

.run_scenario <- function(s, input) {
    tryCatch(
        scanRUNS(
            genoFile     = input,
            mapFile      = if (grepl("\\.ped$", input)) MAP else NULL,
            method       = s$method,
            ROHet        = s$ROHet,
            minSNP       = s$minSNP,
            maxOpp       = s$maxOpp,
            maxMiss      = s$maxMiss,
            minLengthBps = s$minLengthBps,
            maxGap       = s$maxGap,
            windowSize   = s$windowSize,
            threshold    = s$threshold,
            minDensity   = s$minDensity,
            maxOppRun    = s$maxOppRun,
            maxMissRun   = s$maxMissRun,
            nThreads     = 1L,
            verbose      = FALSE
        ),
        error = function(e) list(error = conditionMessage(e))
    )
}

.compare <- function(ped_runs, bed_runs) {
    p <- .sort_runs(ped_runs)
    b <- .sort_runs(bed_runs)
    row.names(p) <- NULL
    row.names(b) <- NULL

    if (nrow(p) != nrow(b)) {
        return(list(pass = FALSE,
                    reason = sprintf("row count: PED=%d  BED=%d", nrow(p), nrow(b)),
                    only_ped = NULL, only_bed = NULL))
    }

    pk <- with(p, paste(id, chrom, from, to, sep = "\t"))
    bk <- with(b, paste(id, chrom, from, to, sep = "\t"))

    only_ped <- setdiff(pk, bk)
    only_bed <- setdiff(bk, pk)

    if (length(only_ped) == 0 && length(only_bed) == 0) {
        # Value comparison (not identical() — types may differ: numeric vs integer)
        mismatches <- which(as.integer(p$nSNP)      != as.integer(b$nSNP) |
                            as.integer(p$lengthBps) != as.integer(b$lengthBps))
        if (length(mismatches) > 0) {
            return(list(pass = FALSE,
                        reason = sprintf("%d rows differ on nSNP/lengthBps", length(mismatches)),
                        only_ped = NULL, only_bed = NULL))
        }
        return(list(pass = TRUE, reason = "exact match", only_ped = NULL, only_bed = NULL))
    }

    list(pass = FALSE,
         reason = sprintf("only_PED=%d  only_BED=%d", length(only_ped), length(only_bed)),
         only_ped = head(only_ped, 5),
         only_bed = head(only_bed, 5))
}

# ---------------------------------------------------------------------------
# Run all scenarios
# ---------------------------------------------------------------------------
results <- vector("list", length(SCENARIOS))

for (i in seq_along(SCENARIOS)) {
    s <- SCENARIOS[[i]]
    cat(sprintf("[%2d/%d] %-45s ", i, length(SCENARIOS), s$label))

    ped_res <- .run_scenario(s, PED)
    bed_res <- .run_scenario(s, BED)

    if (!is.null(ped_res$error) || !is.null(bed_res$error)) {
        err <- if (!is.null(ped_res$error)) paste("PED:", ped_res$error) else
                                              paste("BED:", bed_res$error)
        cat("ERROR  —", err, "\n")
        results[[i]] <- list(label=s$label, status="ERROR", detail=err,
                             n_ped=NA, n_bed=NA)
        next
    }

    cmp <- .compare(ped_res$runs, bed_res$runs)
    n_ped <- nrow(ped_res$runs)
    n_bed <- nrow(bed_res$runs)

    if (cmp$pass) {
        cat(sprintf("PASS   runs=%d\n", n_ped))
    } else {
        cat(sprintf("FAIL   PED=%d  BED=%d  [%s]\n", n_ped, n_bed, cmp$reason))
        if (!is.null(cmp$only_ped) && length(cmp$only_ped) > 0) {
            cat("         only_PED:", paste(head(cmp$only_ped, 3), collapse=" | "), "\n")
        }
        if (!is.null(cmp$only_bed) && length(cmp$only_bed) > 0) {
            cat("         only_BED:", paste(head(cmp$only_bed, 3), collapse=" | "), "\n")
        }
    }

    results[[i]] <- list(label   = s$label,
                         status  = if (cmp$pass) "PASS" else "FAIL",
                         detail  = cmp$reason,
                         n_ped   = n_ped,
                         n_bed   = n_bed)
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_pass  <- sum(sapply(results, function(r) r$status == "PASS"))
n_fail  <- sum(sapply(results, function(r) r$status == "FAIL"))
n_error <- sum(sapply(results, function(r) r$status == "ERROR"))

cat(sprintf("\n%s\n RESULTS: %d PASS  |  %d FAIL  |  %d ERROR  (of %d scenarios)\n%s\n",
    strrep("=", 65), n_pass, n_fail, n_error, length(SCENARIOS), strrep("=", 65)))

if (n_fail > 0 || n_error > 0) {
    cat("\nFailed / errored scenarios:\n")
    for (r in results) {
        if (r$status != "PASS")
            cat(sprintf("  %-45s  [%s]  %s\n", r$label, r$status, r$detail))
    }
}

# ---------------------------------------------------------------------------
# Write report
# ---------------------------------------------------------------------------
RES_DIR <- "Ext_Data/results/ped_vs_bed"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

rpt_path <- file.path(RES_DIR, "ped_vs_bed_report.md")
rpt <- file(rpt_path, "w")
.w <- function(...) cat(..., "\n", file = rpt, sep = "")

.w("# PED (R algorithm) vs BED (C++ engine) — Comparison Report")
.w()
.w(sprintf("**Date:** %s  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
.w(sprintf("**Package:** detectRUNS %s  ", as.character(utils::packageVersion("detectRUNS"))))
.w(sprintf("**Dataset:** test.ped/bed — 20 individuals | 563 SNPs | chr 24  "))
.w()
.w(sprintf("**PASS:** %d  |  **FAIL:** %d  |  **ERROR:** %d  |  **Total:** %d",
           n_pass, n_fail, n_error, length(SCENARIOS)))
.w()
.w("| # | Scenario | Status | PED runs | BED runs | Detail |")
.w("|---|----------|--------|----------|----------|--------|")
for (i in seq_along(results)) {
    r <- results[[i]]
    .w(sprintf("| %d | %s | **%s** | %s | %s | %s |",
               i, r$label, r$status,
               ifelse(is.na(r$n_ped), "—", r$n_ped),
               ifelse(is.na(r$n_bed), "—", r$n_bed),
               r$detail))
}
close(rpt)

cat(sprintf("\nReport: %s\n", rpt_path))
cat(sprintf("Done at %s\n", format(Sys.time(), "%H:%M:%S")))
