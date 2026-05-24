###############################################################################
## Smoke test — all Ext_Data datasets, 2 parameter sets, BED path
## Goal: no crashes, sane run counts, timing overview
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))

N_CORES <- parallel::detectCores(logical = FALSE)
.ts  <- function() format(Sys.time(), "[%H:%M:%S]")
.sec <- function(t) round(as.numeric(t["elapsed"]), 2)

datasets <- list(
    list(name = "SELMOL",       bed = "Ext_Data/SELMOL_codACGT.bed",                 n_chr = 29L),
    list(name = "ADAPTmap",     bed = "Ext_Data/ADAPTmap_genotypeTOP_20161201.bed",   n_chr = 31L),
    list(name = "suini_12",     bed = "Ext_Data/suini_12_plink.bed",                  n_chr = 18L),
    list(name = "Innovagen_HD", bed = "Ext_Data/Innovagen_HD.bed",                    n_chr = 32L)
)

param_sets <- list(
    list(label = "default",
         windowSize = 15L, minSNP = 15L, maxOpp = 1L, maxMiss = 1L,
         threshold = 0.05, minLengthBps = 500000L, maxGap = 5000000L,
         minDensity = 1/1000),
    list(label = "strict",
         windowSize = 20L, minSNP = 20L, maxOpp = 1L, maxMiss = 1L,
         threshold = 0.05, minLengthBps = 1000000L, maxGap = 5000000L,
         minDensity = 1/50)
)

cat(sprintf("\n%s  SMOKE TEST — all Ext_Data datasets (%d cores)\n", .ts(), N_CORES))
cat(strrep("=", 72), "\n")

results <- list()
n_pass <- 0L; n_fail <- 0L

for (ds in datasets) {
    bim  <- sub("\\.bed$", ".bim", ds$bed)
    snps <- nrow(read.table(bim, header = FALSE))
    inds <- nrow(read.table(sub("\\.bed$", ".fam", ds$bed), header = FALSE))
    cat(sprintf("\n%s  %s  (%d SNPs | %d ind)\n", .ts(), ds$name, snps, inds))

    for (ps in param_sets) {
        status <- "PASS"; err_msg <- ""
        t <- tryCatch({
            system.time(
                res <- scanRUNS(
                    genoFile     = ds$bed,
                    method       = "sliding",
                    windowSize   = ps$windowSize,
                    minSNP       = ps$minSNP,
                    maxOpp       = ps$maxOpp,
                    maxMiss      = ps$maxMiss,
                    threshold    = ps$threshold,
                    minLengthBps = ps$minLengthBps,
                    maxGap       = ps$maxGap,
                    minDensity   = ps$minDensity,
                    nThreads     = N_CORES,
                    verbose      = FALSE
                )
            )
        }, error = function(e) { status <<- "FAIL"; err_msg <<- conditionMessage(e); NULL })

        if (status == "PASS") {
            n_roh  <- nrow(res$runs)
            n_ind  <- length(unique(res$runs$id))
            elapsed <- .sec(t)
            cat(sprintf("  [%s] %-8s  runs=%7d  ind_with_roh=%4d  time=%.2fs\n",
                        status, ps$label, n_roh, n_ind, elapsed))
            n_pass <- n_pass + 1L
            results[[length(results)+1]] <- list(
                dataset=ds$name, params=ps$label, status="PASS",
                n_runs=n_roh, n_ind=n_ind, time_s=elapsed)
        } else {
            cat(sprintf("  [%s] %-8s  ERROR: %s\n", status, ps$label, err_msg))
            n_fail <- n_fail + 1L
            results[[length(results)+1]] <- list(
                dataset=ds$name, params=ps$label, status="FAIL",
                n_runs=NA, n_ind=NA, time_s=NA)
        }

        # Also test consecutive (quick sanity)
        t2 <- tryCatch({
            system.time(
                res2 <- scanRUNS(
                    genoFile     = ds$bed,
                    method       = "consecutive",
                    minSNP       = ps$minSNP,
                    maxOpp       = ps$maxOpp,
                    maxMiss      = ps$maxMiss,
                    minLengthBps = ps$minLengthBps,
                    maxGap       = ps$maxGap,
                    minDensity   = ps$minDensity,
                    nThreads     = N_CORES,
                    verbose      = FALSE
                )
            )
        }, error = function(e) { cat(sprintf("  [FAIL] %-8s consec  ERROR: %s\n", ps$label, conditionMessage(e))); NULL })
        if (!is.null(t2))
            cat(sprintf("  [PASS] %-8s consec  runs=%7d  time=%.2fs\n",
                        ps$label, nrow(res2$runs), .sec(t2)))
    }
}

cat(sprintf("\n%s  RESULT: %d PASS  |  %d FAIL  (of %d tests)\n",
            .ts(), n_pass, n_fail, n_pass + n_fail))
cat(sprintf("Done at %s\n", format(Sys.time(), "%H:%M:%S")))
