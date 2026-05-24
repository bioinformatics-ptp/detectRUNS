###############################################################################
## Full PLINK vs detectRUNS comparison
## Datasets: suini_12 (pig, 1208 ind, 54K SNPs) | SELMOL (bovine, 4095 ind, 44K SNPs)
## Metrics: timing (PLINK / 1-thread / all-thread), peak RSS, ROH overlap
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))
suppressPackageStartupMessages(library(data.table))

PLINK   <- "Ext_Data/plink"
N_CORES <- parallel::detectCores(logical = FALSE)
RES_DIR <- "Ext_Data/results/plink_comparison_full"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

.ts  <- function() format(Sys.time(), "[%H:%M:%S]")

# ---------------------------------------------------------------------------
# Peak RSS — macOS /usr/bin/time -l
# ---------------------------------------------------------------------------
rss_mb <- function(cmd) {
    out  <- system(paste("/usr/bin/time -l", cmd, "2>&1"), intern = TRUE)
    line <- grep("maximum resident set size", out, value = TRUE)
    if (!length(line)) return(NA_real_)
    round(as.numeric(gsub("^\\s*(\\d+).*", "\\1", line[1])) / 1024^2, 1)
}

dr_rss_mb <- function(bed, params_str, nthreads) {
    tmp <- tempfile(fileext = ".R")
    writeLines(c(
        "suppressMessages(library(detectRUNS))",
        sprintf("invisible(scanRUNS(genoFile='%s', method='sliding', %s, nThreads=%dL, verbose=FALSE))",
                bed, params_str, nthreads)
    ), tmp)
    on.exit(unlink(tmp))
    rss_mb(paste("Rscript", tmp))
}

# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------
datasets <- list(
    list(
        name   = "suini_12",
        label  = "Pig (suini_12)",
        bed    = "Ext_Data/suini_12_plink.bed",
        n_auto = 18L,
        chr_set = 18L
    ),
    list(
        name   = "SELMOL",
        label  = "Bovine (SELMOL)",
        bed    = "Ext_Data/SELMOL_codACGT.bed",
        n_auto = 29L,
        chr_set = 29L
    )
)

# Parameters — match PLINK defaults closely
P <- list(
    minSNP       = 15L,
    maxOpp       = 1L,
    maxMiss      = 1L,
    minLengthBps = 500000L,
    maxGap       = 5000000L,
    windowSize   = 15L,
    threshold    = 0.05,
    minDensity   = 1 / 1000
)
params_str <- sprintf(
    "windowSize=%dL, minSNP=%dL, maxOpp=%dL, maxMiss=%dL, threshold=%.2f, minLengthBps=%dL, maxGap=%dL, minDensity=1/1000",
    P$windowSize, P$minSNP, P$maxOpp, P$maxMiss, P$threshold, P$minLengthBps, P$maxGap
)

cat(sprintf("\n%s  FULL PLINK vs detectRUNS COMPARISON  (%d cores)\n", .ts(), N_CORES))
cat(strrep("=", 72), "\n")
cat("Parameters:", params_str, "\n\n")

all_results <- list()

for (ds in datasets) {
    bim  <- sub("\\.bed$", ".bim", ds$bed)
    base <- sub("\\.bed$", "", ds$bed)
    snps <- nrow(fread(bim, header = FALSE))
    inds <- nrow(fread(sub("\\.bed$", ".fam", ds$bed), header = FALSE))
    out_dir <- file.path(RES_DIR, ds$name)
    dir.create(out_dir, showWarnings = FALSE)

    cat(sprintf("\n%s  %s  (%d SNPs | %d ind)\n", .ts(), ds$label, snps, inds))
    cat(strrep("-", 60), "\n")

    # ---- 1. PLINK ----
    plink_prefix <- file.path(out_dir, "plink_roh")
    plink_cmd <- paste(
        PLINK,
        "--bfile", base,
        sprintf("--chr-set %d", ds$chr_set),
        "--homozyg",
        "--homozyg-snp",              P$minSNP,
        "--homozyg-kb",               P$minLengthBps / 1000,
        "--homozyg-gap",              P$maxGap / 1000,
        "--homozyg-window-snp",       P$windowSize,
        "--homozyg-window-het",       P$maxOpp,
        "--homozyg-window-missing",   P$maxMiss,
        "--homozyg-window-threshold", P$threshold,
        "--homozyg-density",          round(1 / P$minDensity),
        "--out", plink_prefix,
        "--silent"
    )

    cat("  Running PLINK...\n")
    t_plink  <- system.time(system(plink_cmd))["elapsed"]
    rss_plink <- rss_mb(plink_cmd)
    cat(sprintf("  PLINK:       %.2f s  |  peak RSS: %.0f MB\n", t_plink, rss_plink))

    hom <- paste0(plink_prefix, ".hom")
    plink_roh <- if (file.exists(hom)) {
        df <- read.table(hom, header = TRUE, stringsAsFactors = FALSE)
        df[df$CHR %in% seq_len(ds$n_auto), ]
    } else data.frame()
    cat(sprintf("  PLINK runs:  %d  (in %d individuals)\n",
                nrow(plink_roh), length(unique(plink_roh$IID))))

    plink_keys <- if (nrow(plink_roh) > 0)
        with(plink_roh, paste(IID, CHR, POS1, POS2, sep = "\t")) else character(0)

    # ---- 2. detectRUNS 1 thread ----
    cat("  Running detectRUNS (1 thread)...\n")
    t_dr1 <- system.time(
        res1 <- scanRUNS(genoFile = ds$bed, method = "sliding",
            windowSize = P$windowSize, minSNP = P$minSNP,
            maxOpp = P$maxOpp, maxMiss = P$maxMiss,
            threshold = P$threshold, minLengthBps = P$minLengthBps,
            maxGap = P$maxGap, minDensity = P$minDensity,
            nThreads = 1L, verbose = FALSE)
    )["elapsed"]
    rss_dr1 <- dr_rss_mb(ds$bed, params_str, 1L)
    cat(sprintf("  detectRUNS 1-thread:   %.2f s  |  peak RSS: %.0f MB\n", t_dr1, rss_dr1))

    # ---- 3. detectRUNS all threads ----
    cat(sprintf("  Running detectRUNS (%d threads)...\n", N_CORES))
    t_drN <- system.time(
        resN <- scanRUNS(genoFile = ds$bed, method = "sliding",
            windowSize = P$windowSize, minSNP = P$minSNP,
            maxOpp = P$maxOpp, maxMiss = P$maxMiss,
            threshold = P$threshold, minLengthBps = P$minLengthBps,
            maxGap = P$maxGap, minDensity = P$minDensity,
            nThreads = N_CORES, verbose = FALSE)
    )["elapsed"]
    rss_drN <- dr_rss_mb(ds$bed, params_str, N_CORES)
    cat(sprintf("  detectRUNS %2d-thread:  %.2f s  |  peak RSS: %.0f MB\n",
                N_CORES, t_drN, rss_drN))

    # ---- 4. ROH overlap ----
    dr_runs <- resN$runs[resN$runs$chrom %in% seq_len(ds$n_auto), ]
    dr_keys  <- with(dr_runs, paste(id, chrom, from, to, sep = "\t"))

    n_pl     <- length(plink_keys)
    n_dr     <- length(dr_keys)
    n_match  <- length(intersect(plink_keys, dr_keys))
    n_pl_only <- n_pl - n_match
    n_dr_only <- n_dr - n_match
    pct_match <- if (n_pl > 0) round(n_match / n_pl * 100, 2) else 100.0

    cat(sprintf("\n  ROH OVERLAP (detectRUNS %d-thread vs PLINK):\n", N_CORES))
    cat(sprintf("    PLINK runs:           %6d\n", n_pl))
    cat(sprintf("    detectRUNS runs:      %6d\n", n_dr))
    cat(sprintf("    Exact matches:        %6d  (%.2f%% of PLINK)\n", n_match, pct_match))
    cat(sprintf("    Only in PLINK:        %6d\n", n_pl_only))
    cat(sprintf("    Only in detectRUNS:   %6d\n", n_dr_only))

    status <- if (pct_match == 100) "*** PASS (100%) ***" else sprintf("%.2f%% match", pct_match)
    cat(sprintf("    RESULT: %s\n", status))

    # Save discrepancies if any
    if (n_pl_only > 0 || n_dr_only > 0) {
        only_plink <- plink_roh[!(plink_keys %in% dr_keys), ]
        only_dr    <- dr_runs[!(dr_keys %in% plink_keys), ]
        if (nrow(only_plink) > 0)
            fwrite(only_plink, file.path(out_dir, "only_in_plink.csv"))
        if (nrow(only_dr) > 0)
            fwrite(only_dr, file.path(out_dir, "only_in_detectRUNS.csv"))
    }

    all_results[[ds$name]] <- list(
        dataset   = ds$label,
        n_snps = snps, n_ind = inds,
        t_plink   = round(t_plink, 2),
        t_dr1     = round(t_dr1, 2),
        t_drN     = round(t_drN, 2),
        rss_plink = rss_plink,
        rss_dr1   = rss_dr1,
        rss_drN   = rss_drN,
        n_plink   = n_pl,
        n_dr      = n_dr,
        n_match   = n_match,
        pct_match = pct_match
    )
}

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
cat(sprintf("\n\n%s  SUMMARY TABLE\n", .ts()))
cat(strrep("=", 72), "\n")
cat(sprintf("%-20s %6s %6s | %6s %6s %6s | %7s %7s %7s | %6s %6s %7s\n",
            "Dataset", "SNPs", "Ind",
            "PLINK", "DR_1T", sprintf("DR_%dT", N_CORES),
            "RSS_PL", "RSS_1T", sprintf("RSS_%dT", N_CORES),
            "N_PLINK", "N_DR", "Match%"))
cat(strrep("-", 90), "\n")
for (r in all_results) {
    cat(sprintf("%-20s %6d %6d | %6.2f %6.2f %6.2f | %7.0f %7.0f %7.0f | %6d %6d %7.2f%%\n",
                r$dataset, r$n_snps, r$n_ind,
                r$t_plink, r$t_dr1, r$t_drN,
                r$rss_plink, r$rss_dr1, r$rss_drN,
                r$n_plink, r$n_dr, r$pct_match))
}

fwrite(rbindlist(lapply(all_results, as.data.frame)),
       file.path(RES_DIR, "plink_comparison_full_results.csv"))
cat(sprintf("\nResults saved to %s\n", RES_DIR))
cat(sprintf("Done at %s\n", format(Sys.time(), "%H:%M:%S")))
