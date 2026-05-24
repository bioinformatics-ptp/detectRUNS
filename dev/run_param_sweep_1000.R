###############################################################################
## Parameter sweep: PLINK --homozyg vs detectRUNS BED engine
## Dataset: SELMOL 1000-individual random subset | 44191 SNPs | 29 autosomes
## Methods: BED_slide (vs PLINK) + BED_consec (count only; different algorithm)
## Two sections:
##   1. One-at-a-time sweep of all 8 parameters
##   2. Extreme combinations (stress-test parameter interactions)
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))
suppressPackageStartupMessages(library(data.table))

PLINK   <- "Ext_Data/plink"
BED     <- "Ext_Data/selmol_1000.bed"
BIM     <- "Ext_Data/selmol_1000.bim"
RES_DIR <- "Ext_Data/results/param_sweep_1000"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

N_CORES <- parallel::detectCores(logical = FALSE)

.ts  <- function() format(Sys.time(), "[%H:%M:%S]")
.sec <- function(t) round(as.numeric(t["elapsed"]), 2)

# ---------------------------------------------------------------------------
# PLINK runner
# ---------------------------------------------------------------------------
run_plink <- function(out_prefix, p) {
    cmd <- paste(
        PLINK,
        "--bfile", sub("\\.bed$", "", BED),
        "--chr-set 29 --homozyg",
        "--homozyg-snp",              p$minSNP,
        "--homozyg-kb",               p$minLengthBps / 1000,
        "--homozyg-gap",              p$maxGap / 1000,
        "--homozyg-window-snp",       p$windowSize,
        "--homozyg-window-het",       p$maxOpp,
        "--homozyg-window-missing",   p$maxMiss,
        "--homozyg-window-threshold", p$threshold,
        "--homozyg-density",          round(1 / p$minDensity),
        "--out", out_prefix,
        "--silent 2>/dev/null"
    )
    t <- system.time(system(cmd))
    hom <- paste0(out_prefix, ".hom")
    if (!file.exists(hom)) return(list(runs = data.frame(), t = .sec(t)))
    df <- read.table(hom, header = TRUE, stringsAsFactors = FALSE)
    list(runs = df, t = .sec(t))
}

# ---------------------------------------------------------------------------
# Key extractors for exact-match comparison
# ---------------------------------------------------------------------------
plink_key <- function(df) {
    if (nrow(df) == 0) return(character(0))
    with(df, paste(IID, CHR, POS1, POS2, sep = "\t"))
}
dr_key <- function(dt) {
    if (nrow(dt) == 0) return(character(0))
    with(dt, paste(id, chrom, from, to, sep = "\t"))
}

# ---------------------------------------------------------------------------
# Run one configuration and return comparison row
# ---------------------------------------------------------------------------
run_one <- function(label, p, plink_out) {
    pl    <- run_plink(plink_out, p)
    pk    <- plink_key(pl$runs)

    t_bs  <- system.time(
        res_bs <- scanRUNS(genoFile = BED, method = "sliding",
            windowSize = p$windowSize, minSNP = p$minSNP,
            maxOpp = p$maxOpp, maxMiss = p$maxMiss,
            threshold = p$threshold, minLengthBps = p$minLengthBps,
            maxGap = p$maxGap, minDensity = p$minDensity,
            nThreads = N_CORES, verbose = FALSE)
    )
    t_bc  <- system.time(
        res_bc <- scanRUNS(genoFile = BED, method = "consecutive",
            minSNP = p$minSNP, maxOpp = p$maxOpp, maxMiss = p$maxMiss,
            minLengthBps = p$minLengthBps, maxGap = p$maxGap,
            minDensity = p$minDensity,
            nThreads = N_CORES, verbose = FALSE)
    )

    n_pl  <- length(pk)
    dk_bs <- dr_key(res_bs$runs)
    n_bs  <- length(dk_bs)
    n_bc  <- nrow(res_bc$runs)
    match <- length(intersect(pk, dk_bs))
    pct   <- if (n_pl > 0) round(match / n_pl * 100, 2) else 100.0
    status <- if (pct == 100) "[OK]" else sprintf("[FAIL %.2f%%]", pct)

    cat(sprintf("  %-30s  PLINK:%7d  slide:%7d  consec:%7d  match:%s  "
                , label, n_pl, n_bs, n_bc, status))
    cat(sprintf("t_pl=%.2fs  t_slide=%.2fs  t_consec=%.2fs\n",
                pl$t, .sec(t_bs), .sec(t_bc)))

    list(section = NA_character_, label = label,
         n_plink = n_pl, n_bed_slide = n_bs, n_bed_consec = n_bc,
         match_slide = match, pct_slide = pct,
         t_plink = pl$t, t_bed_slide = .sec(t_bs), t_bed_consec = .sec(t_bc),
         minSNP = p$minSNP, maxOpp = p$maxOpp, maxMiss = p$maxMiss,
         minLengthBps = p$minLengthBps, maxGap = p$maxGap,
         windowSize = p$windowSize, threshold = p$threshold,
         minDensity = p$minDensity)
}

# ---------------------------------------------------------------------------
# Baseline
# ---------------------------------------------------------------------------
BASE <- list(
    minSNP       = 15L,
    maxOpp       = 1L,
    maxMiss      = 1L,
    minLengthBps = 500000L,
    maxGap       = 5000000L,
    windowSize   = 15L,
    threshold    = 0.05,
    minDensity   = 1 / 1000
)

results <- list()

# ===========================================================================
# SECTION 1: One-at-a-time sweep (all 8 parameters)
# ===========================================================================
cat(sprintf("\n%s  SECTION 1 — One-at-a-time sweep\n", .ts()))
cat(sprintf("    Dataset: SELMOL 1000 ind | 44191 SNPs | 29 chr | %d cores\n", N_CORES))
cat(strrep("=", 90), "\n")

sweeps <- list(
    list(param = "minSNP",       values = c(3, 5, 10, 15, 20, 30)),
    list(param = "maxOpp",       values = c(0, 1, 2, 3)),
    list(param = "maxMiss",      values = c(0, 1, 2, 3)),
    list(param = "maxGap",       values = c(1e5, 5e5, 1e6, 2e6, 5e6, 1e7)),
    list(param = "windowSize",   values = c(5, 10, 15, 20, 30)),
    list(param = "threshold",    values = c(0.05, 0.10, 0.20, 0.30, 0.50)),
    list(param = "minLengthBps", values = c(1e5, 2.5e5, 5e5, 1e6, 2e6)),
    list(param = "minDensity",   values = c(1/50, 1/100, 1/500, 1/1000, 1/5000))
)

for (sw in sweeps) {
    param <- sw$param
    cat(sprintf("\n--- %s ---\n", param))
    for (val in sw$values) {
        p <- BASE
        p[[param]] <- val
        label <- sprintf("%s=%s", param, val)
        plink_out <- file.path(RES_DIR, paste0("s1_", gsub("[^A-Za-z0-9]", "_", label)))
        row <- run_one(label, p, plink_out)
        row$section <- "oat"
        results[[length(results) + 1]] <- row
    }
}

# ===========================================================================
# SECTION 2: Extreme combinations
# ===========================================================================
cat(sprintf("\n\n%s  SECTION 2 — Extreme combinations\n", .ts()))
cat(strrep("=", 90), "\n\n")

extreme_combos <- list(
    list(
        label = "all_lenient",
        desc  = "maximise runs: small windows, many hets, short runs",
        p = list(minSNP=3L, maxOpp=3L, maxMiss=3L, minLengthBps=1e5L,
                 maxGap=1e7L, windowSize=5L, threshold=0.50, minDensity=1/5000)
    ),
    list(
        label = "all_strict",
        desc  = "minimise runs: large windows, no hets, long runs",
        p = list(minSNP=30L, maxOpp=0L, maxMiss=0L, minLengthBps=2e6L,
                 maxGap=5e5L, windowSize=30L, threshold=0.05, minDensity=1/50)
    ),
    list(
        label = "boundary_500k",
        desc  = "minLengthBps exactly at 500000 (was the off-by-one boundary)",
        p = modifyList(BASE, list(minLengthBps=500000L))
    ),
    list(
        label = "boundary_1M",
        desc  = "minLengthBps at 1000000",
        p = modifyList(BASE, list(minLengthBps=1000000L))
    ),
    list(
        label = "high_het_tolerance",
        desc  = "maxOpp=3, maxMiss=3, threshold=0.50 — stress het filtering",
        p = modifyList(BASE, list(maxOpp=3L, maxMiss=3L, threshold=0.50))
    ),
    list(
        label = "tiny_window",
        desc  = "windowSize=5, minSNP=3, threshold=0.20",
        p = modifyList(BASE, list(windowSize=5L, minSNP=3L, threshold=0.20))
    ),
    list(
        label = "large_window",
        desc  = "windowSize=30, minSNP=20, threshold=0.05",
        p = modifyList(BASE, list(windowSize=30L, minSNP=20L, threshold=0.05))
    ),
    list(
        label = "tight_gap",
        desc  = "maxGap=100000 — split runs at small gaps",
        p = modifyList(BASE, list(maxGap=1e5L))
    ),
    list(
        label = "short_dense_runs",
        desc  = "short runs allowed, strict density",
        p = modifyList(BASE, list(minLengthBps=1e5L, minSNP=5L,
                                  minDensity=1/100, windowSize=5L))
    ),
    list(
        label = "long_sparse_runs",
        desc  = "only very long runs, relaxed density",
        p = modifyList(BASE, list(minLengthBps=5e6L, minSNP=30L,
                                  minDensity=1/5000, windowSize=30L))
    )
)

for (combo in extreme_combos) {
    cat(sprintf("  # %s\n", combo$desc))
    plink_out <- file.path(RES_DIR, paste0("s2_", combo$label))
    row <- run_one(combo$label, combo$p, plink_out)
    row$section <- "extreme"
    results[[length(results) + 1]] <- row
}

# ===========================================================================
# Summary
# ===========================================================================
cat(sprintf("\n\n%s  SUMMARY\n", .ts()))
cat(strrep("=", 90), "\n")

res_dt <- rbindlist(lapply(results, as.data.frame))

ok_oat  <- sum(res_dt$pct_slide[res_dt$section == "oat"]     == 100)
tot_oat <- sum(res_dt$section == "oat")
ok_ext  <- sum(res_dt$pct_slide[res_dt$section == "extreme"] == 100)
tot_ext <- sum(res_dt$section == "extreme")

cat(sprintf("OAT sweep  — exact match (BED slide vs PLINK): %d / %d at 100%%\n", ok_oat, tot_oat))
cat(sprintf("Extreme    — exact match (BED slide vs PLINK): %d / %d at 100%%\n", ok_ext, tot_ext))

fail <- res_dt[res_dt$pct_slide < 100, ]
if (nrow(fail) > 0) {
    cat("\nFAILURES:\n")
    print(fail[, c("section", "label", "n_plink", "n_bed_slide", "pct_slide",
                   "minSNP", "maxOpp", "maxMiss", "minLengthBps", "maxGap",
                   "windowSize", "threshold", "minDensity")],
          row.names = FALSE)
} else {
    cat("\nAll combinations passed.\n")
}

fwrite(res_dt, file.path(RES_DIR, "param_sweep_1000_results.csv"))
cat(sprintf("\nResults saved to %s\n", file.path(RES_DIR, "param_sweep_1000_results.csv")))
cat(sprintf("Done at %s\n", format(Sys.time(), "%H:%M:%S")))
