###############################################################################
## detectRUNS vs PLINK --homozyg comparison
##
## For each dataset, run PLINK --homozyg and detectRUNS scanRUNS (sliding)
## with identical parameters (lenient set) and compare ROH calls exactly.
##
## Parameters (lenient):
##   minSNP=10 | minLen=100kb | maxOpp=2 | maxMiss=2
##   windowSize=15 | threshold=0.05 | maxGap=1500kb
##
## Comparison: per-run match on (IID, CHR, POS1/from, POS2/to).
## Only autosomes are compared (PLINK default excludes sex chr).
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))

# ---------------------------------------------------------------------------
# Working directory — script lives in dev/, project root is one level up
# ---------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L) {
    project_root <- dirname(dirname(normalizePath(script_path)))
    setwd(project_root)
}
cat("Working directory:", getwd(), "\n")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
PLINK   <- "Ext_Data/plink"
EXT_DIR <- "Ext_Data"
RES_DIR <- file.path(EXT_DIR, "results", "plink_comparison")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

DATASETS <- list(
    ADAPTmap = list(
        bed    = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201.bed"),
        n_auto = 29L   # goat: 29 autosome pairs; chr 30 = X
    ),
    SELMOL = list(
        bed    = file.path(EXT_DIR, "SELMOL_codACGT.bed"),
        n_auto = 29L   # cattle: 29 autosome pairs; chr 30 = X
    ),
    pigData = list(
        bed    = file.path(EXT_DIR, "suini_12_plink.bed"),
        n_auto = 18L   # pig: 18 autosome pairs (all SNPs already filtered 1-18)
    ),
    Innovagen_HD = list(
        bed    = file.path(EXT_DIR, "Innovagen_HD.bed"),
        n_auto = 29L   # bovine: 29 autosome pairs; chr 30 = X, 31 = Y
    )
)

# Lenient parameters
# minDensity = 1/50 SNP/kbp matches PLINK's default --homozyg-density 50 (50 kb/SNP)
P <- list(
    minSNP       = 10,
    maxOpp       = 2,
    maxMiss      = 2,
    minLengthBps = 1e5,
    maxGap       = 1.5e6,
    windowSize   = 15,
    threshold    = 0.05,
    minDensity   = 1/50
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
.section <- function(title) {
    cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
}

.ts <- function() format(Sys.time(), "[%H:%M:%S]")

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
log <- list()

for (dsname in names(DATASETS)) {
    ds      <- DATASETS[[dsname]]
    out_dir <- file.path(RES_DIR, dsname)
    dir.create(out_dir, showWarnings = FALSE)
    prefix  <- tools::file_path_sans_ext(ds$bed)   # BED prefix (no extension)

    .section(sprintf("%s  %s  PLINK vs detectRUNS", .ts(), dsname))

    if (!file.exists(ds$bed)) {
        cat("BED not found — skipping.\n")
        next
    }

    # -----------------------------------------------------------------------
    # 1.  PLINK --homozyg
    # -----------------------------------------------------------------------
    plink_out <- file.path(out_dir, "plink_roh")

    # --homozyg-density N means max N kb per SNP; 1/50 SNP/kbp → 50 kb/SNP → N=50
    plink_density <- round(1 / P$minDensity)   # 50

    plink_cmd <- paste(
        PLINK,
        "--bfile",  prefix,
        "--chr-set", ds$n_auto,
        "--homozyg",
        "--homozyg-snp",              P$minSNP,
        "--homozyg-kb",               P$minLengthBps / 1000,
        "--homozyg-window-snp",       P$windowSize,
        "--homozyg-window-het",       P$maxOpp,
        "--homozyg-window-missing",   P$maxMiss,
        "--homozyg-window-threshold", P$threshold,
        "--homozyg-gap",              P$maxGap / 1000,
        "--homozyg-density",          plink_density,
        "--allow-no-sex",
        "--out", plink_out,
        "--silent"
    )

    cat("Running PLINK...\n")
    t_plink <- system.time(system(plink_cmd))["elapsed"]
    cat(sprintf("  PLINK:  %.1f s\n", t_plink))

    hom_file <- paste0(plink_out, ".hom")
    if (!file.exists(hom_file)) {
        cat("  ERROR: .hom file not produced — check", paste0(plink_out, ".log"), "\n")
        log[[dsname]] <- list(status = "PLINK_ERROR")
        next
    }

    plink_roh <- read.table(hom_file, header = TRUE, stringsAsFactors = FALSE)
    # Restrict to autosomes (PLINK should already do this, but be explicit)
    plink_roh <- plink_roh[plink_roh$CHR %in% seq_len(ds$n_auto), ]
    cat(sprintf("  PLINK:  %d ROH in %d individuals\n",
                nrow(plink_roh), length(unique(plink_roh$IID))))

    # -----------------------------------------------------------------------
    # 2.  detectRUNS
    # -----------------------------------------------------------------------
    cat("Running detectRUNS...\n")
    t_dr <- system.time(
        dr <- tryCatch(
            scanRUNS(genoFile      = ds$bed,
                     method        = "sliding",
                     ROHet         = FALSE,
                     minSNP        = P$minSNP,
                     maxOpp        = P$maxOpp,
                     maxMiss       = P$maxMiss,
                     minLengthBps  = P$minLengthBps,
                     maxGap        = P$maxGap,
                     windowSize    = P$windowSize,
                     threshold     = P$threshold,
                     minDensity    = P$minDensity,
                     verbose       = FALSE),
            error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
        )
    )["elapsed"]
    cat(sprintf("  detectRUNS:  %.1f s\n", t_dr))

    if (is.null(dr)) {
        log[[dsname]] <- list(status = "DR_ERROR")
        next
    }

    # Restrict to autosomes
    dr_runs <- dr$runs
    dr_runs <- dr_runs[dr_runs$chrom %in% seq_len(ds$n_auto), ]
    cat(sprintf("  detectRUNS:  %d ROH in %d individuals\n",
                nrow(dr_runs), length(unique(dr_runs$id))))

    # -----------------------------------------------------------------------
    # 3.  Exact comparison: key = IID_CHR_from_to
    # -----------------------------------------------------------------------
    plink_keys <- with(plink_roh, paste(IID,    CHR,   POS1, POS2, sep = "\t"))
    dr_keys    <- with(dr_runs,   paste(id,  chrom, from,  to,   sep = "\t"))

    n_plink    <- length(plink_keys)
    n_dr       <- length(dr_keys)
    n_match    <- length(intersect(plink_keys, dr_keys))
    only_plink <- setdiff(plink_keys, dr_keys)
    only_dr    <- setdiff(dr_keys,    plink_keys)

    concordant <- length(only_plink) == 0L && length(only_dr) == 0L

    cat(sprintf("\n  PLINK runs:       %d\n", n_plink))
    cat(sprintf("  detectRUNS runs:  %d\n", n_dr))
    cat(sprintf("  Exact matches:    %d\n", n_match))
    cat(sprintf("  Only in PLINK:    %d\n", length(only_plink)))
    cat(sprintf("  Only in detectRUNS: %d\n", length(only_dr)))
    cat(sprintf("  RESULT: %s\n", if (concordant) "*** PASS (100% match) ***" else "!!! FAIL !!!"))

    # Save discrepancy detail tables
    if (length(only_plink) > 0L) {
        disc_p <- plink_roh[plink_keys %in% only_plink,
                            c("IID","CHR","SNP1","SNP2","POS1","POS2","NSNP","KB"), drop=FALSE]
        write.csv(disc_p, file.path(out_dir, "only_in_plink.csv"), row.names = FALSE)
        cat(sprintf("  Saved only_in_plink.csv (%d rows)\n", nrow(disc_p)))
        cat("  First discrepancies (PLINK only):\n")
        print(head(disc_p, 10L), row.names = FALSE)
    }
    if (length(only_dr) > 0L) {
        disc_dr <- dr_runs[dr_keys %in% only_dr, , drop=FALSE]
        write.csv(disc_dr, file.path(out_dir, "only_in_detectRUNS.csv"), row.names = FALSE)
        cat(sprintf("  Saved only_in_detectRUNS.csv (%d rows)\n", nrow(disc_dr)))
        cat("  First discrepancies (detectRUNS only):\n")
        print(head(disc_dr, 10L), row.names = FALSE)
    }

    log[[dsname]] <- list(
        n_plink    = n_plink,
        n_dr       = n_dr,
        n_match    = n_match,
        only_plink = length(only_plink),
        only_dr    = length(only_dr),
        concordant = concordant,
        t_plink    = round(t_plink, 1),
        t_dr       = round(t_dr, 1)
    )
}

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
.section("SUMMARY")
all_pass <- TRUE
for (nm in names(log)) {
    r <- log[[nm]]
    if (!is.null(r$status)) {
        cat(sprintf("  %-15s  [%s]\n", nm, r$status))
        all_pass <- FALSE
        next
    }
    cat(sprintf("  %-15s  PLINK=%d  detectRUNS=%d  match=%d  only_P=%d  only_DR=%d  [%s]\n",
                nm, r$n_plink, r$n_dr, r$n_match,
                r$only_plink, r$only_dr,
                if (r$concordant) "PASS" else "FAIL"))
    if (!r$concordant) all_pass <- FALSE
}
cat(sprintf("\nOverall: %s\n", if (all_pass) "ALL PASS" else "DISCREPANCIES FOUND"))

# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
rpt_path <- file.path(RES_DIR, "plink_comparison_report.md")
rpt <- file(rpt_path, "w")

.w <- function(...) cat(..., "\n", file = rpt, sep = "")

.w("# PLINK vs detectRUNS — ROH Comparison Report")
.w()
.w(sprintf("**Date:** %s  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
.w(sprintf("**detectRUNS:** %s  ", as.character(utils::packageVersion("detectRUNS"))))
.w(sprintf("**PLINK:** v1.9  "))
.w()
.w("## Parameters (lenient)")
.w()
.w(sprintf("| minSNP | minLen | maxOpp | maxMiss | windowSize | threshold | maxGap | minDensity (SNP/kbp) |"))
.w(sprintf("|--------|--------|--------|---------|------------|-----------|--------|---------------------|"))
.w(sprintf("| %d | %s bp | %d | %d | %d | %.2f | %s bp | %.4f (≡ PLINK density %d kb/SNP) |",
    P$minSNP, format(P$minLengthBps, big.mark=","),
    P$maxOpp, P$maxMiss, P$windowSize, P$threshold,
    format(P$maxGap, big.mark=","), P$minDensity, plink_density))
.w()
.w("## Results")
.w()
.w("| Dataset | PLINK ROH | detectRUNS ROH | Exact matches | Only PLINK | Only detectRUNS | PLINK (s) | detectRUNS (s) | Result |")
.w("|---------|-----------|----------------|---------------|------------|-----------------|-----------|----------------|--------|")
for (nm in names(log)) {
    r <- log[[nm]]
    if (!is.null(r$status)) {
        .w(sprintf("| %s | — | — | — | — | — | — | — | **%s** |", nm, r$status))
        next
    }
    .w(sprintf("| %s | %d | %d | %d | %d | %d | %.1f | %.1f | **%s** |",
        nm, r$n_plink, r$n_dr, r$n_match,
        r$only_plink, r$only_dr,
        r$t_plink, r$t_dr,
        if (r$concordant) "PASS" else "FAIL"))
}
.w()
.w(sprintf("**Overall: %s**", if (all_pass) "ALL PASS" else "DISCREPANCIES FOUND — see per-dataset CSV files"))
.w()
.w("## Notes")
.w()
.w("- Comparison restricted to autosomes only (PLINK excludes sex chromosomes by default).")
.w("- Run key: `IID + CHR + POS1(from) + POS2(to)` — exact bp-position match required.")
.w("- Discrepancy files saved in `results/plink_comparison/<dataset>/` when not 100% concordant.")

close(rpt)
cat(sprintf("\nReport: %s\n", rpt_path))
cat(sprintf("Done at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
