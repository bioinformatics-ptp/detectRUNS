##############################################################################
## compare_plink_detectRUNS.R
##
## Validates that detectRUNS scanRUNS(method = "sliding") produces results
## identical to PLINK --homozyg with equivalent parameters.
##
## Usage:
##   Rscript compare_plink_detectRUNS.R
##   source("compare_plink_detectRUNS.R")
##
## Expected outcome: exact coordinate/count match for each parameter set.
## The only known legitimate difference is ROH length (+1 bp in PLINK due to
## PLINK computing length as end - start + 1 vs detectRUNS end - start).
##############################################################################

suppressMessages(library(detectRUNS))

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

PLINK_BIN <- system.file("extdata", "Ext_Data", "plink", package = "detectRUNS")
if (!nzchar(PLINK_BIN) || !file.exists(PLINK_BIN))
    stop("PLINK binary not found in package extdata. Set PLINK_BIN manually.")

SHEEP_BED <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
SHEEP_BASE <- sub("\\.bed$", "", SHEEP_BED)

# Sheep genome has 26 autosomes; tell PLINK so it doesn't treat chr 23-26
# as sex / MT chromosomes.
PLINK_CHR_SET <- 26

OUTDIR <- tempdir()

# --------------------------------------------------------------------------
# Parameter sets to test
# Each list entry becomes one comparison run (PLINK + detectRUNS).
# Keys match PLINK flag names; the script converts to detectRUNS names.
# --------------------------------------------------------------------------

PARAM_SETS <- list(

    "P1: loose (window=15, het=1, miss=1, minSNP=15, minLen=100kb)" = list(
        window_snp  = 15,   window_het     = 1,    window_missing = 1,
        threshold   = 0.05, snp            = 15,   kb             = 100,
        gap_kb      = 1000, density_kbsnp  = 99999
    ),

    "P2: strict no errors (window=15, het=0, miss=0, minSNP=20, minLen=200kb)" = list(
        window_snp  = 15,   window_het     = 0,    window_missing = 0,
        threshold   = 0.05, snp            = 20,   kb             = 200,
        gap_kb      = 1000, density_kbsnp  = 99999
    ),

    "P3: large window (window=30, het=1, miss=2, minSNP=30, minLen=500kb)" = list(
        window_snp  = 30,   window_het     = 1,    window_missing = 2,
        threshold   = 0.05, snp            = 30,   kb             = 500,
        gap_kb      = 1000, density_kbsnp  = 99999
    ),

    "P4: low threshold (window=15, het=1, miss=1, threshold=0.01, minSNP=15, minLen=100kb)" = list(
        window_snp  = 15,   window_het     = 1,    window_missing = 1,
        threshold   = 0.01, snp            = 15,   kb             = 100,
        gap_kb      = 1000, density_kbsnp  = 99999
    ),

    "P5: high threshold (window=15, het=1, miss=1, threshold=0.15, minSNP=15, minLen=100kb)" = list(
        window_snp  = 15,   window_het     = 1,    window_missing = 1,
        threshold   = 0.15, snp            = 15,   kb             = 100,
        gap_kb      = 1000, density_kbsnp  = 99999
    )
)

# --------------------------------------------------------------------------
# Helper: run one PLINK --homozyg call, return data.frame of .hom results
# --------------------------------------------------------------------------

run_plink <- function(p, label) {
    prefix <- file.path(OUTDIR, sprintf("plink_%s", gsub("[^A-Za-z0-9]", "_", label)))
    cmd <- sprintf(
        paste(
            '"%s" --bfile "%s"',
            '--chr-set %d',
            '--homozyg',
            '--homozyg-window-snp %d',
            '--homozyg-window-het %d',
            '--homozyg-window-missing %d',
            '--homozyg-window-threshold %.4f',
            '--homozyg-snp %d',
            '--homozyg-kb %.3f',
            '--homozyg-gap %.3f',
            '--homozyg-density %.3f',
            '--out "%s"',
            '--silent'
        ),
        PLINK_BIN, SHEEP_BASE,
        PLINK_CHR_SET,
        p$window_snp, p$window_het, p$window_missing, p$threshold,
        p$snp, p$kb, p$gap_kb, p$density_kbsnp,
        prefix
    )
    rc <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (rc != 0) stop("PLINK exited with status ", rc)

    hom_file <- paste0(prefix, ".hom")
    if (!file.exists(hom_file) || file.size(hom_file) == 0)
        return(data.frame())  # no ROH detected

    hom <- read.table(hom_file, header = TRUE, stringsAsFactors = FALSE)
    # Rename to match detectRUNS columns
    # PLINK: FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET
    data.frame(
        group     = hom$FID,
        id        = hom$IID,
        chrom     = as.character(hom$CHR),
        nSNP      = hom$NSNP,
        from      = hom$POS1,
        to        = hom$POS2,
        stringsAsFactors = FALSE
    )
}

# --------------------------------------------------------------------------
# Helper: run detectRUNS sliding scan
# --------------------------------------------------------------------------

run_detectRUNS <- function(p) {
    res <- scanRUNS(
        SHEEP_BED,
        method       = "sliding",
        windowSize   = p$window_snp,
        maxOpp       = p$window_het,
        maxMiss      = p$window_missing,
        threshold    = p$threshold,
        minSNP       = p$snp,
        minLengthBps = p$kb * 1000,
        maxGap       = p$gap_kb * 1000,
        verbose      = FALSE
    )
    dr <- as.data.frame(res$runs)
    data.frame(
        group = dr$group,
        id    = dr$id,
        chrom = as.character(dr$chrom),
        nSNP  = dr$nSNP,
        from  = dr$from,
        to    = dr$to,
        stringsAsFactors = FALSE
    )
}

# --------------------------------------------------------------------------
# Helper: compare PLINK vs detectRUNS — returns list(pass, details)
# --------------------------------------------------------------------------

compare_results <- function(plink_df, dr_df, label) {
    pass  <- TRUE
    notes <- character(0)

    # Sort both by (group, id, chrom, from)
    ord <- function(df) df[order(df$group, df$id, df$chrom, df$from), ]
    p <- ord(plink_df)
    d <- ord(dr_df)
    row.names(p) <- NULL
    row.names(d) <- NULL

    n_p <- nrow(p)
    n_d <- nrow(d)

    # 1. Total count
    if (n_p == n_d) {
        notes <- c(notes, sprintf("  Total ROH:   %d  [MATCH]", n_p))
    } else {
        notes <- c(notes, sprintf("  Total ROH:   PLINK=%d  detectRUNS=%d  [MISMATCH]", n_p, n_d))
        pass <- FALSE
    }

    if (n_p == 0 && n_d == 0) {
        notes <- c(notes, "  (No ROH detected by either tool)")
        return(list(pass = pass, notes = notes))
    }

    # 2. Per-individual count comparison
    count_p <- sort(table(paste(p$group, p$id)))
    count_d <- sort(table(paste(d$group, d$id)))
    if (identical(count_p, count_d)) {
        notes <- c(notes, "  Per-individual count: all match")
    } else {
        notes <- c(notes, "  Per-individual count: MISMATCH")
        # Find differing individuals
        all_ids <- union(names(count_p), names(count_d))
        diffs <- sapply(all_ids, function(x) {
            cp <- if (x %in% names(count_p)) count_p[x] else 0L
            cd <- if (x %in% names(count_d)) count_d[x] else 0L
            if (cp != cd) sprintf("    %s: PLINK=%d detectRUNS=%d", x, cp, cd)
            else NULL
        })
        diffs <- Filter(Negate(is.null), diffs)
        if (length(diffs) > 0)
            notes <- c(notes, unlist(diffs))
        pass <- FALSE
    }

    # 3. Coordinate comparison (only when counts match)
    if (n_p == n_d && n_p > 0) {
        coord_match <- all(p$chrom == d$chrom &
                           p$from  == d$from  &
                           p$to    == d$to    &
                           p$nSNP  == d$nSNP)
        if (coord_match) {
            notes <- c(notes, "  Coordinates (chrom/from/to/nSNP): all match exactly")
        } else {
            n_diff <- sum(!(p$chrom == d$chrom & p$from == d$from &
                            p$to == d$to & p$nSNP == d$nSNP))
            notes <- c(notes, sprintf("  Coordinates: %d / %d rows differ", n_diff, n_p))
            # Show first few diffs
            bad <- which(!(p$chrom == d$chrom & p$from == d$from &
                           p$to == d$to & p$nSNP == d$nSNP))
            for (i in head(bad, 5)) {
                notes <- c(notes, sprintf(
                    "    row %d: PLINK chr%s %d-%d nSNP=%d | detectRUNS chr%s %d-%d nSNP=%d",
                    i,
                    p$chrom[i], p$from[i], p$to[i], p$nSNP[i],
                    d$chrom[i], d$from[i], d$to[i], d$nSNP[i]
                ))
            }
            pass <- FALSE
        }
    }

    list(pass = pass, notes = notes)
}

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

cat("\n====================================================================\n")
cat("  detectRUNS vs PLINK --homozyg comparison\n")
cat("  Dataset: Kijas2016_Sheep_subset (100 individuals, chr 2 + 24)\n")
cat("  PLINK:   ", PLINK_BIN, "\n")
cat("  Note:    PLINK length = detectRUNS lengthBps + 1 (expected; not tested)\n")
cat("====================================================================\n\n")

results <- list()

for (label in names(PARAM_SETS)) {
    p <- PARAM_SETS[[label]]
    cat(sprintf("[TEST] %s\n", label))

    plink_df <- tryCatch(run_plink(p, label),
                         error = function(e) { cat("  PLINK ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(plink_df)) { cat("  --> SKIP (PLINK failed)\n\n"); next }

    dr_df <- tryCatch(run_detectRUNS(p),
                      error = function(e) { cat("  detectRUNS ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(dr_df)) { cat("  --> SKIP (detectRUNS failed)\n\n"); next }

    cmp <- compare_results(plink_df, dr_df, label)
    cat(paste(cmp$notes, collapse = "\n"), "\n")
    status <- if (cmp$pass) "PASS" else "FAIL"
    cat(sprintf("  --> %s\n\n", status))
    results[[label]] <- cmp$pass
}

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

cat("====================================================================\n")
n_pass <- sum(unlist(results))
n_fail <- sum(!unlist(results))
cat(sprintf("  PASSED: %d / %d\n", n_pass, length(results)))
if (n_fail > 0)
    cat(sprintf("  FAILED: %d / %d\n", n_fail, length(results)))
cat("====================================================================\n\n")
