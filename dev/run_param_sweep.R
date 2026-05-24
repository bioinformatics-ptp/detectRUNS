###############################################################################
## Exhaustive parameter sweep: PLINK --homozyg vs detectRUNS (ROHom)
##
## Methods tested for every parameter variation:
##   PLINK  — --homozyg (sliding window, reference)
##   BED_slide  — scanRUNS(bed, method="sliding")     new C++ engine
##   PED_slide  — scanRUNS(ped, method="sliding")     old R slidingRUNS
##   BED_consec — scanRUNS(bed, method="consecutive") new C++ engine
##   PED_consec — scanRUNS(ped, method="consecutive") old R consecutiveRUNS
##
## Test data: detectRUNS/tests/testthat/test.{ped,map,bed,bim,fam}
##   563 SNPs | 20 individuals | chr 24 | avg gap ~74.6 kbp
##
## Sliding comparison:  exact key match  (IID + CHR + from + to)
## Consecutive comparison: overlap metrics vs PLINK (any overlap, >=50% overlap)
##
## Parameters swept one-at-a-time from a lenient baseline:
##   minSNP     : 3, 5, 10, 15, 20
##   maxOpp     : 0, 1, 2, 3
##   maxMiss    : 0, 1, 2, 3
##   minLengthBps: 1000, 25000, 50000, 100000, 250000
##   maxGap     : 100000, 500000, 1000000, 2000000, 5000000
##   minDensity : 0, 1/1000, 1/100, 1/50
##   windowSize  (sliding only): 5, 10, 15, 20, 30
##   threshold   (sliding only): 0.05, 0.10, 0.20, 0.50
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))
options(warn = 1)

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L) {
    project_root <- dirname(dirname(normalizePath(script_path)))
    setwd(project_root)
}
cat("Working directory:", getwd(), "\n")

# ---------------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------------
PLINK    <- "Ext_Data/plink"
TEST_DIR <- "detectRUNS/tests/testthat"
BED      <- file.path(TEST_DIR, "test.bed")
PED      <- file.path(TEST_DIR, "test.ped")
MAP      <- file.path(TEST_DIR, "test.map")
BED_PFX  <- file.path(TEST_DIR, "test")     # prefix for --bfile

OUT_DIR  <- "Ext_Data/results/param_sweep"
PLK_DIR  <- file.path(OUT_DIR, "plink_tmp")
dir.create(PLK_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Baseline parameters (lenient; produces runs on this small dataset)
# ---------------------------------------------------------------------------
BASE <- list(
    minSNP       = 5,
    maxOpp       = 1,
    maxMiss      = 1,
    windowSize   = 10,
    threshold    = 0.05,
    minLengthBps = 50000,
    maxGap       = 2e6,
    minDensity   = 0        # 0 = no density filter
)

# PLINK density flag: --homozyg-density N (N kbp/SNP; 0 is invalid → use 10000)
.plink_density <- function(min_dens) {
    if (min_dens <= 0) return(10000)
    as.integer(round(1 / min_dens))
}

# ---------------------------------------------------------------------------
# Parameter sweep scenarios (one-at-a-time variation)
# shared = applies to both sliding and consecutive
# slide_only = applies only to sliding
# ---------------------------------------------------------------------------
SWEEP_SHARED <- list(
    list(param="minSNP",      values=c(3, 5, 10, 15, 20)),
    list(param="maxOpp",      values=c(0, 1, 2, 3)),
    list(param="maxMiss",     values=c(0, 1, 2, 3)),
    list(param="minLengthBps",values=c(1000, 25000, 50000, 100000, 250000)),
    list(param="maxGap",      values=c(100000, 500000, 1e6, 2e6, 5e6)),
    list(param="minDensity",  values=c(0, 1/1000, 1/100, 1/50))
)

SWEEP_SLIDE <- list(
    list(param="windowSize",  values=c(5, 10, 15, 20, 30)),
    list(param="threshold",   values=c(0.05, 0.10, 0.20, 0.50))
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
.section <- function(txt)
    cat("\n", strrep("=", 72), "\n", txt, "\n", strrep("=", 72), "\n", sep="")

.ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Build a full parameter list by overriding BASE with the varied parameter.
.make_params <- function(param, value) {
    p <- BASE
    p[[param]] <- value
    p
}

# Run PLINK --homozyg; return data.frame of .hom rows (or NULL on error).
.run_plink <- function(params, label) {
    out_pfx <- file.path(PLK_DIR, paste0("plink_", label))
    density <- .plink_density(params$minDensity)

    cmd <- paste(
        PLINK,
        "--bfile",                    BED_PFX,
        "--chr-set 24",
        "--allow-no-sex",
        "--allow-extra-chr",
        "--homozyg",
        "--homozyg-snp",              params$minSNP,
        "--homozyg-kb",               params$minLengthBps / 1000,
        "--homozyg-window-snp",       params$windowSize,
        "--homozyg-window-het",       params$maxOpp,
        "--homozyg-window-missing",   params$maxMiss,
        "--homozyg-window-threshold", params$threshold,
        "--homozyg-gap",              params$maxGap / 1000,
        "--homozyg-density",          density,
        "--out",                      out_pfx,
        "--silent"
    )
    ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    hom <- paste0(out_pfx, ".hom")
    if (ret != 0 || !file.exists(hom)) return(NULL)
    df <- tryCatch(
        read.table(hom, header = TRUE, stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df[df$CHR == 24, ]   # restrict to the one chromosome in this dataset
}

# Run scanRUNS; return cleaned data.frame of runs (or NULL on error).
.run_dr <- function(method, file_type, params) {
    input  <- if (file_type == "BED") BED else PED
    map_arg <- if (file_type == "PED") MAP else NULL
    tryCatch(
        {
            res <- scanRUNS(
                genoFile     = input,
                mapFile      = map_arg,
                method       = method,
                ROHet        = FALSE,
                minSNP       = params$minSNP,
                maxOpp       = params$maxOpp,
                maxMiss      = params$maxMiss,
                minLengthBps = params$minLengthBps,
                maxGap       = params$maxGap,
                windowSize   = params$windowSize,
                threshold    = params$threshold,
                minDensity   = params$minDensity,
                nThreads     = 1L,
                verbose      = FALSE
            )
            as.data.frame(res$runs)
        },
        error = function(e) { message("  DR error: ", conditionMessage(e)); NULL }
    )
}

# Exact-key comparison (sliding vs PLINK).
.compare_exact <- function(plink_df, dr_df) {
    if (is.null(plink_df) || is.null(dr_df))
        return(list(n_p=NA, n_dr=NA, n_exact=NA, n_only_P=NA, n_only_DR=NA,
                    pct_p=NA, pct_dr=NA))
    pk <- if (nrow(plink_df) > 0)
            with(plink_df, paste(IID, CHR, POS1, POS2, sep="\t"))
          else character(0)
    dk <- if (nrow(dr_df) > 0)
            with(dr_df, paste(id, chrom, from, to, sep="\t"))
          else character(0)
    n_exact  <- length(intersect(pk, dk))
    n_only_P <- length(setdiff(pk, dk))
    n_only_DR <- length(setdiff(dk, pk))
    list(
        n_p       = length(pk),
        n_dr      = length(dk),
        n_exact   = n_exact,
        n_only_P  = n_only_P,
        n_only_DR = n_only_DR,
        pct_p     = if (length(pk) > 0) round(n_exact / length(pk) * 100, 1) else 100,
        pct_dr    = if (length(dk) > 0) round(n_exact / length(dk) * 100, 1) else 100
    )
}

# Overlap comparison (consecutive vs PLINK sliding).
# For each PLINK run, checks detectRUNS for runs of same ind+chr that overlap.
.compare_overlap <- function(plink_df, dr_df) {
    if (is.null(plink_df) || is.null(dr_df))
        return(list(n_p=NA, n_dr=NA, n_any=NA, n_50=NA, pct_any=NA, pct_50=NA))
    if (nrow(plink_df) == 0)
        return(list(n_p=0, n_dr=if(is.null(dr_df)) NA else nrow(dr_df),
                    n_any=0, n_50=0, pct_any=100, pct_50=100))
    n_any <- 0L; n_50 <- 0L
    for (i in seq_len(nrow(plink_df))) {
        r    <- plink_df[i, ]
        cand <- if (!is.null(dr_df) && nrow(dr_df) > 0)
                    dr_df[dr_df$id == r$IID & dr_df$chrom == as.character(r$CHR), ]
                else data.frame()
        if (nrow(cand) == 0) next
        ov    <- pmax(0L, pmin(cand$to, r$POS2) - pmax(cand$from, r$POS1))
        plen  <- r$POS2 - r$POS1
        if (any(ov > 0))                           n_any <- n_any + 1L
        if (plen > 0 && any(ov / plen >= 0.50))    n_50  <- n_50  + 1L
    }
    list(
        n_p    = nrow(plink_df),
        n_dr   = if (!is.null(dr_df)) nrow(dr_df) else NA,
        n_any  = n_any,
        n_50   = n_50,
        pct_any = round(n_any / nrow(plink_df) * 100, 1),
        pct_50  = round(n_50  / nrow(plink_df) * 100, 1)
    )
}

# Formatting helper (used in console summary and report)
fmt_pct <- function(x) if (is.na(x)) "  n/a" else sprintf("%5.1f%%", x)

# ---------------------------------------------------------------------------
# Accumulator: list of result rows
# ---------------------------------------------------------------------------
results <- list()

.add_row <- function(param, value, plink_n,
                     bs_ex, ps_ex,         # BED_slide, PED_slide exact comparison
                     bc_ov, pc_ov) {        # BED_consec, PED_consec overlap comparison
    results[[length(results) + 1L]] <<- list(
        param        = param,
        value        = format(value, scientific = FALSE),
        plink_n      = plink_n,
        # sliding
        BS_n         = bs_ex$n_dr,
        BS_exact     = bs_ex$n_exact,
        BS_onlyP     = bs_ex$n_only_P,
        BS_onlyDR    = bs_ex$n_only_DR,
        BS_pct_P     = bs_ex$pct_p,
        BS_pct_DR    = bs_ex$pct_dr,
        PS_n         = ps_ex$n_dr,
        PS_exact     = ps_ex$n_exact,
        PS_onlyP     = ps_ex$n_only_P,
        PS_onlyDR    = ps_ex$n_only_DR,
        PS_pct_P     = ps_ex$pct_p,
        PS_pct_DR    = ps_ex$pct_dr,
        # consecutive
        BC_n         = bc_ov$n_dr,
        BC_any_pct   = bc_ov$pct_any,
        BC_50pct     = bc_ov$pct_50,
        PC_n         = pc_ov$n_dr,
        PC_any_pct   = pc_ov$pct_any,
        PC_50pct     = pc_ov$pct_50
    )
}

# ---------------------------------------------------------------------------
# Main sweep loop
# ---------------------------------------------------------------------------

# ---- Shared parameters (all four methods) ----
for (sweep in SWEEP_SHARED) {
    .section(sprintf("%s  Sweeping %s", .ts(), sweep$param))

    for (val in sweep$values) {
        params <- .make_params(sweep$param, val)
        label  <- sprintf("%s_%s", sweep$param, gsub("\\.", "p", format(val, scientific=FALSE)))
        cat(sprintf("  %s = %s\n", sweep$param, format(val, scientific=FALSE)))

        # PLINK (uses windowed params = BASE windowSize/threshold)
        plink_df <- .run_plink(params, label)
        n_p <- if (is.null(plink_df)) NA else nrow(plink_df)
        cat(sprintf("    PLINK: %s runs\n", n_p))

        # detectRUNS — 4 variants
        bs <- .run_dr("sliding",     "BED", params)
        ps <- .run_dr("sliding",     "PED", params)
        bc <- .run_dr("consecutive", "BED", params)
        pc <- .run_dr("consecutive", "PED", params)

        cat(sprintf("    BED_slide=%d  PED_slide=%d  BED_consec=%d  PED_consec=%d\n",
                    if(is.null(bs)) -1L else nrow(bs),
                    if(is.null(ps)) -1L else nrow(ps),
                    if(is.null(bc)) -1L else nrow(bc),
                    if(is.null(pc)) -1L else nrow(pc)))

        bs_ex <- .compare_exact(plink_df, bs)
        ps_ex <- .compare_exact(plink_df, ps)
        bc_ov <- .compare_overlap(plink_df, bc)
        pc_ov <- .compare_overlap(plink_df, pc)

        cat(sprintf("    BED_slide: %d exact (%.1f%% of PLINK)  |  BED_consec: %.1f%% any-overlap  %.1f%% >=50%%-overlap\n",
                    if(is.na(bs_ex$n_exact)) -1L else bs_ex$n_exact,
                    if(is.na(bs_ex$pct_p))   -1   else bs_ex$pct_p,
                    if(is.na(bc_ov$pct_any)) -1   else bc_ov$pct_any,
                    if(is.na(bc_ov$pct_50))  -1   else bc_ov$pct_50))

        .add_row(sweep$param, val, n_p, bs_ex, ps_ex, bc_ov, pc_ov)
    }
}

# ---- Sliding-only parameters (no consecutive) ----
for (sweep in SWEEP_SLIDE) {
    .section(sprintf("%s  Sweeping %s  [sliding only]", .ts(), sweep$param))

    for (val in sweep$values) {
        params <- .make_params(sweep$param, val)
        label  <- sprintf("%s_%s", sweep$param, gsub("\\.", "p", format(val, scientific=FALSE)))
        cat(sprintf("  %s = %s\n", sweep$param, format(val, scientific=FALSE)))

        plink_df <- .run_plink(params, label)
        n_p <- if (is.null(plink_df)) NA else nrow(plink_df)
        cat(sprintf("    PLINK: %s runs\n", n_p))

        bs <- .run_dr("sliding", "BED", params)
        ps <- .run_dr("sliding", "PED", params)

        cat(sprintf("    BED_slide=%d  PED_slide=%d\n",
                    if(is.null(bs)) -1L else nrow(bs),
                    if(is.null(ps)) -1L else nrow(ps)))

        bs_ex <- .compare_exact(plink_df, bs)
        ps_ex <- .compare_exact(plink_df, ps)
        # No consecutive for windowed-only params — fill with NA row
        bc_ov <- list(n_p=n_p, n_dr=NA, n_any=NA, n_50=NA, pct_any=NA, pct_50=NA)
        pc_ov <- bc_ov

        cat(sprintf("    BED_slide: %d exact (%.1f%% of PLINK)\n",
                    if(is.na(bs_ex$n_exact)) -1L else bs_ex$n_exact,
                    if(is.na(bs_ex$pct_p))   -1   else bs_ex$pct_p))

        .add_row(sweep$param, val, n_p, bs_ex, ps_ex, bc_ov, pc_ov)
    }
}


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
.section("SUMMARY")
cat(sprintf("  %-14s  %-12s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n",
    "param", "value", "PLINK",
    "BS_n", "BS_%P", "BS_%DR",
    "PS_n", "PS_%P", "PS_%DR",
    "BC_n", "BC_any", "BC_50%"))
cat(strrep("-", 120), "\n")

last_param <- ""
for (r in results) {

    if (r$param != last_param) { cat("\n"); last_param <- r$param }
    cat(sprintf("  %-14s  %-12s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n",
        r$param, r$value, r$plink_n,
        r$BS_n,   fmt_pct(r$BS_pct_P),  fmt_pct(r$BS_pct_DR),
        r$PS_n,   fmt_pct(r$PS_pct_P),  fmt_pct(r$PS_pct_DR),
        r$BC_n,   fmt_pct(r$BC_any_pct), fmt_pct(r$BC_50pct)))
}

# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
.section("Writing report")
rpt_path <- file.path(OUT_DIR, "param_sweep_report_plink_compat.md")
rpt  <- file(rpt_path, "w")
.w   <- function(...) cat(..., "\n", file = rpt, sep = "")

.w("# detectRUNS — Exhaustive Parameter Sweep Report")
.w()
.w(sprintf("**Date:** %s  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
.w(sprintf("**detectRUNS:** %s  ", as.character(utils::packageVersion("detectRUNS"))))
.w()
.w("## Test data")
.w()
.w("| File | Content |")
.w("|------|---------|")
.w(sprintf("| `test.{ped,map,bed,bim,fam}` | 563 SNPs · 20 individuals · chr 24 · avg gap ~74.6 kbp |"))
.w()
.w("## Baseline parameters")
.w()
.w("| minSNP | maxOpp | maxMiss | windowSize | threshold | minLengthBps | maxGap | minDensity |")
.w("|--------|--------|---------|------------|-----------|--------------|--------|------------|")
.w(sprintf("| %d | %d | %d | %d | %.2f | %s | %s | %s |",
    BASE$minSNP, BASE$maxOpp, BASE$maxMiss,
    BASE$windowSize, BASE$threshold,
    format(BASE$minLengthBps, big.mark=","),
    format(BASE$maxGap, big.mark=","),
    ifelse(BASE$minDensity == 0, "0 (none)", format(BASE$minDensity, digits=4))))
.w()
.w("## Column legend")
.w()
.w("**Sliding columns** (comparison against PLINK --homozyg, exact key match IID+CHR+from+to):")
.w("- `PLINK_n` — total PLINK ROH")
.w("- `BS_n` / `PS_n` — BED_slide / PED_slide run count")
.w("- `BS_%P` / `PS_%P` — % of PLINK runs exactly matched")
.w("- `BS_%DR` / `PS_%DR` — % of detectRUNS runs exactly matched by PLINK")
.w("- `BS_onlyP` / `PS_onlyP` — runs only in PLINK (not in detectRUNS)")
.w("- `BS_onlyDR` / `PS_onlyDR` — runs only in detectRUNS (not in PLINK)")
.w()
.w("**Consecutive columns** (overlap metrics vs PLINK; exact match not expected):")
.w("- `BC_n` / `PC_n` — BED_consec / PED_consec run count")
.w("- `BC_any` / `PC_any` — % of PLINK runs with ANY positional overlap in detectRUNS")
.w("- `BC_50%` / `PC_50%` — % of PLINK runs covered ≥50% by a detectRUNS run")
.w()
.w("*Note: `BED_slide` = new C++ engine; `PED_slide` = old R slidingRUNS;")
.w("`BED_consec` = new C++ consecutive; `PED_consec` = old R consecutiveRUNS*")
.w()

# Group results by parameter
params_order <- unique(sapply(results, `[[`, "param"))

for (pname in params_order) {
    rows <- results[sapply(results, function(r) r$param == pname)]
    .w(sprintf("## `%s`", pname))
    .w()

    # Sliding table
    .w("### Sliding window vs PLINK")
    .w()
    .w("| value | PLINK_n | BS_n | BS_%P | BS_onlyP | BS_onlyDR | BS_%DR | PS_n | PS_%P | PS_onlyP | PS_onlyDR | PS_%DR |")
    .w("|-------|---------|------|-------|----------|-----------|--------|------|-------|----------|-----------|--------|")
    for (r in rows) {
        .w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            r$value, r$plink_n,
            r$BS_n,  fmt_pct(r$BS_pct_P),  r$BS_onlyP,  r$BS_onlyDR,  fmt_pct(r$BS_pct_DR),
            r$PS_n,  fmt_pct(r$PS_pct_P),  r$PS_onlyP,  r$PS_onlyDR,  fmt_pct(r$PS_pct_DR)))
    }
    .w()

    # Consecutive table (skip if all NA — windowed-only sweep)
    has_consec <- any(sapply(rows, function(r) !is.na(r$BC_any_pct)))
    if (has_consec) {
        .w("### Consecutive vs PLINK (overlap metrics)")
        .w()
        .w("| value | PLINK_n | BC_n | BC_any% | BC_50% | PC_n | PC_any% | PC_50% |")
        .w("|-------|---------|------|---------|--------|------|---------|--------|")
        for (r in rows) {
            .w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                r$value, r$plink_n,
                r$BC_n, fmt_pct(r$BC_any_pct), fmt_pct(r$BC_50pct),
                r$PC_n, fmt_pct(r$PC_any_pct), fmt_pct(r$PC_50pct)))
        }
        .w()
    }
}

# Overall summary table at the end
.w("## Overall concordance overview")
.w()
.w("*(Only sliding method; exact match % of PLINK runs)*")
.w()
.w("| param | value | PLINK_n | BED_slide_%P | PED_slide_%P | BED_onlyP | BED_onlyDR | PED_onlyP | PED_onlyDR |")
.w("|-------|-------|---------|-------------|-------------|-----------|------------|-----------|------------|")
for (r in results) {
    if (!is.na(r$BS_pct_P))
        .w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            r$param, r$value, r$plink_n,
            fmt_pct(r$BS_pct_P), fmt_pct(r$PS_pct_P),
            r$BS_onlyP, r$BS_onlyDR,
            r$PS_onlyP, r$PS_onlyDR))
}

close(rpt)
cat(sprintf("\nReport written: %s\n", rpt_path))
cat(sprintf("Done at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
