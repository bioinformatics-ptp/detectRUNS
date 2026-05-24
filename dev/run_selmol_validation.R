###############################################################################
## SELMOL validation: all scanRUNS() arguments vs PLINK --homozyg
##
## Dataset: SELMOL_codACGT — 44,191 SNPs | 4,095 individuals | 29 bovine autosomes
## Avg inter-SNP gap ~56 kbp.
##
## Tests:
##   Part 1 — BED (C++ engine) vs PLINK: one-at-a-time parameter sweep
##              Parameters: minSNP, maxOpp, maxMiss, minLengthBps, maxGap,
##                          minDensity, windowSize, threshold
##   Part 2 — maxOppRun / maxMissRun (run-level filters, no direct PLINK equiv)
##              Tested as BED_slide internal monotonicity + BED vs PED consistency
##   Part 3 — BED vs PED consistency on a 30-individual subset
##              Key scenarios: ROHom sliding + consecutive
##   Part 4 — ROHet (BED only; no PLINK comparison — PLINK has no equivalent)
##   Part 5 — Multi-parameter combination tests (BED vs PLINK)
##
## Sliding comparison:  exact key match  (IID + CHR + from + to)
## Consecutive comparison: overlap metrics (any overlap, >=50% overlap)
##
## Output: Ext_Data/results/selmol_validation/selmol_validation_report.md
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
cat("Working directory:", getwd(), "\n\n")

# ---------------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------------
PLINK   <- "Ext_Data/plink"

OUT_DIR <- "Ext_Data/results/selmol_validation"
PLK_DIR <- file.path(OUT_DIR, "plink_tmp")
dir.create(PLK_DIR, recursive = TRUE, showWarnings = FALSE)

# 200-individual subset (pre-created; used for all Parts to keep runtime tractable)
# Full SELMOL: 44,191 SNPs | 4,095 individuals | 29 bovine autosomes
BED_PFX <- file.path(OUT_DIR, "selmol_200")   # .bed/.bim/.fam  (200 individuals)
PED     <- file.path(OUT_DIR, "selmol_200.ped")
MAP     <- file.path(OUT_DIR, "selmol_200.map")

if (!file.exists(paste0(BED_PFX, ".bed"))) stop("Subset BED not found — run PLINK --keep first.")

# ---------------------------------------------------------------------------
# Bovine-appropriate baseline parameters
# ---------------------------------------------------------------------------
BASE <- list(
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    windowSize   = 10,
    threshold    = 0.05,
    minLengthBps = 500000,
    maxGap       = 5e6,
    minDensity   = 0
)

# PLINK density: --homozyg-density N kbp/SNP; 0 is invalid → use large number
.plink_density <- function(min_dens)
    if (min_dens <= 0) 10000L else as.integer(round(1 / min_dens))

# ---------------------------------------------------------------------------
# Parameter sweep scenarios
# ---------------------------------------------------------------------------
SWEEP_SHARED <- list(
    list(param="minSNP",       values=c(5, 10, 15, 20, 25, 30)),
    list(param="maxOpp",       values=c(0, 1, 2, 3)),
    list(param="maxMiss",      values=c(0, 1, 2, 3)),
    list(param="minLengthBps", values=c(100000, 250000, 500000, 1000000, 2000000)),
    list(param="maxGap",       values=c(100000, 500000, 1e6, 2e6, 5e6, 10e6)),
    list(param="minDensity",   values=c(0, 1/1000, 1/100, 1/50))
)

SWEEP_SLIDE <- list(
    list(param="windowSize",   values=c(5, 10, 15, 20, 25)),
    list(param="threshold",    values=c(0.05, 0.10, 0.20, 0.50))
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
.section <- function(txt, char="=")
    cat("\n", strrep(char, 72), "\n[", format(Sys.time(),"%H:%M:%S"), "]  ", txt,
        "\n", strrep(char, 72), "\n", sep="")

.make_params <- function(param, value) { p <- BASE; p[[param]] <- value; p }
fmt_pct <- function(x) if (is.na(x)) "  n/a" else sprintf("%5.1f%%", x)

# Run PLINK --homozyg on SELMOL; returns full .hom data.frame (all chroms).
.run_plink <- function(params, label) {
    out_pfx <- file.path(PLK_DIR, paste0("plink_", label))
    density <- .plink_density(params$minDensity)
    cmd <- paste(
        PLINK,
        "--bfile",                    BED_PFX,
        "--chr-set 29",
        "--allow-no-sex --allow-extra-chr",
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
    ret <- system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
    hom <- paste0(out_pfx, ".hom")
    if (ret != 0 || !file.exists(hom)) { warning("PLINK failed: ", label); return(NULL) }
    df <- tryCatch(read.table(hom, header=TRUE, stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df
}

# Run PLINK --homozyg on an arbitrary --bfile prefix (used for subset).
.run_plink_pfx <- function(bed_pfx, params, label, n_chr=29) {
    out_pfx <- file.path(PLK_DIR, paste0("plink_", label))
    density <- .plink_density(params$minDensity)
    cmd <- paste(
        PLINK,
        "--bfile",                    bed_pfx,
        paste0("--chr-set ", n_chr),
        "--allow-no-sex --allow-extra-chr",
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
    ret <- system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
    hom <- paste0(out_pfx, ".hom")
    if (ret != 0 || !file.exists(hom)) { warning("PLINK failed: ", label); return(NULL) }
    df <- tryCatch(read.table(hom, header=TRUE, stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df
}

# Run scanRUNS on BED or PED. Returns runs data.frame or NULL.
.run_dr <- function(method, file_type, params, bed_pfx=BED_PFX,
                    ped_file=PED, map_file=MAP, rohet=FALSE,
                    max_opp_run=NULL, max_miss_run=NULL) {
    tryCatch({
        if (file_type == "BED") {
            res <- scanRUNS(
                genoFile     = paste0(bed_pfx, ".bed"),
                method       = method, ROHet = rohet,
                minSNP       = params$minSNP,       maxOpp  = params$maxOpp,
                maxMiss      = params$maxMiss,       minLengthBps = params$minLengthBps,
                maxGap       = params$maxGap,        windowSize   = params$windowSize,
                threshold    = params$threshold,     minDensity   = params$minDensity,
                maxOppRun    = max_opp_run,          maxMissRun   = max_miss_run,
                nThreads     = 4L, verbose = FALSE)
        } else {
            res <- scanRUNS(
                genoFile     = ped_file, mapFile = map_file,
                method       = method, ROHet = rohet,
                minSNP       = params$minSNP,       maxOpp  = params$maxOpp,
                maxMiss      = params$maxMiss,       minLengthBps = params$minLengthBps,
                maxGap       = params$maxGap,        windowSize   = params$windowSize,
                threshold    = params$threshold,     minDensity   = params$minDensity,
                maxOppRun    = max_opp_run,          maxMissRun   = max_miss_run,
                nThreads     = 1L, verbose = FALSE)
        }
        as.data.frame(res$runs)
    }, error = function(e) { message("  DR error [", method, "/", file_type, "]: ", conditionMessage(e)); NULL })
}

# Exact key match: IID + CHR + from + to.
.compare_exact <- function(plink_df, dr_df) {
    if (is.null(plink_df) || is.null(dr_df))
        return(list(n_p=NA, n_dr=NA, n_exact=NA, n_only_P=NA, n_only_DR=NA,
                    pct_p=NA, pct_dr=NA))
    pk <- if (nrow(plink_df) > 0) with(plink_df, paste(IID, CHR, POS1, POS2, sep="\t")) else character(0)
    dk <- if (nrow(dr_df)    > 0) with(dr_df,    paste(id, chrom, from, to, sep="\t"))   else character(0)
    n_exact <- length(intersect(pk, dk))
    list(n_p=length(pk), n_dr=length(dk), n_exact=n_exact,
         n_only_P=length(setdiff(pk, dk)), n_only_DR=length(setdiff(dk, pk)),
         pct_p =if(length(pk)>0) round(n_exact/length(pk)*100,1) else 100,
         pct_dr=if(length(dk)>0) round(n_exact/length(dk)*100,1) else 100)
}

# Overlap metrics for consecutive vs PLINK.
.compare_overlap <- function(plink_df, dr_df) {
    if (is.null(plink_df) || is.null(dr_df))
        return(list(n_p=NA, n_dr=NA, n_any=NA, n_50=NA, pct_any=NA, pct_50=NA))
    if (nrow(plink_df) == 0)
        return(list(n_p=0, n_dr=if(is.null(dr_df)) NA else nrow(dr_df),
                    n_any=0, n_50=0, pct_any=100, pct_50=100))
    n_any <- 0L; n_50 <- 0L
    for (i in seq_len(nrow(plink_df))) {
        r    <- plink_df[i,]
        cand <- if (!is.null(dr_df) && nrow(dr_df) > 0)
                    dr_df[dr_df$id == r$IID & dr_df$chrom == as.character(r$CHR), ]
                else data.frame()
        if (nrow(cand) == 0) next
        ov   <- pmax(0L, pmin(cand$to, r$POS2) - pmax(cand$from, r$POS1))
        plen <- r$POS2 - r$POS1
        if (any(ov > 0))                          n_any <- n_any + 1L
        if (plen > 0 && any(ov / plen >= 0.50))  n_50  <- n_50  + 1L
    }
    list(n_p=nrow(plink_df), n_dr=if(!is.null(dr_df)) nrow(dr_df) else NA,
         n_any=n_any, n_50=n_50,
         pct_any=round(n_any/nrow(plink_df)*100,1),
         pct_50 =round(n_50 /nrow(plink_df)*100,1))
}

# Accumulate results
results_sweep <- list()

.add_row <- function(param, value, plink_n, bs_ex, bc_ov) {
    results_sweep[[length(results_sweep)+1L]] <<- list(
        param=param, value=format(value, scientific=FALSE), plink_n=plink_n,
        BS_n=bs_ex$n_dr, BS_pct_P=bs_ex$pct_p, BS_pct_DR=bs_ex$pct_dr,
        BS_onlyP=bs_ex$n_only_P, BS_onlyDR=bs_ex$n_only_DR,
        BC_n=bc_ov$n_dr, BC_any=bc_ov$pct_any, BC_50=bc_ov$pct_50)
}

# ===========================================================================
# PART 1 — Parameter sweep: BED_slide + BED_consec vs PLINK
# ===========================================================================
.section("PART 1 — Parameter sweep: BED (C++ engine) vs PLINK")
cat("Dataset: SELMOL subset  (44191 SNPs | 200 individuals | 29 bovine autosomes)\n")
cat("Baseline:", paste(names(BASE), unlist(BASE), sep="=", collapse="  "), "\n\n")

for (sweep in SWEEP_SHARED) {
    cat(sprintf("\n--- Sweeping %s ---\n", sweep$param))
    for (val in sweep$values) {
        params <- .make_params(sweep$param, val)
        label  <- sprintf("%s_%s", sweep$param, gsub("\\.", "p", format(val, scientific=FALSE)))
        cat(sprintf("  %s = %s  ", sweep$param, format(val, scientific=FALSE)))

        plink_df <- .run_plink(params, label)
        n_p <- if (is.null(plink_df)) NA else nrow(plink_df)

        bs <- .run_dr("sliding",     "BED", params)
        bc <- .run_dr("consecutive", "BED", params)

        bs_ex <- .compare_exact(plink_df, bs)
        bc_ov <- .compare_overlap(plink_df, bc)

        cat(sprintf("PLINK=%s  BED_slide=%s (%s vs PLINK)  BED_consec=%s (any=%s >=50%%=%s)\n",
            n_p,
            if(is.null(bs)) "ERR" else nrow(bs), fmt_pct(bs_ex$pct_p),
            if(is.null(bc)) "ERR" else nrow(bc),
            fmt_pct(bc_ov$pct_any), fmt_pct(bc_ov$pct_50)))

        .add_row(sweep$param, val, n_p, bs_ex, bc_ov)
    }
}

for (sweep in SWEEP_SLIDE) {
    cat(sprintf("\n--- Sweeping %s [sliding only] ---\n", sweep$param))
    for (val in sweep$values) {
        params <- .make_params(sweep$param, val)
        label  <- sprintf("%s_%s", sweep$param, gsub("\\.", "p", format(val, scientific=FALSE)))
        cat(sprintf("  %s = %s  ", sweep$param, format(val, scientific=FALSE)))

        plink_df <- .run_plink(params, label)
        n_p <- if (is.null(plink_df)) NA else nrow(plink_df)

        bs    <- .run_dr("sliding", "BED", params)
        bs_ex <- .compare_exact(plink_df, bs)
        bc_ov <- list(n_p=n_p, n_dr=NA, n_any=NA, n_50=NA, pct_any=NA, pct_50=NA)

        cat(sprintf("PLINK=%s  BED_slide=%s (%s vs PLINK)\n",
            n_p,
            if(is.null(bs)) "ERR" else nrow(bs), fmt_pct(bs_ex$pct_p)))

        .add_row(sweep$param, val, n_p, bs_ex, bc_ov)
    }
}

# ===========================================================================
# PART 2 — maxOppRun / maxMissRun (run-level filters)
# ===========================================================================
.section("PART 2 — maxOppRun / maxMissRun (run-level filters)")
cat("No direct PLINK equivalent. Test: BED run count decreases as filter tightens.\n\n")

base_runs <- .run_dr("sliding", "BED", BASE)
n_base <- if (is.null(base_runs)) NA else nrow(base_runs)
cat(sprintf("Baseline (no run-level filter): %s runs\n\n", n_base))

results_runlevel <- list()

cat("--- maxOppRun sweep (maxMissRun=NULL) ---\n")
for (val in c(0, 1, 2, 3)) {
    bs <- .run_dr("sliding", "BED", BASE, max_opp_run=val)
    n  <- if (is.null(bs)) NA else nrow(bs)
    cat(sprintf("  maxOppRun=%d : %s runs\n", val, n))
    results_runlevel[[length(results_runlevel)+1L]] <- list(filter="maxOppRun", value=val, n=n)
}

cat("\n--- maxMissRun sweep (maxOppRun=NULL) ---\n")
for (val in c(0, 1, 2, 3)) {
    bs <- .run_dr("sliding", "BED", BASE, max_miss_run=val)
    n  <- if (is.null(bs)) NA else nrow(bs)
    cat(sprintf("  maxMissRun=%d : %s runs\n", val, n))
    results_runlevel[[length(results_runlevel)+1L]] <- list(filter="maxMissRun", value=val, n=n)
}

cat("\n--- maxOppRun + maxMissRun combined (both=0) ---\n")
bs <- .run_dr("sliding", "BED", BASE, max_opp_run=0, max_miss_run=0)
n  <- if (is.null(bs)) NA else nrow(bs)
cat(sprintf("  maxOppRun=0 + maxMissRun=0 : %s runs\n", n))

cat("\n--- consecutive: maxOppRun / maxMissRun ---\n")
for (val in c(0, 1, 2, 3)) {
    bs <- .run_dr("consecutive", "BED", BASE, max_opp_run=val)
    n  <- if (is.null(bs)) NA else nrow(bs)
    cat(sprintf("  consecutive maxOppRun=%d : %s runs\n", val, n))
}

# ===========================================================================
# PART 3 — BED vs PED consistency (30-individual subset)
# ===========================================================================
.section("PART 3 — BED vs PED consistency (200-individual subset)")

sub_pfx <- BED_PFX
sub_ped  <- PED
sub_map  <- MAP

if (file.exists(paste0(sub_pfx, ".bed")) && file.exists(sub_ped)) {
    cat("Using pre-built 200-individual subset.\n\n")

    results_ped <- list()

    .check_ped_bed <- function(label, method, params, rohet=FALSE,
                               max_opp_run=NULL, max_miss_run=NULL) {
        bed_r <- .run_dr(method, "BED", params,
                         bed_pfx=sub_pfx, rohet=rohet,
                         max_opp_run=max_opp_run, max_miss_run=max_miss_run)
        ped_r <- .run_dr(method, "PED", params,
                         ped_file=sub_ped, map_file=sub_map, rohet=rohet,
                         max_opp_run=max_opp_run, max_miss_run=max_miss_run)
        n_bed <- if (is.null(bed_r)) NA else nrow(bed_r)
        n_ped <- if (is.null(ped_r)) NA else nrow(ped_r)

        match_str <- "n/a"
        if (!is.na(n_bed) && !is.na(n_ped) && n_bed > 0) {
            ex <- .compare_exact(
                data.frame(IID=ped_r$id, CHR=ped_r$chrom, POS1=ped_r$from, POS2=ped_r$to),
                bed_r)
            match_str <- sprintf("%s (BED) == %s (PED) → %s match",
                                 n_bed, n_ped, fmt_pct(ex$pct_p))
        } else {
            match_str <- sprintf("BED=%s  PED=%s", n_bed, n_ped)
        }
        cat(sprintf("  %-40s  %s\n", label, match_str))
        results_ped[[length(results_ped)+1L]] <<- list(label=label, n_bed=n_bed, n_ped=n_ped)
    }

    cat("--- Sliding ROHom ---\n")
    .check_ped_bed("baseline",                   "sliding", BASE)
    .check_ped_bed("strict (minSNP=25)",          "sliding", .make_params("minSNP", 25))
    .check_ped_bed("tight gap (maxGap=100k)",     "sliding", .make_params("maxGap", 100000))
    .check_ped_bed("high threshold (thr=0.20)",   "sliding", .make_params("threshold", 0.20))
    .check_ped_bed("large window (W=20)",         "sliding", .make_params("windowSize", 20))
    .check_ped_bed("maxOppRun=1",                 "sliding", BASE, max_opp_run=1)
    .check_ped_bed("maxMissRun=1",                "sliding", BASE, max_miss_run=1)

    cat("\n--- Consecutive ROHom ---\n")
    .check_ped_bed("baseline",                   "consecutive", BASE)
    .check_ped_bed("strict (minSNP=25)",          "consecutive", .make_params("minSNP", 25))
    .check_ped_bed("tight gap (maxGap=100k)",     "consecutive", .make_params("maxGap", 100000))
    .check_ped_bed("maxOppRun=1",                 "consecutive", BASE, max_opp_run=1)

} else {
    cat("WARNING: subset extraction failed — skipping Part 3\n")
}

# ===========================================================================
# PART 4 — ROHet (BED path only; no PLINK equivalent)
# ===========================================================================
.section("PART 4 — ROHet (BED path, sliding + consecutive)")
cat("Checking ROHet runs: counts should be > 0 for lenient params,\n")
cat("consistent direction (more lenient = more runs).\n\n")

results_rohet <- list()

rohet_scenarios <- list(
    list(label="baseline",            params=BASE),
    list(label="minSNP=5",            params=.make_params("minSNP", 5)),
    list(label="minSNP=25",           params=.make_params("minSNP", 25)),
    list(label="maxOpp=3",            params=.make_params("maxOpp", 3)),
    list(label="minLengthBps=100000", params=.make_params("minLengthBps", 100000)),
    list(label="minLengthBps=2Mbp",   params=.make_params("minLengthBps", 2000000))
)

cat("--- Sliding ROHet ---\n")
for (sc in rohet_scenarios) {
    bs <- .run_dr("sliding", "BED", sc$params, rohet=TRUE)
    n  <- if (is.null(bs)) "ERR" else nrow(bs)
    cat(sprintf("  %-30s : %s runs\n", sc$label, n))
    results_rohet[[length(results_rohet)+1L]] <- list(method="sliding", label=sc$label, n=n)
}

cat("\n--- Consecutive ROHet ---\n")
for (sc in rohet_scenarios) {
    bs <- .run_dr("consecutive", "BED", sc$params, rohet=TRUE)
    n  <- if (is.null(bs)) "ERR" else nrow(bs)
    cat(sprintf("  %-30s : %s runs\n", sc$label, n))
    results_rohet[[length(results_rohet)+1L]] <- list(method="consec", label=sc$label, n=n)
}

# ===========================================================================
# PART 5 — Multi-parameter combination tests (BED vs PLINK)
# ===========================================================================
.section("PART 5 — Multi-parameter combination tests (BED vs PLINK)")

combo_scenarios <- list(
    list(label="strict_all",
         params=modifyList(BASE, list(minSNP=25, maxOpp=0, maxMiss=0,
                                      minLengthBps=1000000, windowSize=15))),
    list(label="lenient_all",
         params=modifyList(BASE, list(minSNP=5,  maxOpp=3, maxMiss=3,
                                      minLengthBps=100000,  windowSize=5))),
    list(label="strict_snp_lenient_opp",
         params=modifyList(BASE, list(minSNP=25, maxOpp=3))),
    list(label="large_window_low_threshold",
         params=modifyList(BASE, list(windowSize=25, threshold=0.10))),
    list(label="tight_gap_strict_length",
         params=modifyList(BASE, list(maxGap=500000, minLengthBps=1000000))),
    list(label="high_density_filter",
         params=modifyList(BASE, list(minDensity=1/100, minLengthBps=500000))),
    list(label="typical_bovine_ROH_analysis",
         params=list(minSNP=20, maxOpp=1, maxMiss=1, windowSize=10,
                     threshold=0.05, minLengthBps=500000, maxGap=1e6, minDensity=0))
)

results_combo <- list()
for (sc in combo_scenarios) {
    cat(sprintf("  %-35s  ", sc$label))
    label <- paste0("combo_", sc$label)
    plink_df <- .run_plink(sc$params, label)
    n_p <- if (is.null(plink_df)) NA else nrow(plink_df)

    bs    <- .run_dr("sliding",     "BED", sc$params)
    bc    <- .run_dr("consecutive", "BED", sc$params)
    bs_ex <- .compare_exact(plink_df, bs)
    bc_ov <- .compare_overlap(plink_df, bc)

    cat(sprintf("PLINK=%s  BED_slide=%s (%s)  BED_consec=%s (any=%s)\n",
        n_p,
        if(is.null(bs)) "ERR" else nrow(bs), fmt_pct(bs_ex$pct_p),
        if(is.null(bc)) "ERR" else nrow(bc), fmt_pct(bc_ov$pct_any)))

    results_combo[[length(results_combo)+1L]] <- list(
        label=sc$label, plink_n=n_p,
        bs_n=bs_ex$n_dr, bs_pct=bs_ex$pct_p,
        bc_n=bc_ov$n_dr, bc_any=bc_ov$pct_any, bc_50=bc_ov$pct_50)
}

# ===========================================================================
# CONSOLE SUMMARY
# ===========================================================================
.section("PART 1 SUMMARY — BED_slide vs PLINK (all params)")
cat(sprintf("  %-14s  %-12s  %7s  %7s  %7s  %7s  %7s  %8s\n",
    "param", "value", "PLINK_n", "BS_n", "BS_%P", "BS_%DR", "BC_n", "BC_any%"))
cat(strrep("-", 80), "\n")
last_param <- ""
for (r in results_sweep) {
    if (r$param != last_param) { cat("\n"); last_param <- r$param }
    cat(sprintf("  %-14s  %-12s  %7s  %7s  %7s  %7s  %7s  %8s\n",
        r$param, r$value, r$plink_n,
        r$BS_n,  fmt_pct(r$BS_pct_P),  fmt_pct(r$BS_pct_DR),
        if(is.na(r$BC_n)) "  n/a" else r$BC_n,
        fmt_pct(r$BC_any)))
}

# ===========================================================================
# MARKDOWN REPORT
# ===========================================================================
.section("Writing report")
rpt_path <- file.path(OUT_DIR, "selmol_validation_report.md")
rpt  <- file(rpt_path, "w")
.w   <- function(...) cat(..., "\n", file=rpt, sep="")

.w("# detectRUNS — SELMOL Validation Report")
.w()
.w(sprintf("**Date:** %s  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
.w(sprintf("**detectRUNS:** %s  ", as.character(utils::packageVersion("detectRUNS"))))
.w()
.w("## Dataset")
.w()
.w("| File | Content |")
.w("|------|---------|")
.w("| `selmol_200.{bed,bim,fam}` | 44,191 SNPs · 200 individuals (subset of SELMOL) · 29 bovine autosomes · avg gap ~56 kbp |")
.w()
.w("## Baseline parameters")
.w()
.w("| minSNP | maxOpp | maxMiss | windowSize | threshold | minLengthBps | maxGap | minDensity |")
.w("|--------|--------|---------|------------|-----------|--------------|--------|------------|")
.w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
    BASE$minSNP, BASE$maxOpp, BASE$maxMiss, BASE$windowSize, BASE$threshold,
    format(BASE$minLengthBps, big.mark=","),
    format(BASE$maxGap, big.mark=",", scientific=FALSE), BASE$minDensity))
.w()

# Part 1 tables
.w("## Part 1 — Parameter sweep: BED_slide + BED_consec vs PLINK")
.w()
.w("**Sliding**: exact key match (IID+CHR+from+to). **Consecutive**: overlap metrics.")
.w()
last_p <- ""
for (r in results_sweep) {
    if (r$param != last_p) {
        .w(sprintf("### `%s`", r$param))
        .w()
        .w("| value | PLINK_n | BS_n | BS_%%P | BS_onlyP | BS_onlyDR | BS_%%DR | BC_n | BC_any%% | BC_50%% |")
        .w("|-------|---------|------|-------|----------|-----------|--------|------|---------|--------|")
        last_p <- r$param
    }
    .w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        r$value, r$plink_n,
        r$BS_n,  fmt_pct(r$BS_pct_P),  r$BS_onlyP, r$BS_onlyDR, fmt_pct(r$BS_pct_DR),
        if(is.na(r$BC_n)) "n/a" else r$BC_n,
        fmt_pct(r$BC_any), fmt_pct(r$BC_50)))
}
.w()

# Part 2
.w("## Part 2 — maxOppRun / maxMissRun")
.w()
.w(sprintf("Baseline (no run-level filter): **%s runs**. Runs should decrease as filter tightens.", n_base))
.w()
.w("| filter | value | n_runs |")
.w("|--------|-------|--------|")
for (r in results_runlevel) .w(sprintf("| %s | %s | %s |", r$filter, r$value, r$n))
.w()

# Part 3
.w("## Part 3 — BED vs PED consistency (30-individual subset)")
.w()
if (exists("results_ped") && length(results_ped) > 0) {
    .w("| scenario | BED_n | PED_n |")
    .w("|----------|-------|-------|")
    for (r in results_ped) .w(sprintf("| %s | %s | %s |", r$label, r$n_bed, r$n_ped))
} else {
    .w("*Subset extraction failed — Part 3 not run.*")
}
.w()

# Part 4
.w("## Part 4 — ROHet (BED path only)")
.w()
.w("| method | scenario | n_runs |")
.w("|--------|----------|--------|")
for (r in results_rohet) .w(sprintf("| %s | %s | %s |", r$method, r$label, r$n))
.w()

# Part 5
.w("## Part 5 — Multi-parameter combination tests")
.w()
.w("| scenario | PLINK_n | BED_slide_n | BS_%%P | BED_consec_n | BC_any%% | BC_50%% |")
.w("|----------|---------|------------|-------|-------------|---------|--------|")
for (r in results_combo) {
    .w(sprintf("| %s | %s | %s | %s | %s | %s | %s |",
        r$label, r$plink_n,
        r$bs_n,  fmt_pct(r$bs_pct),
        r$bc_n,  fmt_pct(r$bc_any), fmt_pct(r$bc_50)))
}
.w()
close(rpt)

cat(sprintf("\nReport written: %s\n", rpt_path))
cat(sprintf("Done at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
