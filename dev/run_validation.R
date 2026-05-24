###############################################################################
## detectRUNS — overnight validation across 4 species datasets
##
## Datasets (all in Ext_Data/):
##   ADAPTmap    — goat,   4653 animals, 53K SNPs, 144 breeds
##   SELMOL      — cattle, 4095 animals, 44K SNPs, 5 breeds
##   pigData     — pig,    1208 animals, 54K SNPs (suini_12_plink)
##   Innovagen   — bovine, 1009 animals, 777K SNPs, 1 breed
##
## Scans:
##   4 param sets x 2 types (ROHom/ROHet) x 2 methods (sliding/consecutive)
##   x 4 datasets = 64 scans
##
## Per scan: saveRUNS, summaryRuns (4 classes), tableRuns (5 thresholds),
##           Froh (3 classes), all plots, runsIslands (2 percentiles),
##           reportRUNS, sanity checks
##
## Extras:
##   - BED vs PED cross-format comparison (ADAPTmap, SELMOL, pigData)
##   - Master summary CSV across all scans
##   - Per-step timing and peak memory (RSS via ps)
##   - Standalone validation_report.md written at the end
###############################################################################

suppressPackageStartupMessages({
    library(detectRUNS)
    library(ggplot2)
    library(parallel)
})

# =============================================================================
# Fix working directory — script lives in dev/, project root is one level up
# =============================================================================
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L) {
    project_root <- dirname(dirname(normalizePath(script_path)))
    setwd(project_root)
    cat("Working directory set to:", project_root, "\n")
} else {
    cat("Working directory:", getwd(), "\n")
}

# =============================================================================
# CONFIG
# =============================================================================
EXT_DIR <- "Ext_Data"
RES_DIR <- file.path(EXT_DIR, "results")
START_TIME <- Sys.time()

DATASETS <- list(
    ADAPTmap = list(
        bed      = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201.bed"),
        ped      = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201_auto.ped"),
        map      = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201_auto.map"),
        species  = "goat",
        n_breeds = 144
    ),
    SELMOL = list(
        bed      = file.path(EXT_DIR, "SELMOL_codACGT.bed"),
        ped      = file.path(EXT_DIR, "SELMOL_codACGT_auto.ped"),
        map      = file.path(EXT_DIR, "SELMOL_codACGT_auto.map"),
        species  = "cattle",
        n_breeds = 5
    ),
    pigData = list(
        bed      = file.path(EXT_DIR, "suini_12_plink.bed"),
        ped      = file.path(EXT_DIR, "suini_12_plink.ped"),
        map      = file.path(EXT_DIR, "suini_12_plink.map"),
        species  = "pig",
        n_breeds = NULL
    ),
    Innovagen_HD = list(
        bed      = file.path(EXT_DIR, "Innovagen_HD.bed"),
        ped      = NULL,
        map      = NULL,
        species  = "bovine",
        n_breeds = 1
    )
)

# 4 parameter sets: very_lenient → lenient → strict → very_strict
PARAMS <- list(
    very_lenient = list(
        minSNP       = 5,
        maxOpp       = 3,
        maxMiss      = 3,
        minLengthBps = 5e4,
        maxGap       = 2e6,
        windowSize   = 10,
        threshold    = 0.10
    ),
    lenient = list(
        minSNP       = 10,
        maxOpp       = 2,
        maxMiss      = 2,
        minLengthBps = 1e5,
        maxGap       = 1.5e6,
        windowSize   = 15,
        threshold    = 0.05
    ),
    strict = list(
        minSNP       = 20,
        maxOpp       = 1,
        maxMiss      = 1,
        minLengthBps = 5e5,
        maxGap       = 1e6,
        windowSize   = 20,
        threshold    = 0.05
    ),
    very_strict = list(
        minSNP       = 30,
        maxOpp       = 0,
        maxMiss      = 1,
        minLengthBps = 1e6,
        maxGap       = 5e5,
        windowSize   = 25,
        threshold    = 0.05
    )
)

TABLRUNS_THRESHOLDS  <- c(0.10, 0.25, 0.50, 0.75, 0.90)
SUMMARY_CLASSES      <- c(1, 2, 4, 8)
FROH_CLASSES         <- c(1, 2, 4)
ISLANDS_PERCENTILES  <- c(0.95, 0.99)
ISLANDS_NPERMS       <- 200L

# Master log of all scans for end-of-run summary
master_log <- data.frame(
    dataset        = character(),
    type           = character(),
    method         = character(),
    params         = character(),
    n_runs         = integer(),
    n_indiv        = integer(),
    n_groups       = integer(),
    froh_mean      = numeric(),
    froh_min       = numeric(),
    froh_max       = numeric(),
    scan_s         = numeric(),   # scanRUNS() time
    total_s        = numeric(),   # full scan block time
    peak_mem_mb    = numeric(),   # RSS peak during scan block
    status         = character(),
    stringsAsFactors = FALSE
)

# Sanity check log
check_log <- data.frame(
    dataset = character(),
    tag     = character(),
    check   = character(),
    result  = character(),
    detail  = character(),
    stringsAsFactors = FALSE
)

# Step-level timing log (scan, summaryRuns, tableRuns, Froh, plots, islands)
step_log <- data.frame(
    dataset = character(),
    tag     = character(),
    step    = character(),
    elapsed_s = numeric(),
    stringsAsFactors = FALSE
)

# =============================================================================
# Helpers
# =============================================================================
.section <- function(title) {
    cat("\n", strrep("=", 70), "\n")
    cat(" ", title, "\n")
    cat(strrep("=", 70), "\n")
}

.ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Current process RSS in MB (macOS/Linux via ps; NA if unavailable)
.mem_mb <- function() {
    tryCatch({
        pid  <- Sys.getpid()
        info <- system(sprintf("ps -o rss= -p %d", pid), intern = TRUE)
        as.numeric(trimws(info[1L])) / 1024
    }, error = function(e) NA_real_)
}

# Time a block and record the step
.timed <- function(log, dsname, tag, step_name, expr) {
    t0  <- proc.time()["elapsed"]
    val <- tryCatch(expr, error = function(e) { cat("  ERROR [", step_name, "]:", conditionMessage(e), "\n"); NULL })
    dt  <- proc.time()["elapsed"] - t0
    cat(sprintf("  %s  %.1f s\n", step_name, dt))
    log <<- rbind(log, data.frame(dataset=dsname, tag=tag, step=step_name,
                                  elapsed_s=round(dt,1), stringsAsFactors=FALSE))
    invisible(val)
}

.pdf_check <- function(path, fn, w = 12, h = 7, keep = TRUE) {
    ok <- tryCatch({
        grDevices::pdf(path, width = w, height = h)
        tryCatch(fn(), finally = grDevices::dev.off())
        file.exists(path) && file.info(path)$size > 500
    }, error = function(e) {
        try(grDevices::dev.off(), silent = TRUE)
        cat("    PDF ERROR:", conditionMessage(e), "\n")
        FALSE
    })
    if (!keep && ok) file.remove(path)
    invisible(ok)
}

.check <- function(log, dsname, tag, check_name, expr, detail = "") {
    result <- tryCatch({
        if (isTRUE(expr)) "PASS" else "FAIL"
    }, error = function(e) "ERROR")
    if (result != "PASS")
        cat(sprintf("  [%s] %s: %s %s\n", result, check_name,
                    if (nchar(detail) > 0) detail else "", ""))
    rbind(log, data.frame(
        dataset = dsname, tag = tag, check = check_name,
        result = result, detail = as.character(detail),
        stringsAsFactors = FALSE
    ))
}

# =============================================================================
# Main loop
# =============================================================================
for (dsname in names(DATASETS)) {
    ds  <- DATASETS[[dsname]]
    out <- file.path(RES_DIR, dsname, "detectRUNS")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)

    .section(paste("DATASET:", dsname, "|", ds$species))
    cat("BED:", ds$bed, "\n")

    if (!file.exists(ds$bed)) {
        cat("  BED file not found — skipping dataset.\n")
        next
    }

    many_breeds  <- !is.null(ds$n_breeds) && ds$n_breeds > 20
    scan_results <- list()

    for (type in c("ROHom", "ROHet")) {
        for (method in c("sliding", "consecutive")) {
            for (param_name in names(PARAMS)) {
                p        <- PARAMS[[param_name]]
                tag      <- paste(type, method, param_name, sep = "_")
                scan_out <- file.path(out, tag)
                t0       <- proc.time()["elapsed"]

                cat(sprintf("\n%s --- %s | %s | %s | %s ---\n",
                            .ts(), dsname, type, method, param_name))

                ROHet_flag <- (type == "ROHet")
                mem_before <- .mem_mb()
                t_scan     <- proc.time()["elapsed"]

                scan <- tryCatch({
                    if (method == "sliding") {
                        scanRUNS(ds$bed, method = "sliding", ROHet = ROHet_flag,
                                 minSNP = p$minSNP, maxOpp = p$maxOpp,
                                 maxMiss = p$maxMiss, minLengthBps = p$minLengthBps,
                                 maxGap = p$maxGap, windowSize = p$windowSize,
                                 threshold = p$threshold)
                    } else {
                        scanRUNS(ds$bed, method = "consecutive", ROHet = ROHet_flag,
                                 minSNP = p$minSNP, maxOpp = p$maxOpp,
                                 maxMiss = p$maxMiss, minLengthBps = p$minLengthBps,
                                 maxGap = p$maxGap)
                    }
                }, error = function(e) {
                    cat("  SCAN ERROR:", conditionMessage(e), "\n")
                    NULL
                })

                scan_s     <- round(proc.time()["elapsed"] - t_scan, 1)
                mem_after  <- .mem_mb()
                peak_mem   <- if (!is.na(mem_after)) round(mem_after, 1) else NA_real_
                step_log   <- rbind(step_log, data.frame(dataset=dsname, tag=tag,
                                    step="scanRUNS", elapsed_s=scan_s,
                                    stringsAsFactors=FALSE))
                cat(sprintf("  scanRUNS: %.1f s  |  RSS: %.0f MB (delta +%.0f MB)\n",
                            scan_s, if (!is.na(mem_after)) mem_after else 0,
                            if (!is.na(mem_before) && !is.na(mem_after))
                                mem_after - mem_before else 0))

                if (is.null(scan)) {
                    master_log <- rbind(master_log, data.frame(
                        dataset=dsname, type=type, method=method, params=param_name,
                        n_runs=NA, n_indiv=NA, n_groups=NA,
                        froh_mean=NA, froh_min=NA, froh_max=NA,
                        scan_s=scan_s, total_s=scan_s, peak_mem_mb=peak_mem,
                        status="SCAN_ERROR", stringsAsFactors=FALSE))
                    next
                }

                print(scan)
                scan_results[[tag]] <- scan

                # Save RUNS object
                saveRUNS(scan, file.path(out, paste0(tag, ".roh")))
                cat(sprintf("  Saved %s.roh  (scan: %.1f s)\n", tag, scan_s))

                n_runs   <- nrow(scan$runs)
                n_indiv  <- length(unique(scan$runs$id))
                n_groups <- length(unique(scan$runs$group))

                # Sanity checks
                check_log <- .check(check_log, dsname, tag,
                    "runs_nonnegative", n_runs >= 0)
                check_log <- .check(check_log, dsname, tag,
                    "very_lenient_has_runs",
                    !(param_name == "very_lenient" && n_runs == 0),
                    detail = sprintf("n_runs=%d", n_runs))
                check_log <- .check(check_log, dsname, tag,
                    "strict_le_lenient",
                    !(param_name == "very_strict" &&
                      !is.null(scan_results[[sub("very_strict", "very_lenient", tag)]]) &&
                      n_runs > nrow(scan_results[[sub("very_strict","very_lenient",tag)]]$runs)),
                    detail = "very_strict should have fewer runs than very_lenient")

                if (n_runs == 0L) {
                    cat("  No runs detected — skipping downstream.\n")
                    master_log <- rbind(master_log, data.frame(
                        dataset=dsname, type=type, method=method, params=param_name,
                        n_runs=0L, n_indiv=0L, n_groups=0L,
                        froh_mean=0, froh_min=0, froh_max=0,
                        scan_s=scan_s, total_s=scan_s, peak_mem_mb=peak_mem,
                        status="NO_RUNS", stringsAsFactors=FALSE))
                    next
                }

                dir.create(scan_out, showWarnings = FALSE)

                # -------------------------------------------------------
                # summaryRuns — multiple classes
                # -------------------------------------------------------
                cat("  summaryRuns...\n")
                for (cls in SUMMARY_CLASSES) {
                    sm <- tryCatch(
                        summaryRuns(scan, Class = cls, snpInRuns = (cls == 2)),
                        error = function(e) { cat("  summaryRuns ERROR (Class=",cls,"): ",conditionMessage(e),"\n"); NULL }
                    )
                    if (!is.null(sm)) {
                        write.csv(sm$summary_ROH_count,
                                  file.path(scan_out, sprintf("summary_count_class%d.csv", cls)),
                                  row.names = FALSE)
                        write.csv(sm$summary_ROH_mean_chr,
                                  file.path(scan_out, sprintf("summary_mean_chr_class%d.csv", cls)),
                                  row.names = FALSE)
                        write.csv(sm$result_Froh_class,
                                  file.path(scan_out, sprintf("Froh_class%d.csv", cls)),
                                  row.names = FALSE)
                    }
                }
                cat("  summaryRuns OK\n")

                # -------------------------------------------------------
                # tableRuns — multiple thresholds
                # -------------------------------------------------------
                cat("  tableRuns...\n")
                any_table <- FALSE
                for (thr in TABLRUNS_THRESHOLDS) {
                    tbl <- tryCatch(
                        tableRuns(scan, threshold = thr),
                        error = function(e) { cat("  tableRuns ERROR:", conditionMessage(e), "\n"); NULL }
                    )
                    if (!is.null(tbl) && nrow(tbl) > 0) {
                        write.csv(tbl,
                                  file.path(scan_out, sprintf("tableRuns_pct%02d.csv", as.integer(thr*100))),
                                  row.names = FALSE)
                        cat(sprintf("  tableRuns %.0f%%: %d regions\n", thr*100, nrow(tbl)))
                        any_table <- TRUE
                    }
                }
                check_log <- .check(check_log, dsname, tag,
                    "tableRuns_produced_output", any_table || param_name == "very_strict",
                    detail = "no regions at any threshold")
                cat("  tableRuns OK\n")

                # -------------------------------------------------------
                # Froh — genome-wide + multiple classes
                # -------------------------------------------------------
                cat("  Froh...\n")
                froh_gw <- tryCatch(
                    Froh_inbreeding(scan, genome_wide = TRUE),
                    error = function(e) { cat("  Froh ERROR:", conditionMessage(e), "\n"); NULL }
                )
                froh_mean <- froh_min <- froh_max <- NA
                if (!is.null(froh_gw)) {
                    write.csv(froh_gw, file.path(scan_out, "Froh_genomewide.csv"), row.names = FALSE)
                    froh_mean <- mean(froh_gw$Froh_genome)
                    froh_min  <- min(froh_gw$Froh_genome)
                    froh_max  <- max(froh_gw$Froh_genome)
                    cat(sprintf("  Froh: mean=%.4f  min=%.4f  max=%.4f\n",
                                froh_mean, froh_min, froh_max))
                    check_log <- .check(check_log, dsname, tag,
                        "Froh_in_01", all(froh_gw$Froh_genome >= 0 & froh_gw$Froh_genome <= 1),
                        detail = sprintf("range=[%.4f, %.4f]", froh_min, froh_max))
                }
                for (cls in FROH_CLASSES) {
                    fc <- tryCatch(
                        Froh_inbreedingClass(scan, Class = cls),
                        error = function(e) NULL
                    )
                    if (!is.null(fc))
                        write.csv(fc, file.path(scan_out, sprintf("Froh_byClass%d.csv", cls)), row.names = FALSE)
                }
                cat("  Froh OK\n")

                # -------------------------------------------------------
                # Plots
                # -------------------------------------------------------
                cat("  Plots...\n")

                .pdf_check(file.path(scan_out, "violin_sum.pdf"),
                    function() plot_ViolinRuns(scan, method = "sum"),
                    keep = !many_breeds)

                .pdf_check(file.path(scan_out, "violin_mean.pdf"),
                    function() plot_ViolinRuns(scan, method = "mean"),
                    keep = !many_breeds)

                for (style in c("MeanClass", "MeanChr", "RunsPCT")) {
                    .pdf_check(file.path(scan_out, paste0("dist_", style, ".pdf")),
                        function() plot_DistributionRuns(scan, style = style),
                        keep = !many_breeds)
                }

                .pdf_check(file.path(scan_out, "manhattan.pdf"),
                    function() plot_manhattanRuns(scan), w = 16, h = 6,
                    keep = TRUE)

                .pdf_check(file.path(scan_out, "snps_in_runs.pdf"),
                    function() plot_SnpsInRuns(scan), w = 16, h = 8,
                    keep = !many_breeds)

                .pdf_check(file.path(scan_out, "stacked_runs.pdf"),
                    function() plot_StackedRuns(scan), w = 16, h = 8,
                    keep = !many_breeds)

                .pdf_check(file.path(scan_out, "plot_runs.pdf"),
                    function() plot_Runs(scan), w = 16, h = 10,
                    keep = !many_breeds)

                if (type == "ROHom") {
                    for (style in c("FrohBoxPlot", "ChrBarPlot", "ChrBoxPlot")) {
                        local({
                            s <- style
                            .pdf_check(file.path(scan_out, paste0("inbreeding_", s, ".pdf")),
                                function() plot_InbreedingChr(scan, style = s),
                                keep = !many_breeds)
                        })
                    }
                }
                cat("  Plots OK\n")

                # -------------------------------------------------------
                # runsIslands — sliding ROHom lenient only, 2 percentiles
                # -------------------------------------------------------
                isl_list <- list()
                if (method == "sliding" && type == "ROHom" && param_name == "lenient") {
                    for (pct in ISLANDS_PERCENTILES) {
                        pct_tag <- sprintf("p%02d", as.integer(pct * 100))
                        cat(sprintf("  runsIslands (%s, n_perm=%d)...\n", pct_tag, ISLANDS_NPERMS))
                        isl <- tryCatch(
                            runsIslands(scan, bed_path = ds$bed,
                                        n_perm = ISLANDS_NPERMS,
                                        percentile = pct, seed = 42, verbose = FALSE),
                            error = function(e) {
                                cat("  runsIslands ERROR:", conditionMessage(e), "\n"); NULL
                            }
                        )
                        if (!is.null(isl)) {
                            isl_list[[pct_tag]] <- isl
                            n_isl_snps <- nrow(isl$islands)
                            n_isl_chr  <- length(unique(isl$islands$CHR))
                            cat(sprintf("  Islands (%s): %d SNPs across %d chr\n",
                                        pct_tag, n_isl_snps, n_isl_chr))
                            .pdf_check(file.path(scan_out, sprintf("islands_%s.pdf", pct_tag)),
                                function() print(plot(isl)), w = 16, h = 5,
                                keep = TRUE)
                            write.csv(as.data.frame(summary(isl)),
                                      file.path(scan_out, sprintf("islands_%s_regions.csv", pct_tag)),
                                      row.names = FALSE)
                            check_log <- .check(check_log, dsname, tag,
                                sprintf("islands_%s_p99_subset_p95", pct_tag),
                                !(pct == 0.99 && !is.null(isl_list[["p95"]]) &&
                                  n_isl_snps > nrow(isl_list[["p95"]]$islands)),
                                detail = "p99 should have fewer island SNPs than p95")
                        }
                    }
                }

                # -------------------------------------------------------
                # reportRUNS — markdown + HTML for small-breed datasets
                # -------------------------------------------------------
                cat("  reportRUNS...\n")
                isl_obj <- if (length(isl_list) > 0) isl_list[[length(isl_list)]] else NULL
                rep_out <- tryCatch(
                    reportRUNS(scan,
                               islands    = isl_obj,
                               format     = "markdown",
                               output_dir = scan_out,
                               prefix     = paste0(dsname, "_", tag),
                               overwrite  = TRUE),
                    error = function(e) { cat("  reportRUNS ERROR:", conditionMessage(e), "\n"); NULL }
                )
                if (!is.null(rep_out))
                    cat("  Report:", basename(rep_out$report_file), "\n")

                # HTML report for datasets with few breeds (faster to render)
                if (!many_breeds) {
                    tryCatch(
                        reportRUNS(scan,
                                   islands    = isl_obj,
                                   format     = "html",
                                   output_dir = scan_out,
                                   prefix     = paste0(dsname, "_", tag, "_html"),
                                   overwrite  = TRUE),
                        error = function(e) cat("  HTML report skipped:", conditionMessage(e), "\n")
                    )
                }

                # Record in master log
                total_s  <- round(proc.time()["elapsed"] - t0, 1)
                peak_mem <- round(max(c(peak_mem, .mem_mb()), na.rm = TRUE), 1)
                master_log <- rbind(master_log, data.frame(
                    dataset=dsname, type=type, method=method, params=param_name,
                    n_runs=n_runs, n_indiv=n_indiv, n_groups=n_groups,
                    froh_mean=round(froh_mean,4), froh_min=round(froh_min,4),
                    froh_max=round(froh_max,4),
                    scan_s=scan_s, total_s=total_s, peak_mem_mb=peak_mem,
                    status="OK", stringsAsFactors=FALSE))

            } # param_name
        } # method
    } # type

    # -------------------------------------------------------------------
    # Cross-scan comparisons — lenient vs strict (ROHom sliding)
    # -------------------------------------------------------------------
    s_vl <- scan_results[["ROHom_sliding_very_lenient"]]
    s_le <- scan_results[["ROHom_sliding_lenient"]]
    s_st <- scan_results[["ROHom_sliding_strict"]]
    s_vs <- scan_results[["ROHom_sliding_very_strict"]]

    if (!is.null(s_le) && !is.null(s_st)) {
        .section(paste(dsname, "— parameter comparison (ROHom sliding)"))
        for (nm in c("very_lenient","lenient","strict","very_strict")) {
            sc <- scan_results[[paste0("ROHom_sliding_", nm)]]
            if (!is.null(sc)) {
                f <- tryCatch(Froh_inbreeding(sc, genome_wide=TRUE), error=function(e) NULL)
                froh_s <- if (!is.null(f)) sprintf("%.4f", mean(f$Froh_genome)) else "NA"
                cat(sprintf("  %-14s  runs=%7d  Froh=%s\n", nm, nrow(sc$runs), froh_s))
            }
        }
        # Monotonicity check: very_lenient >= lenient >= strict >= very_strict
        counts <- sapply(c("very_lenient","lenient","strict","very_strict"), function(nm) {
            sc <- scan_results[[paste0("ROHom_sliding_", nm)]]
            if (!is.null(sc)) nrow(sc$runs) else NA
        })
        is_mono <- all(diff(na.omit(counts)) <= 0)
        cat(sprintf("  Run count monotonicity (vl>=l>=s>=vs): %s\n",
                    if (is_mono) "PASS" else "FAIL"))
        check_log <- .check(check_log, dsname, "ROHom_sliding",
            "param_monotonicity",
            is_mono,
            detail = paste(counts, collapse=" >= "))
    }

    # -------------------------------------------------------------------
    # Cross-format comparison: BED vs PED (ROHom sliding lenient)
    # -------------------------------------------------------------------
    if (!is.null(ds$ped) && file.exists(ds$ped) && !is.null(ds$map)) {
        .section(paste(dsname, "— BED vs PED cross-format comparison"))
        p <- PARAMS[["lenient"]]
        cat("  Scanning PED (sliding, ROHom, lenient)...\n")
        scan_ped <- tryCatch(
            scanRUNS(ds$ped, mapFile = ds$map, method = "sliding", ROHet = FALSE,
                     minSNP = p$minSNP, maxOpp = p$maxOpp, maxMiss = p$maxMiss,
                     minLengthBps = p$minLengthBps, maxGap = p$maxGap,
                     windowSize = p$windowSize, threshold = p$threshold),
            error = function(e) { cat("  PED scan ERROR:", conditionMessage(e), "\n"); NULL }
        )
        scan_bed <- scan_results[["ROHom_sliding_lenient"]]
        if (!is.null(scan_ped) && !is.null(scan_bed)) {
            n_bed <- nrow(scan_bed$runs)
            n_ped <- nrow(scan_ped$runs)
            ratio <- n_ped / max(n_bed, 1)
            cat(sprintf("  BED runs: %d  |  PED runs: %d  |  ratio=%.3f\n",
                        n_bed, n_ped, ratio))
            check_log <- .check(check_log, dsname, "BED_vs_PED",
                "run_count_within_5pct",
                abs(ratio - 1) < 0.05,
                detail = sprintf("BED=%d PED=%d ratio=%.3f", n_bed, n_ped, ratio))
            saveRUNS(scan_ped, file.path(out, "ROHom_sliding_lenient_PED.roh"))
            cat("  Saved PED scan for comparison.\n")
        }
    }

    cat(sprintf("\n%s COMPLETE. Output in: %s\n", dsname, out))

} # dsname

# =============================================================================
# Final summary + markdown report
# =============================================================================
.section("ALL DATASETS DONE")

dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
total_elapsed <- difftime(Sys.time(), START_TIME, units = "mins")
cat(sprintf("Total elapsed: %.1f minutes\n\n", total_elapsed))

# Write CSVs
write.csv(master_log, file.path(RES_DIR, "master_summary.csv"),    row.names = FALSE)
write.csv(check_log,  file.path(RES_DIR, "sanity_checks.csv"),     row.names = FALSE)
write.csv(step_log,   file.path(RES_DIR, "step_timing.csv"),       row.names = FALSE)

# Sanity check counts
n_pass  <- sum(check_log$result == "PASS",  na.rm = TRUE)
n_fail  <- sum(check_log$result == "FAIL",  na.rm = TRUE)
n_error <- sum(check_log$result == "ERROR", na.rm = TRUE)
cat(sprintf("Sanity checks: %d PASS  |  %d FAIL  |  %d ERROR\n", n_pass, n_fail, n_error))

# -----------------------------------------------------------------------
# Write standalone markdown report
# -----------------------------------------------------------------------
report_path <- file.path(RES_DIR, "validation_report.md")
rpt <- file(report_path, "w")

.w <- function(...) cat(..., "\n", file = rpt, sep = "")

.w("# detectRUNS Validation Report")
.w()
.w(sprintf("**Date:** %s  ", format(START_TIME, "%Y-%m-%d %H:%M:%S")))
.w(sprintf("**Total elapsed:** %.1f minutes  ", as.numeric(total_elapsed)))
.w(sprintf("**detectRUNS version:** %s  ",
           as.character(utils::packageVersion("detectRUNS"))))
.w(sprintf("**R version:** %s  ", R.version$version.string))
.w(sprintf("**Platform:** %s  ", .Platform$OS.type))
.w(sprintf("**Cores available:** %d  ", parallel::detectCores()))
.w()

.w("---")
.w()
.w("## Datasets")
.w()
.w("| Dataset | Species | Animals | BED SNPs | Breeds |")
.w("|---------|---------|---------|----------|--------|")
ds_info <- list(
    ADAPTmap    = c("goat",   "4,653", "53K",  "144"),
    SELMOL      = c("cattle", "4,095", "44K",  "5"),
    pigData     = c("pig",    "1,208", "54K",  "varies"),
    Innovagen_HD = c("bovine","1,009", "777K", "1")
)
for (nm in names(ds_info)) {
    d <- ds_info[[nm]]
    .w(sprintf("| %s | %s | %s | %s | %s |", nm, d[1], d[2], d[3], d[4]))
}
.w()

.w("## Parameter Sets")
.w()
.w("| Set | minSNP | maxOpp | maxMiss | minLen (bp) | maxGap (bp) | windowSize | threshold |")
.w("|-----|--------|--------|---------|-------------|-------------|------------|-----------|")
for (pnm in names(PARAMS)) {
    p <- PARAMS[[pnm]]
    .w(sprintf("| %s | %d | %d | %d | %s | %s | %d | %.2f |",
               pnm, p$minSNP, p$maxOpp, p$maxMiss,
               format(p$minLengthBps, big.mark=",", scientific=FALSE),
               format(p$maxGap,       big.mark=",", scientific=FALSE),
               p$windowSize, p$threshold))
}
.w()

.w("---")
.w()
.w("## Scan Results")
.w()
ok_scans <- master_log[master_log$status == "OK", ]
.w(sprintf("**Total scans attempted:** %d  ", nrow(master_log)))
.w(sprintf("**Successful:** %d  ", nrow(ok_scans)))
.w(sprintf("**Failed / no runs:** %d  ", nrow(master_log) - nrow(ok_scans)))
.w()
.w("| Dataset | Type | Method | Params | Runs | Individuals | Groups | Froh mean | Froh min | Froh max | Scan (s) | Total (s) | Mem (MB) | Status |")
.w("|---------|------|--------|--------|------|-------------|--------|-----------|----------|----------|----------|-----------|----------|--------|")
for (i in seq_len(nrow(master_log))) {
    r <- master_log[i, ]
    .w(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
               r$dataset, r$type, r$method, r$params,
               ifelse(is.na(r$n_runs),    "-", format(r$n_runs,   big.mark=",")),
               ifelse(is.na(r$n_indiv),   "-", format(r$n_indiv,  big.mark=",")),
               ifelse(is.na(r$n_groups),  "-", as.character(r$n_groups)),
               ifelse(is.na(r$froh_mean), "-", sprintf("%.4f", r$froh_mean)),
               ifelse(is.na(r$froh_min),  "-", sprintf("%.4f", r$froh_min)),
               ifelse(is.na(r$froh_max),  "-", sprintf("%.4f", r$froh_max)),
               ifelse(is.na(r$scan_s),    "-", sprintf("%.1f", r$scan_s)),
               ifelse(is.na(r$total_s),   "-", sprintf("%.1f", r$total_s)),
               ifelse(is.na(r$peak_mem_mb),"-",sprintf("%.0f", r$peak_mem_mb)),
               r$status))
}
.w()

.w("---")
.w()
.w("## Step Timing Summary")
.w()
if (nrow(step_log) > 0) {
    step_agg <- aggregate(elapsed_s ~ step, data = step_log, FUN = function(x)
        c(mean=round(mean(x),1), min=round(min(x),1), max=round(max(x),1), total=round(sum(x),1)))
    .w("| Step | Mean (s) | Min (s) | Max (s) | Total (s) |")
    .w("|------|----------|---------|---------|-----------|")
    for (i in seq_len(nrow(step_agg))) {
        v <- step_agg$elapsed_s[i, ]
        .w(sprintf("| %s | %.1f | %.1f | %.1f | %.1f |",
                   step_agg$step[i], v["mean"], v["min"], v["max"], v["total"]))
    }
}
.w()

.w("---")
.w()
.w("## Memory Usage")
.w()
ok_mem <- ok_scans[!is.na(ok_scans$peak_mem_mb), ]
if (nrow(ok_mem) > 0) {
    .w(sprintf("- **Overall peak RSS:** %.0f MB", max(ok_mem$peak_mem_mb)))
    .w(sprintf("- **Mean peak RSS per scan:** %.0f MB", mean(ok_mem$peak_mem_mb)))
    .w(sprintf("- **Min peak RSS per scan:** %.0f MB", min(ok_mem$peak_mem_mb)))
    if (nrow(ok_mem) > 0) {
        heaviest <- ok_mem[which.max(ok_mem$peak_mem_mb), ]
        .w(sprintf("- **Most memory-intensive scan:** %s / %s / %s / %s (%.0f MB)",
                   heaviest$dataset, heaviest$type, heaviest$method,
                   heaviest$params, heaviest$peak_mem_mb))
    }
} else {
    .w("Memory data not available on this platform.")
}
.w()

.w("---")
.w()
.w("## Sanity Checks")
.w()
.w(sprintf("**PASS:** %d  |  **FAIL:** %d  |  **ERROR:** %d",
           n_pass, n_fail, n_error))
.w()
if (n_fail > 0 || n_error > 0) {
    .w("### Failed / Errored Checks")
    .w()
    .w("| Result | Dataset | Tag | Check | Detail |")
    .w("|--------|---------|-----|-------|--------|")
    bad <- check_log[check_log$result != "PASS", ]
    for (i in seq_len(nrow(bad))) {
        .w(sprintf("| %s | %s | %s | %s | %s |",
                   bad$result[i], bad$dataset[i], bad$tag[i],
                   bad$check[i], bad$detail[i]))
    }
    .w()
} else {
    .w("All checks passed.")
    .w()
}

.w("---")
.w()
.w("## Output Files")
.w()
.w("| Dataset | Files produced |")
.w("|---------|----------------|")
for (dsname in names(DATASETS)) {
    out_d  <- file.path(RES_DIR, dsname, "detectRUNS")
    nfiles <- length(list.files(out_d, recursive = TRUE))
    .w(sprintf("| %s | %d |", dsname, nfiles))
}
.w()
.w(sprintf("_Report generated: %s_", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

close(rpt)
cat(sprintf("\nValidation report written: %s\n", report_path))

# Console summary
cat("\nMaster scan summary (n_runs | scan_s | total_s | peak_mem_mb | status):\n")
print(master_log[, c("dataset","type","method","params","n_runs",
                     "scan_s","total_s","peak_mem_mb","status")],
      row.names = FALSE)

cat(sprintf("\nDone at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))