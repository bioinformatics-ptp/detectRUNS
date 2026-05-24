###############################################################################
## detectRUNS — deep validation across 4 species datasets
##
## Datasets (all in Ext_Data/):
##   ADAPTmap  — goat,   4653 animals, 53K SNPs, 144 breeds
##   SELMOL    — cattle, 4095 animals, 44K SNPs, 5 breeds
##   pigData   — pig,    1208 animals, 62K SNPs
##   Innovagen — bovine, 1009 animals, 777K SNPs, 1 breed
##
## For each dataset: ROHom + ROHet × sliding + consecutive × lenient + strict
## = 8 scans per dataset × 4 datasets = 32 total scans.
## Each scan: saveRUNS, summaryRuns, tableRuns, Froh, all plots, markdown report.
## runsIslands run on BED-based scans only.
###############################################################################

suppressPackageStartupMessages({
    library(detectRUNS)
    library(ggplot2)
})

# =============================================================================
# CONFIG
# =============================================================================
EXT_DIR <- "Ext_Data"
RES_DIR <- file.path(EXT_DIR, "results")

DATASETS <- list(
    ADAPTmap = list(
        bed     = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201.bed"),
        ped     = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201_auto.ped"),
        map     = file.path(EXT_DIR, "ADAPTmap_genotypeTOP_20161201_auto.map"),
        species = "goat",
        n_breeds = 144   # many breeds — group plots generated-then-removed
    ),
    SELMOL = list(
        bed     = file.path(EXT_DIR, "SELMOL_codACGT.bed"),
        ped     = file.path(EXT_DIR, "SELMOL_codACGT_auto.ped"),
        map     = file.path(EXT_DIR, "SELMOL_codACGT_auto.map"),
        species = "cattle",
        n_breeds = 5
    ),
    pigData = list(
        bed     = file.path(EXT_DIR, "pigData.bed"),
        ped     = file.path(EXT_DIR, "pigData.ped"),
        map     = file.path(EXT_DIR, "pigData.map"),
        species = "pig",
        n_breeds = NULL   # derived from FAM
    ),
    Innovagen_HD = list(
        bed     = file.path(EXT_DIR, "Innovagen_HD.bed"),
        ped     = NULL,   # BED only — too large to convert
        map     = NULL,
        species = "bovine",
        n_breeds = 1
    )
)

# Parameter sets: lenient (catch short ROH) and strict (confident long ROH)
PARAMS <- list(
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
    )
)

# =============================================================================
# Helpers
# =============================================================================
.section <- function(title) {
    cat("\n", strrep("=", 70), "\n")
    cat(" ", title, "\n")
    cat(strrep("=", 70), "\n")
}

.pdf <- function(path, fn, w = 12, h = 7) {
    grDevices::pdf(path, width = w, height = h)
    tryCatch(fn(), finally = grDevices::dev.off())
    invisible(path)
}

# Save PDF, verify it opened cleanly, optionally delete
.pdf_check <- function(path, fn, w = 12, h = 7, keep = TRUE) {
    ok <- tryCatch({
        .pdf(path, fn, w, h)
        file.exists(path) && file.info(path)$size > 1000
    }, error = function(e) {
        cat("    ERROR:", conditionMessage(e), "\n")
        FALSE
    })
    if (!keep && ok) file.remove(path)
    invisible(ok)
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

    # Many-breeds flag: plots generated then removed after validation
    many_breeds <- !is.null(ds$n_breeds) && ds$n_breeds > 20

    # ------------------------------------------------------------------
    # Scans: ROHom + ROHet × sliding + consecutive × lenient + strict
    # ------------------------------------------------------------------
    scan_results <- list()

    for (type in c("ROHom", "ROHet")) {
        for (method in c("sliding", "consecutive")) {
            for (param_name in c("lenient", "strict")) {
                p   <- PARAMS[[param_name]]
                tag <- paste(type, method, param_name, sep = "_")

                cat(sprintf("\n--- %s | %s | %s | %s ---\n",
                            dsname, type, method, param_name))

                ROHet_flag <- (type == "ROHet")

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

                if (is.null(scan)) next
                print(scan)

                # Save object
                roh_path <- file.path(out, paste0(tag, ".roh"))
                saveRUNS(scan, roh_path)
                cat("  Saved:", basename(roh_path), "\n")

                scan_results[[tag]] <- scan

                # Skip further analysis if no runs detected
                if (nrow(scan$runs) == 0L) {
                    cat("  No runs detected — skipping downstream.\n")
                    next
                }

                scan_out <- file.path(out, tag)
                dir.create(scan_out, showWarnings = FALSE)

                # ------------------------------------------------------
                # summaryRuns
                # ------------------------------------------------------
                cat("  summaryRuns (Class=2)...\n")
                summ2 <- tryCatch(
                    summaryRuns(scan, Class = 2, snpInRuns = TRUE),
                    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
                )
                if (!is.null(summ2)) {
                    write.csv(summ2$summary_ROH_count,
                              file.path(scan_out, "summary_count.csv"), row.names = FALSE)
                    write.csv(summ2$summary_ROH_mean_chr,
                              file.path(scan_out, "summary_mean_chr.csv"), row.names = FALSE)
                    n_grp <- length(unique(summ2$summary_ROH_count$group))
                    cat(sprintf("  summaryRuns OK — %d groups\n", n_grp))
                }

                cat("  summaryRuns (Class=4)...\n")
                summ4 <- tryCatch(
                    summaryRuns(scan, Class = 4, snpInRuns = FALSE),
                    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
                )
                if (!is.null(summ4))
                    write.csv(summ4$result_Froh_class,
                              file.path(scan_out, "Froh_class4.csv"), row.names = FALSE)

                # ------------------------------------------------------
                # tableRuns
                # ------------------------------------------------------
                for (thr in c(0.25, 0.50, 0.75)) {
                    tbl <- tryCatch(
                        tableRuns(scan, threshold = thr),
                        error = function(e) { cat("  tableRuns ERROR:", conditionMessage(e), "\n"); NULL }
                    )
                    if (!is.null(tbl) && nrow(tbl) > 0) {
                        write.csv(tbl, file.path(scan_out,
                                  sprintf("tableRuns_pct%02d.csv", as.integer(thr * 100))),
                                  row.names = FALSE)
                        cat(sprintf("  tableRuns %.0f%%: %d regions\n", thr * 100, nrow(tbl)))
                    }
                }
                cat("  tableRuns OK\n")

                # ------------------------------------------------------
                # Froh
                # ------------------------------------------------------
                cat("  Froh_inbreeding...\n")
                froh_gw <- tryCatch(
                    Froh_inbreeding(scan, genome_wide = TRUE),
                    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
                )
                if (!is.null(froh_gw)) {
                    write.csv(froh_gw, file.path(scan_out, "Froh_genomewide.csv"),
                              row.names = FALSE)
                    cat(sprintf("  Froh genome-wide: mean=%.4f  range=[%.4f, %.4f]\n",
                                mean(froh_gw$Froh_genome),
                                min(froh_gw$Froh_genome),
                                max(froh_gw$Froh_genome)))
                }

                froh_cls <- tryCatch(
                    Froh_inbreedingClass(scan, Class = 2),
                    error = function(e) { cat("  Froh_class ERROR:", conditionMessage(e), "\n"); NULL }
                )
                if (!is.null(froh_cls))
                    write.csv(froh_cls, file.path(scan_out, "Froh_byClass.csv"), row.names = FALSE)

                # ------------------------------------------------------
                # Plots
                # ------------------------------------------------------
                cat("  Generating plots...\n")

                # Violin
                .pdf_check(file.path(scan_out, "violin_sum.pdf"),
                    function() plot_ViolinRuns(scan, method = "sum"),
                    keep = !many_breeds)

                .pdf_check(file.path(scan_out, "violin_mean.pdf"),
                    function() plot_ViolinRuns(scan, method = "mean"),
                    keep = !many_breeds)

                # Distribution
                for (style in c("MeanClass", "MeanChr", "RunsPCT")) {
                    .pdf_check(file.path(scan_out, paste0("dist_", style, ".pdf")),
                        function() plot_DistributionRuns(scan, style = style),
                        keep = !many_breeds)
                }

                # Manhattan
                .pdf_check(file.path(scan_out, "manhattan.pdf"),
                    function() plot_manhattanRuns(scan), w = 16, h = 6,
                    keep = TRUE)   # always keep manhattan

                # SNPs in runs
                .pdf_check(file.path(scan_out, "snps_in_runs.pdf"),
                    function() plot_SnpsInRuns(scan), w = 16, h = 8,
                    keep = !many_breeds)

                # Stacked runs
                .pdf_check(file.path(scan_out, "stacked_runs.pdf"),
                    function() plot_StackedRuns(scan), w = 16, h = 8,
                    keep = !many_breeds)

                # Chromosome runs
                .pdf_check(file.path(scan_out, "plot_runs.pdf"),
                    function() plot_Runs(scan), w = 16, h = 10,
                    keep = !many_breeds)

                # Inbreeding plots (ROHom only)
                if (type == "ROHom") {
                    for (style in c("FrohBoxPlot", "ChrBarPlot", "ChrBoxPlot")) {
                        .pdf_check(file.path(scan_out, paste0("inbreeding_", style, ".pdf")),
                            function() plot_InbreedingChr(scan, style = style),
                            keep = !many_breeds)
                    }
                }

                cat("  Plots done.\n")

                # ------------------------------------------------------
                # runsIslands (BED sliding only — needs snp_freq)
                # ------------------------------------------------------
                if (method == "sliding" && !is.null(scan$snp_freq)) {
                    cat("  runsIslands (n_perm=200)...\n")
                    isl <- tryCatch(
                        runsIslands(scan, n_perm = 200, percentile = 0.99, seed = 42,
                                    verbose = FALSE),
                        error = function(e) { cat("  runsIslands ERROR:", conditionMessage(e), "\n"); NULL }
                    )
                    if (!is.null(isl)) {
                        cat(sprintf("  Islands: %d SNPs across %d chr\n",
                                    nrow(isl$islands),
                                    length(unique(isl$islands$CHR))))
                        print(summary(isl))
                        .pdf_check(file.path(scan_out, "islands.pdf"),
                            function() print(plot(isl)), w = 16, h = 5,
                            keep = TRUE)
                        write.csv(as.data.frame(summary(isl)),
                                  file.path(scan_out, "islands_regions.csv"),
                                  row.names = FALSE)
                    }
                }

                # ------------------------------------------------------
                # Markdown report
                # ------------------------------------------------------
                cat("  reportRUNS (markdown)...\n")
                isl_obj <- if (exists("isl") && !is.null(isl)) isl else NULL
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

                # Clean up isl between iterations
                if (exists("isl")) rm(isl)
            }
        }
    }

    # ------------------------------------------------------------------
    # Cross-scan comparison: lenient vs strict (ROHom sliding)
    # ------------------------------------------------------------------
    s_len  <- scan_results[["ROHom_sliding_lenient"]]
    s_str  <- scan_results[["ROHom_sliding_strict"]]
    if (!is.null(s_len) && !is.null(s_str)) {
        .section(paste(dsname, "— lenient vs strict comparison"))
        cat("ROHom sliding lenient:", nrow(s_len$runs), "runs\n")
        cat("ROHom sliding strict :", nrow(s_str$runs), "runs\n")

        f_len <- Froh_inbreeding(s_len, genome_wide = TRUE)
        f_str <- Froh_inbreeding(s_str, genome_wide = TRUE)
        cat(sprintf("Froh (lenient): mean=%.4f  |  Froh (strict): mean=%.4f\n",
                    mean(f_len$Froh_genome), mean(f_str$Froh_genome)))
    }

    cat(sprintf("\n%s COMPLETE. Output in: %s\n", dsname, out))
}

# =============================================================================
# Done
# =============================================================================
.section("ALL DATASETS DONE")
for (dsname in names(DATASETS)) {
    out <- file.path(RES_DIR, dsname, "detectRUNS")
    files <- list.files(out, recursive = TRUE)
    cat(sprintf("  %s: %d files\n", dsname, length(files)))
}
