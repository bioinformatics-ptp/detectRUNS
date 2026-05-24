###############################################################################
## detectRUNS — comprehensive function exercise script
##
## Change the paths in the CONFIG section below, then source the whole file.
## The script runs every exported function, saves all plots and tables,
## and generates HTML reports for ROHom, ROHet, and PED-format scans.
###############################################################################

# =============================================================================
# CONFIG — edit these paths before running
# =============================================================================

BED_FILE   <- "Ext_Data/SELMOL_codACGT.bed"   # BED/BIM/FAM prefix
PED_FILE   <- "Ext_Data/pigData.ped"           # PED file
MAP_FILE   <- "Ext_Data/pigData.map"           # MAP file matching PED
OUT_DIR    <- "dev/test_output"                # all results land here

# =============================================================================
# Setup
# =============================================================================
suppressPackageStartupMessages({
    library(detectRUNS)
    library(ggplot2)
})

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", normalizePath(OUT_DIR), "\n")

# Thin wrapper: open PDF, call fn (which prints to device), close
.pdf <- function(name, fn, w = 10, h = 6) {
    path <- file.path(OUT_DIR, paste0(name, ".pdf"))
    grDevices::pdf(path, width = w, height = h)
    tryCatch(fn(), finally = grDevices::dev.off())
    cat("  saved:", basename(path), "\n")
    invisible(path)
}

.section <- function(title) {
    cat("\n", strrep("=", 70), "\n")
    cat(" ", title, "\n")
    cat(strrep("=", 70), "\n")
}

# =============================================================================
# 1. Scan — BED, sliding, ROHom
# =============================================================================
.section("1. scanRUNS | BED | sliding | ROHom")

roh_slide <- scanRUNS(
    BED_FILE,
    method       = "sliding",
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 1e5,
    maxGap       = 1e6,
    windowSize   = 15,
    threshold    = 0.05
)
print(roh_slide)
cat("Total runs:", nrow(roh_slide$runs), "\n")

# =============================================================================
# 2. Scan — BED, consecutive, ROHom
# =============================================================================
.section("2. scanRUNS | BED | consecutive | ROHom")

roh_cons <- scanRUNS(
    BED_FILE,
    method       = "consecutive",
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 1e5,
    maxGap       = 1e6
)
print(roh_cons)
cat("Total runs:", nrow(roh_cons$runs), "\n")

# =============================================================================
# 3. Scan — BED, sliding, ROHet
# =============================================================================
.section("3. scanRUNS | BED | sliding | ROHet")

roh_het <- scanRUNS(
    BED_FILE,
    method       = "sliding",
    ROHet        = TRUE,
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 1e5,
    maxGap       = 1e6,
    windowSize   = 15,
    threshold    = 0.05
)
print(roh_het)
cat("ROHet runs:", nrow(roh_het$runs), "\n")

# =============================================================================
# 4. Scan — PED format
# =============================================================================
roh_ped <- NULL
if (file.exists(PED_FILE)) {
    .section("4. scanRUNS | PED | sliding | ROHom")
    roh_ped <- scanRUNS(
        PED_FILE,
        mapFile      = MAP_FILE,
        method       = "sliding",
        minSNP       = 10,
        maxOpp       = 1,
        maxMiss      = 1,
        minLengthBps = 5e4,
        maxGap       = 1e6,
        windowSize   = 15,
        threshold    = 0.05
    )
    print(roh_ped)
    cat("PED runs:", nrow(roh_ped$runs), "\n")
}

# =============================================================================
# 5. saveRUNS / loadRUNS round-trip
# =============================================================================
.section("5. saveRUNS / loadRUNS")

roh_file <- file.path(OUT_DIR, "roh_slide.roh")
saveRUNS(roh_slide, roh_file)
roh2 <- loadRUNS(roh_file)
cat("Round-trip OK:", isTRUE(all.equal(
    roh_slide$runs, roh2$runs, check.attributes = FALSE)), "\n")
cat("method preserved:", roh2$method, "\n")
cat("type preserved  :", roh2$type, "\n")

# =============================================================================
# 6. summaryRuns
# =============================================================================
.section("6. summaryRuns | ROHom sliding")

summ <- summaryRuns(roh_slide, Class = 2, snpInRuns = TRUE)
cat("summaryRuns elements:", paste(names(summ), collapse = ", "), "\n")
cat("Groups summarised   :", unique(summ$summary_ROH_count$group), "\n")
cat("SNPinRun rows       :", nrow(summ$SNPinRun), "\n")

summ2 <- summaryRuns(roh_slide, Class = 4, snpInRuns = FALSE)
cat("Class=4 Froh_class cols:", ncol(summ2$result_Froh_class), "\n")

summ_het <- summaryRuns(roh_het, Class = 2, snpInRuns = FALSE)
cat("ROHet summaryRuns OK:", !is.null(summ_het), "\n")

# =============================================================================
# 7. tableRuns
# =============================================================================
.section("7. tableRuns — multiple thresholds")

for (thr in c(0.25, 0.50, 0.75)) {
    tbl <- tableRuns(roh_slide, threshold = thr)
    cat(sprintf("  threshold=%.0f%%: %d common ROH regions\n", thr*100, nrow(tbl)))
}

# =============================================================================
# 8. Froh functions
# =============================================================================
.section("8. Froh_inbreeding / Froh_inbreedingClass")

froh_gw  <- Froh_inbreeding(roh_slide, genome_wide = TRUE)
froh_chr <- Froh_inbreeding(roh_slide, genome_wide = FALSE)
cat("Genome-wide Froh: mean =", round(mean(froh_gw$Froh_genome), 4),
    "  range [", round(min(froh_gw$Froh_genome), 4),
    ",", round(max(froh_gw$Froh_genome), 4), "]\n")
cat("Froh_chr rows:", nrow(froh_chr), "\n")

for (cls in c(2, 4, 8)) {
    f <- Froh_inbreedingClass(roh_slide, Class = cls)
    cat(sprintf("  Froh_inbreedingClass(Class=%d): %d rows, %d cols\n",
                cls, nrow(f), ncol(f)))
}

# =============================================================================
# 9. Plots — ROHom sliding
# =============================================================================
.section("9. Plots | ROHom | sliding")

cat("  plot_ViolinRuns (sum)...\n")
.pdf("01_violin_sum", function() plot_ViolinRuns(roh_slide, method = "sum"))

cat("  plot_ViolinRuns (mean)...\n")
.pdf("02_violin_mean", function() plot_ViolinRuns(roh_slide, method = "mean"))

cat("  plot_DistributionRuns (MeanClass)...\n")
.pdf("03_dist_MeanClass",
     function() plot_DistributionRuns(roh_slide, style = "MeanClass"))

cat("  plot_DistributionRuns (MeanChr)...\n")
.pdf("04_dist_MeanChr",
     function() plot_DistributionRuns(roh_slide, style = "MeanChr"))

cat("  plot_DistributionRuns (RunsPCT)...\n")
.pdf("05_dist_RunsPCT",
     function() plot_DistributionRuns(roh_slide, style = "RunsPCT"))

cat("  plot_InbreedingChr (FrohBoxPlot)...\n")
.pdf("06_froh_boxplot",
     function() plot_InbreedingChr(roh_slide, style = "FrohBoxPlot"))

cat("  plot_InbreedingChr (ChrBarPlot)...\n")
.pdf("07_froh_chr_barplot",
     function() plot_InbreedingChr(roh_slide, style = "ChrBarPlot"))

cat("  plot_InbreedingChr (ChrBoxPlot)...\n")
.pdf("08_froh_chr_boxplot",
     function() plot_InbreedingChr(roh_slide, style = "ChrBoxPlot"))

cat("  plot_manhattanRuns...\n")
.pdf("09_manhattan_sliding",
     function() plot_manhattanRuns(roh_slide), w = 14, h = 6)

cat("  plot_SnpsInRuns...\n")
.pdf("10_snps_in_runs",
     function() plot_SnpsInRuns(roh_slide), w = 14, h = 8)

cat("  plot_StackedRuns...\n")
.pdf("11_stacked_runs",
     function() plot_StackedRuns(roh_slide), w = 14, h = 8)

cat("  plot_Runs (all chromosomes)...\n")
.pdf("12_plot_runs_AllChromosomes",
     function() plot_Runs(roh_slide), w = 16, h = 10)

# =============================================================================
# 10. Plots — ROHom consecutive
# =============================================================================
.section("10. Plots | ROHom | consecutive")

cat("  plot_manhattanRuns (consecutive)...\n")
.pdf("13_manhattan_consecutive",
     function() plot_manhattanRuns(roh_cons), w = 14, h = 6)

cat("  plot_ViolinRuns (consecutive)...\n")
.pdf("14_violin_sum_consecutive",
     function() plot_ViolinRuns(roh_cons, method = "sum"))

# =============================================================================
# 11. Plots — ROHet
# =============================================================================
.section("11. Plots | ROHet")

cat("  plot_ViolinRuns (ROHet)...\n")
.pdf("15_rohet_violin_sum",
     function() plot_ViolinRuns(roh_het, method = "sum"))

cat("  plot_DistributionRuns (ROHet, MeanClass)...\n")
.pdf("16_rohet_dist_MeanClass",
     function() plot_DistributionRuns(roh_het, style = "MeanClass"))

cat("  plot_manhattanRuns (ROHet)...\n")
.pdf("17_rohet_manhattan",
     function() plot_manhattanRuns(roh_het), w = 14, h = 6)

# =============================================================================
# 12. runsIslands — two methods and percentiles
# =============================================================================
.section("12. runsIslands | sliding | p99 | n_perm=500")

isl_slide_p99 <- runsIslands(
    roh_slide,
    n_perm     = 500,
    percentile = 0.99,
    seed       = 42
)
print(isl_slide_p99)
cat("Islands (sliding p99):", nrow(isl_slide_p99$islands), "\n")

.section("12b. runsIslands | consecutive | p95 | n_perm=500")

isl_cons_p95 <- runsIslands(
    roh_cons,
    n_perm     = 500,
    percentile = 0.95,
    seed       = 42
)
cat("Islands (consecutive p95):", nrow(isl_cons_p95$islands), "\n")

cat("  summary.RunsIslands...\n")
print(summary(isl_slide_p99))

cat("  plot.RunsIslands (sliding p99)...\n")
.pdf("18_islands_sliding_p99",
     function() print(plot(isl_slide_p99)), w = 16, h = 5)

cat("  plot.RunsIslands (consecutive p95)...\n")
.pdf("19_islands_consecutive_p95",
     function() print(plot(isl_cons_p95)), w = 16, h = 5)

# =============================================================================
# 13. as_RUNS / as.data.frame
# =============================================================================
.section("13. as_RUNS / as.data.frame.RUNS")

df <- as.data.frame(roh_slide)
cat("as.data.frame: rows =", nrow(df), "  cols =", ncol(df), "\n")

roh_rebuilt <- as_RUNS(roh_slide$runs, bedFile = BED_FILE,
                      method = "sliding", type = "ROHom")
cat("as_RUNS S3 class:", inherits(roh_rebuilt, "RUNS"), "\n")
cat("runs identical  :", isTRUE(all.equal(
    roh_slide$runs, roh_rebuilt$runs, check.attributes = FALSE)), "\n")

# =============================================================================
# 14. reportRUNS — ROHom sliding, all formats
# =============================================================================
.section("14. reportRUNS | ROHom | markdown")

out_md <- reportRUNS(
    roh_slide,
    islands    = isl_slide_p99,
    format     = "markdown",
    output_dir = OUT_DIR,
    prefix     = "report_rohom_md",
    overwrite  = TRUE
)
cat("Markdown report:", basename(out_md$report_file), "\n")
cat("Plots generated:", length(out_md$plots), "\n")

.section("14b. reportRUNS | ROHom | html (self-contained)")

out_html <- reportRUNS(
    roh_slide,
    islands    = isl_slide_p99,
    format     = "html",
    output_dir = OUT_DIR,
    prefix     = "report_rohom_html",
    overwrite  = TRUE
)
cat("HTML report  :", basename(out_html$report_file), "\n")
cat("File size    :", round(file.info(out_html$report_file)$size / 1024), "KB\n")
cat("Base64 plots :", length(out_html$plots), "\n")

# =============================================================================
# 15. reportRUNS — ROHet
# =============================================================================
.section("15. reportRUNS | ROHet | html")

out_het_html <- reportRUNS(
    roh_het,
    format     = "html",
    output_dir = OUT_DIR,
    prefix     = "report_rohet_html",
    overwrite  = TRUE
)
cat("ROHet HTML   :", basename(out_het_html$report_file), "\n")

# =============================================================================
# 16. reportRUNS — consecutive
# =============================================================================
.section("16. reportRUNS | consecutive | html")

out_cons_html <- reportRUNS(
    roh_cons,
    islands    = isl_cons_p95,
    format     = "html",
    output_dir = OUT_DIR,
    prefix     = "report_cons_html",
    overwrite  = TRUE
)
cat("Consecutive HTML:", basename(out_cons_html$report_file), "\n")

# =============================================================================
# 17. reportRUNS — PED format (if available)
# =============================================================================
if (!is.null(roh_ped)) {
    .section("17. reportRUNS | PED scan | html")
    out_ped <- reportRUNS(
        roh_ped,
        format     = "html",
        output_dir = OUT_DIR,
        prefix     = "report_ped_html",
        overwrite  = TRUE
    )
    cat("PED HTML report:", basename(out_ped$report_file), "\n")
}

# =============================================================================
# Done
# =============================================================================
.section("DONE")
cat("All results saved in:", normalizePath(OUT_DIR), "\n\n")
all_files <- list.files(OUT_DIR, recursive = FALSE)
cat("Files generated (", length(all_files), "):\n")
cat(paste(" ", all_files), sep = "\n")
