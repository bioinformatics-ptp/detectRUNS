###############################################################################
## detectRUNS — comprehensive function exercise script
##
## Change the paths in the CONFIG section below, then source the whole file.
## The script runs every exported function with both methods and both ROH types,
## saves all plots and tables, and generates a full HTML report.
###############################################################################

# =============================================================================
# CONFIG — edit these paths before running
# =============================================================================

BED_FILE   <- "Ext_Data/SELMOL_codACGT.bed"    # BED/BIM/FAM (no extension needed)
PED_FILE   <- "Ext_Data/pigData.ped"            # PED + MAP
MAP_FILE   <- "Ext_Data/pigData.map"            # MAP for PED format
OUT_DIR    <- "dev/test_output"                 # all results go here
N_CORES    <- parallel::detectCores() - 1L      # threads for BED scans

# =============================================================================
# Setup
# =============================================================================
suppressPackageStartupMessages({
    library(detectRUNS)
    library(ggplot2)
})

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", normalizePath(OUT_DIR), "\n")

.save_pdf <- function(p, name, w = 10, h = 6) {
    path <- file.path(OUT_DIR, paste0(name, ".pdf"))
    if (inherits(p, "gg") || inherits(p, "ggplot")) {
        ggplot2::ggsave(path, plot = p, width = w, height = h)
    } else {
        grDevices::pdf(path, width = w, height = h)
        tryCatch(print(p), finally = grDevices::dev.off())
    }
    cat("  saved:", path, "\n")
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

roh_bed_slide <- scanRUNS(
    BED_FILE,
    method        = "sliding",
    minSNP        = 15,
    maxOpp        = 1,
    maxMiss       = 1,
    minLengthBps  = 1e5,
    maxGap        = 1e6,
    windowSize    = 15,
    threshold     = 0.05,
    nCores        = N_CORES
)
print(roh_bed_slide)
cat("Runs:", nrow(roh_bed_slide$runs), "\n")

# =============================================================================
# 2. Scan — BED, consecutive, ROHom
# =============================================================================
.section("2. scanRUNS | BED | consecutive | ROHom")

roh_bed_cons <- scanRUNS(
    BED_FILE,
    method        = "consecutive",
    minSNP        = 15,
    maxOpp        = 1,
    maxMiss       = 1,
    minLengthBps  = 1e5,
    maxGap        = 1e6,
    nCores        = N_CORES
)
print(roh_bed_cons)
cat("Runs:", nrow(roh_bed_cons$runs), "\n")

# =============================================================================
# 3. Scan — BED, sliding, ROHet
# =============================================================================
.section("3. scanRUNS | BED | sliding | ROHet")

roh_bed_rohet <- scanRUNS(
    BED_FILE,
    method        = "sliding",
    ROHet         = TRUE,
    minSNP        = 15,
    maxOpp        = 1,
    maxMiss       = 1,
    minLengthBps  = 1e5,
    maxGap        = 1e6,
    windowSize    = 15,
    threshold     = 0.05,
    nCores        = N_CORES
)
print(roh_bed_rohet)
cat("ROHet runs:", nrow(roh_bed_rohet$runs), "\n")

# =============================================================================
# 4. Scan — PED format (if file available)
# =============================================================================
if (file.exists(PED_FILE)) {
    .section("4. scanRUNS | PED | sliding | ROHom")
    roh_ped_slide <- scanRUNS(
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
    print(roh_ped_slide)
    cat("PED runs:", nrow(roh_ped_slide$runs), "\n")
} else {
    roh_ped_slide <- NULL
    cat("Skipping PED scan: file not found\n")
}

# =============================================================================
# 5. saveROH / loadROH round-trip
# =============================================================================
.section("5. saveROH / loadROH")

roh_file <- file.path(OUT_DIR, "roh_bed_slide.roh")
saveROH(roh_bed_slide, roh_file)
roh_reloaded <- loadROH(roh_file)
cat("Round-trip OK:", isTRUE(all.equal(
    roh_bed_slide$runs, roh_reloaded$runs, check.attributes = FALSE)), "\n")

# =============================================================================
# 6. summaryRuns
# =============================================================================
.section("6. summaryRuns | ROHom")

summ <- summaryRuns(roh_bed_slide, Class = 2, snpInRuns = TRUE)
cat("summaryRuns elements:", paste(names(summ), collapse = ", "), "\n")
cat("\nMean run length per group (Mb):\n")
print(summ$summary_ROH_mean[, c("group", "mean_run_Mbps")])

# with snpInRuns on ROHet
summ_het <- summaryRuns(roh_bed_rohet, Class = 2, snpInRuns = FALSE)
cat("\nROHet summaryRuns OK:", !is.null(summ_het), "\n")

# =============================================================================
# 7. tableRuns
# =============================================================================
.section("7. tableRuns")

tbl <- tableRuns(roh_bed_slide, threshold = 0.50)
cat("Common ROH regions (>= 50%):", nrow(tbl), "\n")
if (nrow(tbl) > 0) print(head(tbl, 5))

tbl25 <- tableRuns(roh_bed_slide, threshold = 0.25)
cat("Common ROH regions (>= 25%):", nrow(tbl25), "\n")

# =============================================================================
# 8. Froh functions
# =============================================================================
.section("8. Froh_inbreeding / Froh_inbreedingClass")

froh <- Froh_inbreeding(roh_bed_slide)
cat("Froh columns:", paste(names(froh), collapse = ", "), "\n")
cat("Froh range: [", round(min(froh$Froh_genome), 4),
    ",", round(max(froh$Froh_genome), 4), "]\n")

froh_cls <- Froh_inbreedingClass(roh_bed_slide, Class = 2)
cat("Froh by class columns:", paste(names(froh_cls), collapse = ", "), "\n")

# =============================================================================
# 9. Plots — ROHom sliding
# =============================================================================
.section("9. Plots | ROHom | sliding")

cat("  violin (sum)...\n")
p <- plot_DistRuns(roh_bed_slide, style = "violin", plotType = "sum")
.save_pdf(p, "01_violin_sum")

cat("  violin (mean)...\n")
p <- plot_DistRuns(roh_bed_slide, style = "violin", plotType = "mean")
.save_pdf(p, "02_violin_mean")

cat("  histogram (RunsPCT)...\n")
p <- plot_DistRuns(roh_bed_slide, style = "histogram", plotType = "RunsPCT")
.save_pdf(p, "03_hist_runspct")

cat("  class distribution (MeanClass)...\n")
p <- plot_DistRuns(roh_bed_slide, style = "ChrBarPlot", plotType = "MeanClass")
.save_pdf(p, "04_dist_meanclass")

cat("  Froh distribution...\n")
p <- plot_Froh(roh_bed_slide)
.save_pdf(p, "05_froh")

cat("  Froh barplot...\n")
p <- plot_Froh(roh_bed_slide, plotType = "BarPlot")
.save_pdf(p, "06_froh_barplot")

cat("  Froh boxplot...\n")
p <- plot_Froh(roh_bed_slide, plotType = "BoxPlot")
.save_pdf(p, "07_froh_boxplot")

cat("  SNPs in runs Manhattan...\n")
grps <- unique(roh_bed_slide$runs$group)
for (g in grps) {
    p <- plot_manhattanRuns(roh_bed_slide, group = g)
    .save_pdf(p, paste0("08_manhattan_sliding - ", g))
}

cat("  Stacked runs per group...\n")
for (g in grps) {
    p <- plot_StackedRuns(roh_bed_slide, group = g)
    .save_pdf(p, paste0("09_stacked_runs_", g))
}

cat("  All-chromosome run plot...\n")
p <- plot_Runs(roh_bed_slide)
.save_pdf(p, "10_plot_runs_AllChromosomes", w = 14, h = 10)

# =============================================================================
# 10. Plots — ROHom consecutive
# =============================================================================
.section("10. Plots | ROHom | consecutive")

for (g in unique(roh_bed_cons$runs$group)) {
    p <- plot_manhattanRuns(roh_bed_cons, group = g)
    .save_pdf(p, paste0("11_manhattan_cons - ", g))
}

# =============================================================================
# 11. Plots — ROHet
# =============================================================================
.section("11. Plots | ROHet")

p <- plot_DistRuns(roh_bed_rohet, style = "violin", plotType = "sum")
.save_pdf(p, "12_rohet_violin_sum")

p <- plot_DistRuns(roh_bed_rohet, style = "ChrBarPlot", plotType = "MeanClass")
.save_pdf(p, "13_rohet_dist")

# =============================================================================
# 12. rohIslands
# =============================================================================
.section("12. rohIslands | sliding | p99 | n_perm=500")

islands_slide <- rohIslands(
    roh_bed_slide,
    n_perm     = 500,
    percentile = 0.99,
    seed       = 42
)
print(islands_slide)
cat("ROH islands:", sum(islands_slide$island), "\n")

.section("12b. rohIslands | consecutive | p95 | n_perm=500")

islands_cons <- rohIslands(
    roh_bed_cons,
    n_perm     = 500,
    percentile = 0.95,
    seed       = 42
)
cat("ROH islands (cons, p95):", sum(islands_cons$island), "\n")

cat("  ROH island Manhattan plot (sliding)...\n")
p <- plot(islands_slide)
.save_pdf(p, "14_islands_slide_p99", w = 14, h = 5)

p <- plot(islands_cons)
.save_pdf(p, "15_islands_cons_p95", w = 14, h = 5)

# =============================================================================
# 13. Froh by chromosome
# =============================================================================
.section("13. Froh_inbreeding chromosome-by-chromosome")

froh_chr <- Froh_inbreeding(roh_bed_slide, genome_wide = FALSE)
cat("Froh_chr rows:", nrow(froh_chr), "\n")
print(head(froh_chr, 5))

p <- plot_Froh(roh_bed_slide, genome_wide = FALSE, plotType = "BoxPlot")
.save_pdf(p, "16_froh_chr_boxplot")

p <- plot_Froh(roh_bed_slide, genome_wide = FALSE, plotType = "BarPlot")
.save_pdf(p, "17_froh_chr_barplot")

# =============================================================================
# 14. as_ROH and as.data.frame
# =============================================================================
.section("14. as_ROH / as.data.frame.ROH")

df <- as.data.frame(roh_bed_slide)
cat("as.data.frame rows:", nrow(df), "cols:", ncol(df), "\n")

roh_rebuilt <- as_ROH(roh_bed_slide$runs, roh_bed_slide$snp_map,
                      roh_bed_slide$sample_info)
cat("as_ROH rebuild OK:", inherits(roh_rebuilt, "ROH"), "\n")

# =============================================================================
# 15. reportRUNS — Markdown + HTML
# =============================================================================
.section("15. reportRUNS | markdown")

out_md <- reportRUNS(
    roh_bed_slide,
    islands    = islands_slide,
    format     = "markdown",
    output_dir = OUT_DIR,
    prefix     = "test_report_md",
    overwrite  = TRUE,
    verbose    = TRUE
)
cat("Markdown report:", out_md$report_file, "\n")
cat("Plots:", length(out_md$plots), "\n")

.section("15b. reportRUNS | html (self-contained)")

out_html <- reportRUNS(
    roh_bed_slide,
    islands    = islands_slide,
    format     = "html",
    output_dir = OUT_DIR,
    prefix     = "test_report_html",
    overwrite  = TRUE,
    verbose    = TRUE
)
cat("HTML report:", out_html$report_file, "\n")
cat("File size:", round(file.info(out_html$report_file)$size / 1024), "KB\n")

# =============================================================================
# 16. reportRUNS — ROHet
# =============================================================================
.section("16. reportRUNS | ROHet | html")

out_het <- reportRUNS(
    roh_bed_rohet,
    format     = "html",
    output_dir = OUT_DIR,
    prefix     = "test_report_rohet",
    overwrite  = TRUE,
    verbose    = TRUE
)
cat("ROHet HTML report:", out_het$report_file, "\n")

# =============================================================================
# 17. PED format report (if available)
# =============================================================================
if (!is.null(roh_ped_slide)) {
    .section("17. reportRUNS | PED scan | html")
    out_ped <- reportRUNS(
        roh_ped_slide,
        format     = "html",
        output_dir = OUT_DIR,
        prefix     = "test_report_ped",
        overwrite  = TRUE,
        verbose    = TRUE
    )
    cat("PED HTML report:", out_ped$report_file, "\n")
}

# =============================================================================
# Done
# =============================================================================
.section("DONE")
cat("All results saved in:", normalizePath(OUT_DIR), "\n")
cat("\nFiles generated:\n")
cat(paste(" ", list.files(OUT_DIR, recursive = FALSE)), sep = "\n")
