##############################################################################
## detectRUNS — integration test script
##
## Run this script interactively or with: source("test_detectRUNS.R")
##
## Dataset: Kijas2016 sheep subset (100 individuals, 4841 SNPs, 2 breeds)
##          Available in both BED and PED format inside inst/extdata/
##############################################################################

VERBOSE_TESTS <- TRUE   # set FALSE to disable the wrapper and run raw

# --------------------------------------------------------------------------
# Helper: try_block(label, { ... code ... })
# --------------------------------------------------------------------------
try_block <- function(label, expr) {
  if (!VERBOSE_TESTS) return(invisible(force(expr)))
  cat(sprintf("\n[TEST] %s\n", label))
  result <- tryCatch(
    withCallingHandlers(
      { val <- expr; cat("  --> PASS\n"); val },
      warning = function(w) {
        cat(sprintf("  --> WARN: %s\n", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat(sprintf("  --> FAIL: %s\n", conditionMessage(e)))
      NULL
    }
  )
  invisible(result)
}

# --------------------------------------------------------------------------
# 0. Load package
# --------------------------------------------------------------------------
try_block("Load detectRUNS", {
  library(detectRUNS)
})

# --------------------------------------------------------------------------
# 1. File paths — sheep dataset only
# --------------------------------------------------------------------------
sheep_ped <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
sheep_map <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
sheep_bed <- system.file("extdata", "Kijas2016_Sheep_subset.bed", package = "detectRUNS")
sheep_bim <- system.file("extdata", "Kijas2016_Sheep_subset.bim", package = "detectRUNS")
sheep_fam <- system.file("extdata", "Kijas2016_Sheep_subset.fam", package = "detectRUNS")

cat("\n=== File check ===\n")
for (f in c(sheep_ped, sheep_map, sheep_bed, sheep_bim, sheep_fam)) {
  cat(sprintf("  %s  %s\n", if (file.exists(f)) "[OK]" else "[MISSING]", f))
}

##############################################################################
## SECTION 1 — scanRUNS (BED engine)
##############################################################################
cat("\n\n========== SECTION 1: scanRUNS — BED engine ==========\n")

# 1a. Explicit bed/bim/fam paths, consecutive
result_cons_bed <- try_block("BED consecutive — explicit bed/bim/fam", {
  scanRUNS(
    sheep_bed, sheep_bim, sheep_fam,
    method       = "consecutive",
    minSNP       = 15,
    maxOpp       = 0,
    maxMiss      = 0,
    minLengthBps = 100000,
    maxGap       = 1e6,
    ROHet        = FALSE,
    verbose      = TRUE
  )
})

# 1b. Single .bed path auto-detect, sliding
result_slid_bed <- try_block("BED sliding — auto-detect from .bed path", {
  scanRUNS(
    sheep_bed,
    method       = "sliding",
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 100000,
    maxGap       = 1e6,
    windowSize   = 15,
    threshold    = 0.05,
    ROHet        = FALSE,
    verbose      = TRUE
  )
})

# 1c. Check return structure
try_block("BED result structure", {
  stopifnot(is.list(result_cons_bed))
  stopifnot(all(c("runs", "summary", "snp_freq", "chrom_map") %in% names(result_cons_bed)))
  stopifnot(data.table::is.data.table(result_cons_bed$runs))
  stopifnot(data.table::is.data.table(result_cons_bed$summary))
  cat(sprintf("  runs: %d rows | summary: %d individuals\n",
              nrow(result_cons_bed$runs), nrow(result_cons_bed$summary)))
  cat(sprintf("  runs columns: %s\n", paste(names(result_cons_bed$runs), collapse=", ")))
  cat(sprintf("  summary columns: %s\n", paste(names(result_cons_bed$summary), collapse=", ")))
})

# 1d. Check group stratification
try_block("BED group stratification", {
  groups <- unique(result_cons_bed$runs$group)
  cat(sprintf("  Groups in runs: %s\n", paste(sort(groups), collapse=", ")))
  cat(sprintf("  Groups in summary: %s\n",
              paste(sort(unique(result_cons_bed$summary$group)), collapse=", ")))
  n_zero <- sum(result_cons_bed$summary$n_ROH == 0)
  cat(sprintf("  Individuals with 0 ROH in summary: %d\n", n_zero))
})

##############################################################################
## SECTION 2 — scanRUNS (PED engine)
##############################################################################
cat("\n\n========== SECTION 2: scanRUNS — PED engine ==========\n")

# 2a. Auto-detect from .ped path, consecutive
result_cons_ped <- try_block("PED consecutive — auto-detect from .ped path", {
  scanRUNS(
    sheep_ped,
    method       = "consecutive",
    minSNP       = 15,
    maxOpp       = 0,
    maxMiss      = 0,
    minLengthBps = 100000,
    maxGap       = 1e6,
    ROHet        = FALSE,
    verbose      = TRUE
  )
})

# 2b. Explicit mapFile, sliding
result_slid_ped <- try_block("PED sliding — explicit mapFile", {
  scanRUNS(
    sheep_ped,
    mapFile      = sheep_map,
    method       = "sliding",
    minSNP       = 15,
    maxOpp       = 1,
    maxMiss      = 1,
    minLengthBps = 100000,
    maxGap       = 1e6,
    windowSize   = 15,
    threshold    = 0.05,
    ROHet        = FALSE,
    verbose      = TRUE
  )
})

# 2c. PED result structure
try_block("PED result structure", {
  stopifnot(is.list(result_cons_ped))
  stopifnot(all(c("runs", "summary", "snp_freq", "chrom_map") %in% names(result_cons_ped)))
  stopifnot(is.null(result_cons_ped$snp_freq))
  stopifnot(is.null(result_cons_ped$chrom_map))
  cat(sprintf("  runs: %d rows | summary: %d individuals\n",
              nrow(result_cons_ped$runs), nrow(result_cons_ped$summary)))
})

# 2d. BED vs PED run count comparison (same dataset, same parameters)
try_block("BED vs PED consecutive — run count comparison", {
  n_bed <- nrow(result_cons_bed$runs)
  n_ped <- nrow(result_cons_ped$runs)
  cat(sprintf("  BED runs: %d  |  PED runs: %d\n", n_bed, n_ped))
  if (n_bed != n_ped)
    message("  NOTE: counts differ — expected if fam group labels differ")
})

##############################################################################
## SECTION 3 — Deprecated functions (should warn but still work)
##############################################################################
cat("\n\n========== SECTION 3: Deprecated functions ==========\n")

result_old_cons <- try_block("consecutiveRUNS.run — deprecated wrapper", {
  suppressWarnings(
    consecutiveRUNS.run(
      sheep_ped, sheep_map,
      minSNP       = 15,
      maxOppRun    = 0,
      maxMissRun   = 0,
      minLengthBps = 100000,
      maxGap       = 1e6,
      ROHet        = FALSE
    )
  )
})

result_old_slid <- try_block("slidingRUNS.run — deprecated wrapper", {
  suppressWarnings(
    slidingRUNS.run(
      sheep_ped, sheep_map,
      windowSize    = 15,
      threshold     = 0.05,
      minSNP        = 15,
      maxOppWindow  = 1,
      maxMissWindow = 1,
      minLengthBps  = 100000,
      maxGap        = 1e6,
      ROHet         = FALSE
    )
  )
})

try_block("Deprecated functions return plain data.frame", {
  stopifnot(is.data.frame(result_old_cons))
  stopifnot(!data.table::is.data.table(result_old_cons))
  cat(sprintf("  consecutiveRUNS.run: %d rows\n", nrow(result_old_cons)))
  cat(sprintf("  slidingRUNS.run:     %d rows\n", nrow(result_old_slid)))
})

##############################################################################
## SECTION 4 — Save / Load binary ROH
##############################################################################
cat("\n\n========== SECTION 4: saveROH / loadROH ==========\n")

roh_file <- tempfile(fileext = ".roh")

try_block("saveROH — write BED consecutive results", {
  saveROH(result_cons_bed, roh_file)
  cat(sprintf("  Written: %s (%.1f KB)\n", roh_file, file.size(roh_file)/1024))
})

result_loaded <- try_block("loadROH — read back and compare", {
  loaded <- loadROH(roh_file)
  stopifnot(data.table::is.data.table(loaded$runs))
  n_orig <- nrow(result_cons_bed$runs)
  n_load <- nrow(loaded$runs)
  cat(sprintf("  Original: %d rows  |  Loaded: %d rows\n", n_orig, n_load))
  stopifnot(n_orig == n_load)
  loaded
})

##############################################################################
## SECTION 5 — Stats functions
##############################################################################
cat("\n\n========== SECTION 5: Stats functions ==========\n")

runs_sheep <- result_cons_ped$runs   # data.table input

try_block("Froh_inbreeding — genome-wide (data.table input)", {
  froh <- Froh_inbreeding(runs = runs_sheep, mapFile = sheep_map, genome_wide = TRUE)
  cat(sprintf("  Froh rows: %d | columns: %s\n", nrow(froh), paste(names(froh), collapse=", ")))
  cat(sprintf("  Groups: %s\n", paste(sort(unique(froh$group)), collapse=", ")))
  print(head(froh[order(froh$Froh_genome, decreasing=TRUE), ], 3))
})

try_block("Froh_inbreeding — chromosome-wide", {
  froh_chr <- Froh_inbreeding(runs = runs_sheep, mapFile = sheep_map, genome_wide = FALSE)
  cat(sprintf("  Froh_chr rows: %d\n", nrow(froh_chr)))
})

try_block("Froh_inbreedingClass (data.table input)", {
  froh_class <- Froh_inbreedingClass(runs = runs_sheep, mapFile = sheep_map, Class = 2)
  cat(sprintf("  Froh_class rows: %d | columns: %s\n",
              nrow(froh_class), paste(names(froh_class), collapse=", ")))
})

try_block("summaryRuns (data.table input, snpInRuns=FALSE)", {
  summ <- summaryRuns(
    runs         = runs_sheep,
    mapFile      = sheep_map,
    genotypeFile = sheep_ped,
    Class        = 2,
    snpInRuns    = FALSE
  )
  cat(sprintf("  Summary list elements: %s\n", paste(names(summ), collapse=", ")))
  cat("  summary_ROH_count_chr:\n")
  print(summ$summary_ROH_count_chr)
})

try_block("summaryRuns (snpInRuns=TRUE)", {
  summ_snp <- summaryRuns(
    runs         = runs_sheep,
    mapFile      = sheep_map,
    genotypeFile = sheep_ped,
    Class        = 2,
    snpInRuns    = TRUE
  )
  cat(sprintf("  SNPinRun rows: %d\n", nrow(summ_snp$SNPinRun)))
})

try_block("tableRuns (data.table input)", {
  tbl <- tableRuns(
    runs         = runs_sheep,
    genotypeFile = sheep_ped,
    mapFile      = sheep_map,
    threshold    = 0.5
  )
  cat(sprintf("  tableRuns rows: %d | columns: %s\n",
              nrow(tbl), paste(names(tbl), collapse=", ")))
  print(head(tbl, 3))
})

##############################################################################
## SECTION 6 — Plot functions
##############################################################################
cat("\n\n========== SECTION 6: Plot functions ==========\n")

runs_chr1 <- runs_sheep[runs_sheep$chrom == "1", ]

try_block("plot_Runs (chr1 only)", {
  plot_Runs(runs = runs_chr1, suppressInds = TRUE, savePlots = FALSE)
})

try_block("plot_StackedRuns (chr1 only)", {
  plot_StackedRuns(runs = runs_chr1, savePlots = FALSE)
})

try_block("plot_ViolinRuns — sum", {
  plot_ViolinRuns(runs = runs_sheep, method = "sum", savePlots = FALSE)
})

try_block("plot_ViolinRuns — mean", {
  plot_ViolinRuns(runs = runs_sheep, method = "mean", savePlots = FALSE)
})

try_block("plot_PatternRuns — sum", {
  plot_PatternRuns(runs = runs_sheep, mapFile = sheep_map, method = "sum", savePlots = FALSE)
})

try_block("plot_PatternRuns — mean", {
  plot_PatternRuns(runs = runs_sheep, mapFile = sheep_map, method = "mean", savePlots = FALSE)
})

try_block("plot_DistributionRuns — All styles", {
  plot_DistributionRuns(runs = runs_sheep, mapFile = sheep_map,
                        style = "All", savePlots = FALSE, Class = 2)
})

try_block("plot_InbreedingChr — All styles", {
  plot_InbreedingChr(runs = runs_sheep, mapFile = sheep_map,
                     style = "All", savePlots = FALSE)
})

try_block("plot_SnpsInRuns (chr1 only)", {
  plot_SnpsInRuns(runs = runs_chr1, genotypeFile = sheep_ped,
                  mapFile = sheep_map, savePlots = FALSE)
})

try_block("plot_manhattanRuns", {
  plot_manhattanRuns(runs = runs_sheep, genotypeFile = sheep_ped,
                     mapFile = sheep_map, savePlots = FALSE,
                     plotTitle = "Sheep ROHom")
})

##############################################################################
## SECTION 7 — Edge cases
##############################################################################
cat("\n\n========== SECTION 7: Edge cases ==========\n")

try_block("scanRUNS — invalid file gives clear error", {
  tryCatch(
    scanRUNS("nonexistent_file.bed"),
    error = function(e) cat(sprintf("  Got expected error: %s\n", conditionMessage(e)))
  )
})

try_block("scanRUNS — non-string genoFile gives clear error", {
  tryCatch(
    scanRUNS(123),
    error = function(e) cat(sprintf("  Got expected error: %s\n", conditionMessage(e)))
  )
})

try_block("saveROH / loadROH round-trip — row count and column names match", {
  if (!is.null(result_loaded)) {
    orig_cols <- names(result_cons_bed$runs)
    load_cols <- names(result_loaded$runs)
    stopifnot(identical(orig_cols, load_cols))
    cat(sprintf("  Columns match: %s\n", paste(orig_cols, collapse=", ")))
  }
})

##############################################################################
## SECTION 8 — ROH object API
##############################################################################
cat("\n\n========== SECTION 8: ROH object API ==========\n")

try_block("print(res) works on ROH object", {
  print(result_cons_ped)
})

try_block("as.data.frame(res) returns the runs data.frame", {
  df <- as.data.frame(result_cons_ped)
  stopifnot(is.data.frame(df))
  stopifnot(all(c("group", "id", "chrom", "nSNP", "from", "to", "lengthBps") %in% names(df)))
  cat(sprintf("  as.data.frame rows: %d\n", nrow(df)))
})

try_block("Froh_inbreeding(res) works without mapFile", {
  froh_roh <- Froh_inbreeding(result_cons_ped)
  cat(sprintf("  Froh rows: %d | columns: %s\n", nrow(froh_roh), paste(names(froh_roh), collapse=", ")))
  stopifnot("Froh_genome" %in% names(froh_roh))
})

try_block("summaryRuns(res) works without mapFile/genotypeFile", {
  summ_roh <- summaryRuns(result_cons_ped, Class = 2, snpInRuns = FALSE)
  cat(sprintf("  Summary list elements: %s\n", paste(names(summ_roh), collapse=", ")))
  stopifnot("result_Froh_genome_wide" %in% names(summ_roh))
})

try_block("as_ROH() with readExternalRuns output works", {
  runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package = "detectRUNS")
  ext_runs <- readExternalRuns(runsFile, program = "detectRUNS")
  roh_obj  <- as_ROH(ext_runs, mapFile = sheep_map, genotypeFile = sheep_ped,
                     method = "sliding", type = "ROHom")
  stopifnot(inherits(roh_obj, "ROH"))
  cat(sprintf("  as_ROH: %d runs, method=%s, type=%s\n",
              nrow(roh_obj$runs), roh_obj$method, roh_obj$type))
  print(roh_obj)
})

try_block("Backward compat: Froh_inbreeding(res$runs, mapFile=mapFile) still works", {
  froh_bc <- Froh_inbreeding(result_cons_ped$runs, mapFile = sheep_map)
  cat(sprintf("  Backward compat Froh rows: %d\n", nrow(froh_bc)))
  stopifnot("Froh_genome" %in% names(froh_bc))
})

cat("\n\n========== ALL TESTS COMPLETE ==========\n\n")
