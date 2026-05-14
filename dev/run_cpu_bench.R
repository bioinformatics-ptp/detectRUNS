## CPU scaling benchmark: scanRUNS timing only (no downstream analysis)
## SELMOL + Innovagen_HD | ROHom | sliding+consecutive | very_lenient/lenient/strict | 1/10/14 threads

suppressPackageStartupMessages(library(detectRUNS))

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0L)
    setwd(dirname(dirname(normalizePath(script_path))))
cat("Working dir:", getwd(), "\n\n")

DATASETS <- list(
    SELMOL       = "Ext_Data/SELMOL_codACGT.bed",
    Innovagen_HD = "Ext_Data/Innovagen_HD.bed"
)

PARAMS <- list(
    very_lenient = list(minSNP=5,  maxOpp=3, maxMiss=3, minLen=5e4, maxGap=2e6,   win=10, thr=0.10),
    lenient      = list(minSNP=10, maxOpp=2, maxMiss=2, minLen=1e5, maxGap=1.5e6, win=15, thr=0.05),
    strict       = list(minSNP=20, maxOpp=1, maxMiss=1, minLen=5e5, maxGap=1e6,   win=20, thr=0.05)
)

METHODS    <- c("sliding", "consecutive")
CPU_COUNTS <- c(1L, 10L, 14L)

results <- data.frame(
    dataset   = character(),
    method    = character(),
    params    = character(),
    n_threads = integer(),
    n_runs    = integer(),
    scan_s    = numeric(),
    stringsAsFactors = FALSE
)

for (dsname in names(DATASETS)) {
    bed <- DATASETS[[dsname]]
    if (!file.exists(bed)) { cat("SKIP:", dsname, "(file not found)\n"); next }
    cat(sprintf("=== %s ===\n", dsname))

    for (method in METHODS) {
        for (pname in names(PARAMS)) {
            p <- PARAMS[[pname]]
            for (nt in CPU_COUNTS) {
                cat(sprintf("  %-12s | %-14s | %2d thr  ... ", method, pname, nt))
                flush(stdout())

                t0 <- proc.time()["elapsed"]
                scan <- tryCatch({
                    if (method == "sliding") {
                        scanRUNS(bed, method = "sliding", ROHet = FALSE,
                                 minSNP = p$minSNP, maxOpp = p$maxOpp, maxMiss = p$maxMiss,
                                 minLengthBps = p$minLen, maxGap = p$maxGap,
                                 windowSize = p$win, threshold = p$thr,
                                 nThreads = nt, verbose = FALSE)
                    } else {
                        scanRUNS(bed, method = "consecutive", ROHet = FALSE,
                                 minSNP = p$minSNP, maxOpp = p$maxOpp, maxMiss = p$maxMiss,
                                 minLengthBps = p$minLen, maxGap = p$maxGap,
                                 nThreads = nt, verbose = FALSE)
                    }
                }, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
                dt <- round(proc.time()["elapsed"] - t0, 2)

                nr <- if (!is.null(scan)) nrow(scan$runs) else NA_integer_
                cat(sprintf("%9s runs  %6.2f s\n", format(nr, big.mark=","), dt))

                results <- rbind(results, data.frame(
                    dataset   = dsname,
                    method    = method,
                    params    = pname,
                    n_threads = nt,
                    n_runs    = nr,
                    scan_s    = dt,
                    stringsAsFactors = FALSE
                ))
                rm(scan); gc(verbose = FALSE)
            }
        }
    }
    cat("\n")
}

write.csv(results, "dev/bench_cpu.csv", row.names = FALSE)
cat("Saved dev/bench_cpu.csv\n\n")

# =============================================================================
# HTML report — pivot table: rows = dataset/method/params, cols = 1/10/14 threads
# =============================================================================
html_path <- "dev/bench_cpu.html"

# Build pivot: for each (dataset, method, params) row, show time for each CPU count
row_keys <- unique(results[, c("dataset","method","params")])
row_keys <- row_keys[order(row_keys$dataset, row_keys$method,
                           match(row_keys$params, names(PARAMS))), ]

cpu_labels <- paste0(CPU_COUNTS, "-CPU")

html_rows <- character(0)
prev_ds <- ""

for (i in seq_len(nrow(row_keys))) {
    ds  <- row_keys$dataset[i]
    mth <- row_keys$method[i]
    prm <- row_keys$params[i]

    # Dataset header row
    if (ds != prev_ds) {
        html_rows <- c(html_rows,
            sprintf('<tr class="ds-header"><td colspan="8"><b>%s</b></td></tr>', ds))
        prev_ds <- ds
    }

    sub <- results[results$dataset == ds & results$method == mth & results$params == prm, ]
    sub <- sub[order(sub$n_threads), ]

    nr <- if (nrow(sub) > 0 && !is.na(sub$n_runs[1])) format(sub$n_runs[1], big.mark=",") else "—"

    cells <- vapply(CPU_COUNTS, function(nt) {
        r <- sub[sub$n_threads == nt, ]
        if (nrow(r) == 0 || is.na(r$scan_s)) return("<td>—</td>")
        sprintf("<td>%.2f s</td>", r$scan_s)
    }, character(1))

    # Speedup 1→14
    t1  <- sub$scan_s[sub$n_threads == 1L]
    t14 <- sub$scan_s[sub$n_threads == 14L]
    spd <- if (length(t1) == 1 && length(t14) == 1 && !is.na(t1) && !is.na(t14) && t14 > 0)
               sprintf("%.1f×", t1 / t14) else "—"

    html_rows <- c(html_rows, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="runs">%s</td>%s<td class="spd">%s</td></tr>',
        mth, prm, nr, paste(cells, collapse=""), spd))
}

html <- sprintf('<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>detectRUNS CPU Benchmark</title>
<style>
  body { font-family: Arial, sans-serif; font-size: 13px; margin: 30px; background: #f9f9f9; }
  h1 { font-size: 18px; margin-bottom: 4px; }
  p.sub { color: #555; margin-top: 2px; margin-bottom: 20px; }
  table { border-collapse: collapse; background: #fff; box-shadow: 0 1px 4px rgba(0,0,0,.1); }
  th, td { border: 1px solid #ddd; padding: 7px 12px; text-align: left; }
  th { background: #2c3e50; color: #fff; font-size: 12px; }
  tr:hover td { background: #f0f4f8; }
  .ds-header td { background: #ecf0f1; font-weight: bold; font-size: 13px; padding: 6px 10px; }
  .runs { color: #27ae60; font-weight: bold; }
  .spd  { color: #e74c3c; font-weight: bold; }
  td:nth-child(4), td:nth-child(5), td:nth-child(6) { text-align: right; font-family: monospace; }
  .note { margin-top: 16px; font-size: 11px; color: #888; }
</style>
</head>
<body>
<h1>detectRUNS — CPU Scaling Benchmark (ROHom)</h1>
<p class="sub">Generated: %s &nbsp;|&nbsp; scanRUNS() wall-clock time only, no downstream analysis</p>
<table>
<thead>
  <tr>
    <th>Method</th><th>Params</th><th>n_runs</th>
    <th>1-CPU</th><th>10-CPU</th><th>14-CPU</th>
    <th>Speedup (1→14)</th>
  </tr>
</thead>
<tbody>
%s
</tbody>
</table>
<p class="note">
  Datasets: SELMOL (4,095 animals, 44K SNPs, 5 breeds) &nbsp;|&nbsp;
  Innovagen_HD (1,009 animals, 777K SNPs, 1 breed)<br>
  n_runs shown is from the 1-CPU run (identical across thread counts if algorithm is deterministic).
</p>
</body>
</html>',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(html_rows, collapse="\n"))

writeLines(html, html_path)
cat(sprintf("HTML report written: %s\n", html_path))
cat(sprintf("Done at %s\n", format(Sys.time(), "%H:%M:%S")))
