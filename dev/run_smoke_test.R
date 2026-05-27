###############################################################################
## Smoke test — all Ext_Data datasets, 2 parameter sets
## BED path: 1 / 10 / N_CORES threads
## PED path: 1 thread (legacy R engine, no OpenMP)
## Goal: no crashes, sane run counts, BED vs PED consistency check
###############################################################################

suppressPackageStartupMessages(library(detectRUNS))

N_CORES <- parallel::detectCores(logical = FALSE)
.ts  <- function() format(Sys.time(), "[%H:%M:%S]")
.sec <- function(t) round(as.numeric(t["elapsed"]), 2)
.mem <- function() {
    tryCatch({
        info <- system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)
        round(as.numeric(trimws(info[1L])) / 1024, 1)
    }, error = function(e) NA_real_)
}

datasets <- list(
    list(name = "SELMOL",       bed = "Ext_Data/SELMOL_codACGT.bed"),
    list(name = "ADAPTmap",     bed = "Ext_Data/ADAPTmap_genotypeTOP_20161201.bed"),
    list(name = "suini_12",     bed = "Ext_Data/suini_12_plink.bed"),
    list(name = "Innovagen_HD", bed = "Ext_Data/Innovagen_HD.bed")
)

param_sets <- list(
    list(label = "default",
         windowSize = 15L, minSNP = 15L, maxOpp = 1L, maxMiss = 1L,
         threshold = 0.05, minLengthBps = 500000L, maxGap = 5000000L,
         minDensity = 1/1000),
    list(label = "strict",
         windowSize = 20L, minSNP = 20L, maxOpp = 1L, maxMiss = 1L,
         threshold = 0.05, minLengthBps = 1000000L, maxGap = 5000000L,
         minDensity = 1/50)
)

THREAD_COUNTS <- c(1L, 10L, N_CORES)

cat(sprintf("\n%s  SMOKE TEST — all Ext_Data datasets (max %d cores)\n", .ts(), N_CORES))
cat(strrep("=", 80), "\n")

# Collect rows for final summary table
summary_rows <- list()

.scan_once <- function(bed, method, ps, nthreads) {
    status <- "PASS"; err_msg <- ""
    args <- list(
        genoFile     = bed,
        method       = method,
        minSNP       = ps$minSNP,
        maxOpp       = ps$maxOpp,
        maxMiss      = ps$maxMiss,
        minLengthBps = ps$minLengthBps,
        maxGap       = ps$maxGap,
        minDensity   = ps$minDensity,
        nThreads     = nthreads,
        verbose      = FALSE
    )
    if (method == "sliding") {
        args$windowSize <- ps$windowSize
        args$threshold  <- ps$threshold
    }
    res <- NULL
    t <- tryCatch(
        system.time(res <- do.call(scanRUNS, args)),
        error = function(e) { status <<- "FAIL"; err_msg <<- conditionMessage(e); NULL }
    )
    list(status = status, err = err_msg, res = res, t = t)
}

for (ds in datasets) {
    bim  <- sub("\\.bed$", ".bim", ds$bed)
    snps <- nrow(read.table(bim, header = FALSE))
    inds <- nrow(read.table(sub("\\.bed$", ".fam", ds$bed), header = FALSE))

    cat(sprintf("\n%s  %s  (%d SNPs | %d ind)\n", .ts(), ds$name, snps, inds))
    cat(sprintf("  %-8s  %-8s  %10s  %8s  %8s  %8s  %8s\n",
                "params", "method", "n_runs", "1 CPU(s)", "10 CPU(s)",
                sprintf("%dCPU(s)", N_CORES), "mem(MB)"))
    cat(sprintf("  %s\n", strrep("-", 70)))

    for (ps in param_sets) {
        for (method in c("sliding", "consecutive")) {
            times <- numeric(length(THREAD_COUNTS))
            n_runs <- NA_integer_

            for (k in seq_along(THREAD_COUNTS)) {
                r <- .scan_once(ds$bed, method, ps, THREAD_COUNTS[k])
                if (r$status == "PASS") {
                    times[k] <- .sec(r$t)
                    if (k == length(THREAD_COUNTS)) {
                        n_runs <- nrow(r$res$runs)
                        mem_mb <- .mem()
                    }
                } else {
                    times[k] <- NA_real_
                    cat(sprintf("  ERROR [%s/%s/%d cpu]: %s\n",
                                ps$label, method, THREAD_COUNTS[k], r$err))
                }
                gc(verbose = FALSE)
            }

            mem_str <- if (!is.na(mem_mb)) sprintf("%.0f", mem_mb) else "N/A"
            cat(sprintf("  %-8s  %-8s  %10s  %8s  %8s  %8s  %8s\n",
                        ps$label, method,
                        format(n_runs, big.mark = ","),
                        ifelse(is.na(times[1]), "FAIL", sprintf("%.1f", times[1])),
                        ifelse(is.na(times[2]), "FAIL", sprintf("%.1f", times[2])),
                        ifelse(is.na(times[3]), "FAIL", sprintf("%.1f", times[3])),
                        mem_str))

            summary_rows[[length(summary_rows) + 1]] <- list(
                dataset   = ds$name,
                snps      = snps,
                animals   = inds,
                params    = ps$label,
                method    = method,
                n_runs    = n_runs,
                t_1cpu    = times[1],
                t_10cpu   = times[2],
                t_Ncpu    = times[3],
                peak_mem_mb = mem_mb
            )
        }
    }
}

cat(sprintf("\n%s  DONE\n", .ts()))
cat(sprintf("Machine: %d physical cores\n\n", N_CORES))

# Final CSV
out <- do.call(rbind, lapply(summary_rows, as.data.frame))
names(out)[names(out) == "t_Ncpu"] <- sprintf("t_%dcpu", N_CORES)
write.csv(out, "dev/smoke_test_results.csv", row.names = FALSE)
cat("Results saved to dev/smoke_test_results.csv\n")

# =============================================================================
# PED/MAP PATH — legacy R engine, single-threaded
# Note: _auto files contain autosomes only, so n_runs < BED (which includes
# sex chromosomes). We verify no crashes and report timing + run counts.
# =============================================================================
ped_datasets <- list(
    list(name = "SELMOL_auto",
         ped = "Ext_Data/SELMOL_codACGT_auto.ped",
         map = "Ext_Data/SELMOL_codACGT_auto.map"),
    list(name = "ADAPTmap_auto",
         ped = "Ext_Data/ADAPTmap_genotypeTOP_20161201_auto.ped",
         map = "Ext_Data/ADAPTmap_genotypeTOP_20161201_auto.map"),
    list(name = "suini_12",
         ped = "Ext_Data/suini_12_plink.ped",
         map = "Ext_Data/suini_12_plink.map")
)

cat(sprintf("\n%s  PED/MAP PATH — legacy R engine (1 thread)\n", .ts()))
cat(strrep("=", 80), "\n")
cat(sprintf("  %-14s  %-8s  %-12s  %10s  %8s  %8s\n",
            "dataset", "params", "method", "n_runs", "time(s)", "mem(MB)"))
cat(sprintf("  %s\n", strrep("-", 70)))

ped_rows <- list()
for (ds in ped_datasets) {
    if (!file.exists(ds$ped)) { cat("  SKIP:", ds$name, "(not found)\n"); next }
    snps <- nrow(read.table(ds$map, header = FALSE))
    inds <- nrow(read.table(ds$ped, header = FALSE))
    cat(sprintf("\n  %s  (%d SNPs | %d ind)\n", ds$name, snps, inds))

    for (ps in param_sets) {
        for (method in c("sliding", "consecutive")) {
            args <- list(
                genoFile     = ds$ped,
                mapFile      = ds$map,
                method       = method,
                minSNP       = ps$minSNP,
                maxOpp       = ps$maxOpp,
                maxMiss      = ps$maxMiss,
                minLengthBps = ps$minLengthBps,
                maxGap       = ps$maxGap,
                minDensity   = ps$minDensity,
                nThreads     = 1L,
                verbose      = FALSE
            )
            if (method == "sliding") {
                args$windowSize <- ps$windowSize
                args$threshold  <- ps$threshold
            }
            res <- NULL; status <- "PASS"
            t <- tryCatch(
                system.time(res <- do.call(scanRUNS, args)),
                error = function(e) { status <<- "FAIL"; cat("  ERROR:", conditionMessage(e), "\n"); NULL }
            )
            nr  <- if (!is.null(res)) nrow(res$runs) else NA_integer_
            sec <- if (!is.null(t))   round(t["elapsed"], 2) else NA_real_
            mem <- .mem()
            cat(sprintf("  %-14s  %-8s  %-12s  %10s  %8s  %8s\n",
                        ds$name, ps$label, method,
                        format(nr, big.mark = ","),
                        ifelse(is.na(sec), "FAIL", sprintf("%.1f", sec)),
                        ifelse(is.na(mem), "N/A",  sprintf("%.0f", mem))))
            ped_rows[[length(ped_rows) + 1]] <- list(
                dataset = ds$name, snps = snps, animals = inds,
                params = ps$label, method = method, path = "PED",
                n_runs = nr, t_1cpu = sec, peak_mem_mb = mem, status = status)
            gc(verbose = FALSE)
        }
    }
}

cat(sprintf("\n%s  PED section DONE\n", .ts()))

# Append PED rows to CSV
ped_out <- do.call(rbind, lapply(ped_rows, function(r) {
    data.frame(dataset=r$dataset, snps=r$snps, animals=r$animals,
               params=r$params, method=r$method, path=r$path,
               n_runs=r$n_runs, t_1cpu=r$t_1cpu, t_10cpu=NA, peak_mem_mb=r$peak_mem_mb,
               stringsAsFactors=FALSE)
}))
write.csv(ped_out, "dev/smoke_test_ped_results.csv", row.names = FALSE)
cat("PED results saved to dev/smoke_test_ped_results.csv\n")
