## Quick speed benchmark — BED path, 1 thread vs all threads
## Two parameter sets per dataset; reports wall-clock seconds and runs found.

library(detectRUNS)

EXT   <- "Ext_Data"
CORES <- parallel::detectCores()

datasets <- list(
  list(name = "SELMOL",    bed = file.path(EXT, "SELMOL_codACGT.bed")),
  list(name = "ADAPTmap",  bed = file.path(EXT, "ADAPTmap_genotypeTOP_20161201.bed")),
  list(name = "suini_12",  bed = file.path(EXT, "suini_12_plink.bed")),
  list(name = "Innovagen_HD", bed = file.path(EXT, "Innovagen_HD.bed"))
)

param_sets <- list(
  list(label = "default",  windowSize = 15L, minSNP = 15L, maxOpp = 1L, maxMiss = 1L,
       threshold = 0.05, minLengthBps = 500000, maxGap = 5000000),
  list(label = "relaxed",  windowSize = 10L, minSNP = 10L, maxOpp = 2L, maxMiss = 2L,
       threshold = 0.05, minLengthBps = 100000, maxGap = 10000000)
)

cat(sprintf("Machine: %d cores\n\n", CORES))
cat(sprintf("%-14s %-10s %6s %6s %8s\n", "Dataset", "Params", "1-CPU", "All-CPU", "N runs"))
cat(strrep("-", 52), "\n")

for (ds in datasets) {
  bim <- sub("\\.bed$", ".bim", ds$bed)
  fam <- sub("\\.bed$", ".fam", ds$bed)
  snps <- nrow(read.table(bim, header = FALSE))
  inds <- nrow(read.table(fam, header = FALSE))
  cat(sprintf("\n%s  (%d SNPs | %d ind)\n", ds$name, snps, inds))

  for (ps in param_sets) {
    # 1 thread
    t1 <- system.time(
      res <- scanRUNS(
        genoFile     = ds$bed,
        windowSize   = ps$windowSize,
        minSNP       = ps$minSNP,
        maxOpp       = ps$maxOpp,
        maxMiss      = ps$maxMiss,
        threshold    = ps$threshold,
        minLengthBps = ps$minLengthBps,
        maxGap       = ps$maxGap,
        nThreads     = 1L,
        verbose      = FALSE
      )
    )["elapsed"]
    nruns <- nrow(res$runs)

    # all threads
    tN <- system.time(
      scanRUNS(
        genoFile     = ds$bed,
        windowSize   = ps$windowSize,
        minSNP       = ps$minSNP,
        maxOpp       = ps$maxOpp,
        maxMiss      = ps$maxMiss,
        threshold    = ps$threshold,
        minLengthBps = ps$minLengthBps,
        maxGap       = ps$maxGap,
        nThreads     = CORES,
        verbose      = FALSE
      )
    )["elapsed"]

    cat(sprintf("  %-12s %-10s %5.1fs %5.1fs   %7d runs\n",
                "", ps$label, t1, tN, nruns))
  }
}

cat("\n", strrep("-", 52), "\n")
cat(sprintf("All-CPU = %d threads\n", CORES))
