suppressPackageStartupMessages(library(detectRUNS))

BED   <- "Ext_Data/ADAPTmap_genotypeTOP_20161201.bed"
N_PER <- 20L   # quick probe; extrapolate to 100 below
N_THR <- parallel::detectCores(logical = FALSE)

rss_mb <- function() {
    as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
}

cat("=========================================\n")
cat(sprintf(" Dataset : ADAPTmap\n"))
cat(sprintf(" SNPs    : 53,347  |  Samples : 4,653\n"))
cat(sprintf(" Threads : %d physical cores\n", N_THR))
cat(sprintf(" Perms   : %d (probe; ×5 → 100 perm estimate)\n", N_PER))
cat("=========================================\n\n")

# ---- 1. scanRUNS -----------------------------------------------------------
cat("[1] scanRUNS ...\n")
t0  <- proc.time()
runs <- scanRUNS(
    genoFile     = BED,
    ROHet        = FALSE,
    method       = "sliding",
    minSNP       = 20,
    maxOpp       = 1,
    maxMiss      = 1,
    maxGap       = 10e6,
    windowSize   = 20,
    threshold    = 0.05,
    minLengthBps = 1e6,
    minDensity   = 1/1e3,
    nThreads     = N_THR,
    verbose      = FALSE
)
dt_scan <- (proc.time() - t0)[["elapsed"]]
rss_post_scan <- rss_mb()
cat(sprintf("   done: %.1f s  |  %d runs  |  RSS %.0f MB\n\n",
            dt_scan, nrow(runs$runs), rss_post_scan))

# ---- 2. runsIslands (20 perm probe) ----------------------------------------
cat(sprintf("[2] runsIslands (%d perms) ...\n", N_PER))
rss_pre <- rss_mb()
t1 <- proc.time()
islands <- runsIslands(
    roh      = runs,
    n_perm   = N_PER,
    nThreads = N_THR,
    verbose  = FALSE
)
dt_isl   <- (proc.time() - t1)[["elapsed"]]
rss_peak <- rss_mb()

cat(sprintf("   done: %.1f s  |  island SNPs: %d\n", dt_isl, sum(islands$is_island)))
cat(sprintf("   RSS before: %.0f MB  |  after: %.0f MB  |  delta: %.0f MB\n",
            rss_pre, rss_peak, rss_peak - rss_pre))
cat(sprintf("\n   *** Estimated time for 100 perms: ~%.0f s (%.1f min) ***\n",
            dt_isl * 5, dt_isl * 5 / 60))
cat("\nBENCH OK\n")
