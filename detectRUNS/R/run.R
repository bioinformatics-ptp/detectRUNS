###########################################################
### Compute Genomic Runs in R (homozygosity/heterozygosity)
###########################################################


# ---------------------------------------------------------------------------
# Private PED scanning helpers (called by scanRUNS and deprecated wrappers)
# ---------------------------------------------------------------------------

#' @keywords internal
.consecutive_ped <- function(ped_file, map_df, ROHet,
                              maxOppRun, maxMissRun,
                              minSNP, minLengthBps, maxGap,
                              nCores = 1L, verbose = FALSE) {

    # map_df must be a plain data.frame — no data.table, safe to fork
    map_df <- as.data.frame(map_df)
    lines  <- readLines(ped_file)
    n_snp  <- nrow(map_df)
    N      <- length(lines)

    .process_one <- function(oneLine) {
        geno <- as.character(strsplit(oneLine, " ")[[1]])
        if (length(geno) - 6 != n_snp * 2)
            stop("Number of markers differ in mapFile and genotype: are those the same dataset?")
        animal <- list(FID = geno[1], IID = geno[2])
        geno   <- pedConvertCpp(geno[7:length(geno)])
        consecutiveRunsCpp(
            geno, animal,
            mapFile             = map_df,
            ROHet               = ROHet,
            minSNP              = minSNP,
            maxOppositeGenotype = maxOppRun,
            maxMiss             = maxMissRun,
            minLengthBps        = minLengthBps,
            maxGap              = maxGap
        )
    }

    # mclapply (fork) is efficient here: map_df is already a plain data.frame,
    # no data.table pointers, copy-on-write — no serialisation overhead.
    # Falls back to progress-bar lapply on Windows (no fork support).
    # Parallel path: process in chunks so the progress bar advances between batches.
    if (verbose)
        pb <- utils::txtProgressBar(min = 0, max = N, style = 3,
                                    char = "#", width = 50)
    results <- vector("list", N)

    if (nCores > 1L && .Platform$OS.type == "unix") {
        chunk_size <- max(nCores * 4L, 20L)
        chunks     <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
        done       <- 0L
        for (ch in chunks) {
            results[ch] <- parallel::mclapply(lines[ch], .process_one,
                                              mc.cores = nCores)
            done <- done + length(ch)
            if (verbose) utils::setTxtProgressBar(pb, done)
        }
    } else {
        for (i in seq_len(N)) {
            results[[i]] <- .process_one(lines[[i]])
            if (verbose) utils::setTxtProgressBar(pb, i)
        }
    }
    if (verbose) { close(pb); cat("\n") }

    RUNs <- do.call(rbind, results)
    if (is.null(RUNs))
        RUNs <- data.frame(group = character(), id = character(),
                           chrom = character(), nSNP = integer(),
                           from  = integer(),   to   = integer(),
                           lengthBps = integer())
    row.names(RUNs) <- NULL
    RUNs
}


#' @keywords internal
.sliding_ped <- function(ped_file, map_df, gaps,
                          windowSize, threshold,
                          minSNP, ROHet,
                          maxOppWindow, maxMissWindow,
                          maxGap, minLengthBps, minDensity,
                          maxOppRun, maxMissRun,
                          nCores = 1L, verbose = FALSE) {

    parameters <- list(
        windowSize    = windowSize,
        threshold     = threshold,
        minSNP        = minSNP,
        ROHet         = ROHet,
        maxOppWindow  = maxOppWindow,
        maxMissWindow = maxMissWindow,
        maxGap        = maxGap,
        minLengthBps  = minLengthBps,
        minDensity    = minDensity,
        maxOppRun     = maxOppRun,
        maxMissRun    = maxMissRun
    )

    # map_df must be a plain data.frame — no data.table, safe to fork
    map_df <- as.data.frame(map_df)
    lines  <- readLines(ped_file)
    n_snp  <- nrow(map_df)
    N      <- length(lines)

    # Validate marker count on the first line before spawning parallel workers
    # so the error propagates directly rather than being wrapped by mclapply.
    first_geno <- as.character(strsplit(lines[1L], " ")[[1]])
    if (length(first_geno) - 6L != n_snp * 2L)
        stop("Number of markers differ in mapFile and genotype: are those the same dataset?")

    .process_one <- function(oneLine) {
        geno <- as.character(strsplit(oneLine, " ")[[1]])
        animal <- list(FID = geno[1], IID = geno[2])
        geno   <- pedConvertCpp(geno[7:length(geno)])
        slidingRuns(geno, animal, map_df, gaps, parameters)
    }

    # mclapply (fork) is efficient here: map_df is already a plain data.frame,
    # no data.table pointers, copy-on-write — no serialisation overhead.
    # Falls back to progress-bar lapply on Windows (no fork support).
    # Parallel path: process in chunks so the progress bar advances between batches.
    if (verbose)
        pb <- utils::txtProgressBar(min = 0, max = N, style = 3,
                                    char = "#", width = 50)
    results <- vector("list", N)

    if (nCores > 1L && .Platform$OS.type == "unix") {
        chunk_size <- max(nCores * 4L, 20L)
        chunks     <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
        done       <- 0L
        for (ch in chunks) {
            results[ch] <- parallel::mclapply(lines[ch], .process_one,
                                              mc.cores = nCores)
            done <- done + length(ch)
            if (verbose) utils::setTxtProgressBar(pb, done)
        }
    } else {
        for (i in seq_len(N)) {
            results[[i]] <- .process_one(lines[[i]])
            if (verbose) utils::setTxtProgressBar(pb, i)
        }
    }
    if (verbose) { close(pb); cat("\n") }

    RUNs <- do.call(rbind, results)
    if (is.null(RUNs))
        RUNs <- data.frame(group = character(), id = character(),
                           chrom = character(), nSNP = integer(),
                           from  = integer(),   to   = integer(),
                           lengthBps = integer())
    row.names(RUNs) <- NULL
    RUNs
}


# ---------------------------------------------------------------------------
# Deprecated public wrappers (kept for backward compatibility)
# ---------------------------------------------------------------------------

#' Main function to detect RUNS (ROHom/ROHet) using sliding windows (a la Plink)
#'
#' \strong{Deprecated.}
#'
#' This function is deprecated. Use \code{\link{scanRUNS}} instead.
#'
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param windowSize the size of sliding window (number of SNP loci) (default = 15)
#' @param threshold the threshold of overlapping windows of the same state
#' (homozygous/heterozygous) to call a SNP in a RUN (default = 0.05)
#' @param minSNP minimum n. of SNP in a RUN (default = 3)
#' @param ROHet should we look for ROHet or ROHom? (default = FALSE)
#' @param maxOppWindow max n. of homozygous/heterozygous SNP in the
#' sliding window (default = 1)
#' @param maxMissWindow max. n. of missing SNP in the sliding window (default = 1)
#' @param maxGap max distance between consecutive SNP to be still considered a
#' potential run (default = 10^6 bps)
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param minDensity minimum n. of SNP per kbps (defaults to 0.1 = 1 SNP every 10 kbps)
#' @param maxOppRun max n. of opposite genotype SNPs in the run (optional)
#' @param maxMissRun max n. of missing SNPs in the run (optional)
#' @param nCores number of cores to use (deprecated, ignored)
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity.
#' @export
#'
#' @import ggplot2
#' @import utils
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#' \dontrun{
#' runs <- slidingRUNS.run(genotypeFile, mapFile, windowSize = 15, threshold = 0.1,
#' minSNP = 15, ROHet = FALSE,  maxOppWindow = 1, maxMissWindow = 1, maxGap=10^6,
#' minLengthBps = 100000,  minDensity = 1/10000)
#' }
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.sliding.csv", package="detectRUNS")
#' colClasses <- c(rep("character", 3), rep("numeric", 4))
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses)
#'
slidingRUNS.run <- function(genotypeFile, mapFile,
                             windowSize = 15, threshold = 0.05,
                             minSNP = 3, ROHet = FALSE,
                             maxOppWindow = 1, maxMissWindow = 1,
                             maxGap = 10^6, minLengthBps = 1000, minDensity = 1/1000,
                             maxOppRun = NULL, maxMissRun = NULL,
                             nCores = NULL) {

    .Deprecated(
        new = "scanRUNS",
        msg = paste0("'slidingRUNS.run' is deprecated.\n",
                     "Use scanRUNS(genoFile, method = \"sliding\") instead.\n",
                     "Note: minDensity is not supported by scanRUNS (BED engine).")
    )

    if (is.null(nCores) || identical(nCores, 0L) || identical(nCores, 0))
        nCores <- max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)
    nCores <- as.integer(nCores)

    if (!is.logical(ROHet) || length(ROHet) != 1L)
        stop(paste("Unknown ROHet value:", ROHet, ". Must be TRUE or FALSE."))

    if (!file.exists(genotypeFile))
        stop(paste("File not found:", genotypeFile))

    map_df <- readMapFile(mapFile)
    colnames(map_df) <- c("Chrom", "SNP", "bps")
    gaps <- diff(map_df$bps)

    .sliding_ped(
        ped_file      = genotypeFile,
        map_df        = map_df,
        gaps          = gaps,
        windowSize    = windowSize,
        threshold     = threshold,
        minSNP        = minSNP,
        ROHet         = ROHet,
        maxOppWindow  = maxOppWindow,
        maxMissWindow = maxMissWindow,
        maxGap        = maxGap,
        minLengthBps  = minLengthBps,
        minDensity    = minDensity,
        maxOppRun     = maxOppRun,
        maxMissRun    = maxMissRun,
        nCores        = nCores
    )
}


#' Main function to detect genomic RUNS (ROHom/ROHet) using the consecutive method
#'
#' \strong{Deprecated.}
#'
#' This function is deprecated. Use \code{\link{scanRUNS}} instead.
#'
#' @param genotypeFile genotype (.ped) file path
#' @param mapFile map file (.map) file path
#' @param ROHet should we look for ROHet or ROHom? (default = FALSE)
#' @param maxOppRun max n. of opposite genotype SNPs in the run (default = 0)
#' @param maxMissRun max n. of missing SNPs in the run (default = 0)
#' @param minSNP minimum n. of SNP in a RUN (default = 15)
#' @param minLengthBps minimum length of run in bps (defaults to 1000 bps = 1 kbps)
#' @param maxGap max distance between consecutive SNP in a window to be still
#' considered a potential run (defaults to 10^6)
#' @param nCores number of cores to use (deprecated, ignored)
#'
#' @return A dataframe with RUNs of Homozygosity or Heterozygosity.
#' @export
#'
#' @import ggplot2
#' @import utils
#'
#' @examples
#' # getting map and ped paths
#' genotypeFile <- system.file("extdata", "Kijas2016_Sheep_subset.ped", package = "detectRUNS")
#' mapFile <- system.file("extdata", "Kijas2016_Sheep_subset.map", package = "detectRUNS")
#' \dontrun{
#' runs <- consecutiveRUNS.run(genotypeFile, mapFile, minSNP = 15, ROHet = FALSE,
#' maxOppRun = 0, maxMissRun = 0, maxGap=10^6, minLengthBps = 100000)
#' }
#' runsFile <- system.file("extdata", "Kijas2016_Sheep_subset.consecutive.csv", package="detectRUNS")
#' colClasses <- c(rep("character", 3), rep("numeric", 4))
#' runs <- read.csv2(runsFile, header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses)
#'
consecutiveRUNS.run <- function(genotypeFile, mapFile,
                                 ROHet = FALSE,
                                 maxOppRun = 0, maxMissRun = 0,
                                 minSNP = 15, minLengthBps = 1000,
                                 maxGap = 10^6,
                                 nCores = NULL) {

    .Deprecated(
        new = "scanRUNS",
        msg = paste0("'consecutiveRUNS.run' is deprecated.\n",
                     "Use scanRUNS(genoFile, method = \"consecutive\") instead.")
    )

    if (is.null(nCores) || identical(nCores, 0L) || identical(nCores, 0))
        nCores <- max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)
    nCores <- as.integer(nCores)

    if (!is.logical(ROHet) || length(ROHet) != 1L)
        stop(paste("Unknown ROHet value:", ROHet, ". Must be TRUE or FALSE."))

    if (!file.exists(genotypeFile))
        stop(paste("File not found:", genotypeFile))

    map_df <- readMapFile(mapFile)
    colnames(map_df) <- c("Chrom", "SNP", "bps")

    .consecutive_ped(
        ped_file     = genotypeFile,
        map_df       = map_df,
        ROHet        = ROHet,
        maxOppRun    = maxOppRun,
        maxMissRun   = maxMissRun,
        minSNP       = minSNP,
        minLengthBps = minLengthBps,
        maxGap       = maxGap,
        nCores       = nCores
    )
}
