###############################################################################
## Deprecated public wrappers — removed from detectRUNS 1.0.0
##
## These functions are preserved here for reference only.
## Use scanRUNS() instead.
###############################################################################

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
    if (!is.na(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)))
        nCores <- min(nCores, 2L)
    nCores <- as.integer(nCores)

    if (!is.logical(ROHet) || length(ROHet) != 1L)
        stop(paste("Unknown ROHet value:", ROHet, ". Must be TRUE or FALSE."))

    if (!file.exists(genotypeFile))
        stop(paste("File not found:", genotypeFile))

    map_df <- detectRUNS:::readMapFile(mapFile)
    colnames(map_df) <- c("Chrom", "SNP", "bps")
    gaps <- diff(map_df$bps)

    detectRUNS:::.sliding_ped(
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
    if (!is.na(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)))
        nCores <- min(nCores, 2L)
    nCores <- as.integer(nCores)

    if (!is.logical(ROHet) || length(ROHet) != 1L)
        stop(paste("Unknown ROHet value:", ROHet, ". Must be TRUE or FALSE."))

    if (!file.exists(genotypeFile))
        stop(paste("File not found:", genotypeFile))

    map_df <- detectRUNS:::readMapFile(mapFile)
    colnames(map_df) <- c("Chrom", "SNP", "bps")

    detectRUNS:::.consecutive_ped(
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
