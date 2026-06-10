
# TODO: benchmark tableRuns() (R, Stats.R) vs tableRunsCpp() (C++ wrapper) and
#       decide which to keep; remove the other in a future release.
#       Also remove slidingRUNS.run() and consecutiveRUNS.run() once usage drops.

# Deprecated aliases — old ROH-naming kept for backward compatibility

#' Deprecated: use \code{\link{saveRUNS}} instead
#' @param ... passed to \code{saveRUNS}
#' @export
#' @keywords internal
saveROH <- function(...) {
  .Deprecated("saveRUNS")
  saveRUNS(...)
}

#' Deprecated: use \code{\link{loadRUNS}} instead
#' @param ... passed to \code{loadRUNS}
#' @export
#' @keywords internal
loadROH <- function(...) {
  .Deprecated("loadRUNS")
  loadRUNS(...)
}

#' Deprecated: use \code{\link{as_RUNS}} instead
#' @param ... passed to \code{as_RUNS}
#' @export
#' @keywords internal
as_ROH <- function(...) {
  .Deprecated("as_RUNS")
  as_RUNS(...)
}

#' Deprecated: use \code{\link{runsIslands}} instead
#' @param ... passed to \code{runsIslands}
#' @export
#' @keywords internal
rohIslands <- function(...) {
  .Deprecated("runsIslands")
  runsIslands(...)
}

#' Deprecated: use \code{\link{scanRUNS}} instead
#'
#' @param genotypeFile Path to the PED genotype file (was first argument).
#' @param mapFile Path to the MAP file (was second argument).
#' @param windowSize Sliding-window size in SNPs.
#' @param threshold Bjelland coverage-ratio threshold.
#' @param minSNP Minimum SNPs in a run.
#' @param ROHet If \code{TRUE}, detect heterozygosity runs.
#' @param maxOppWindow Max opposite genotypes in a window (was \code{maxOppWindow}).
#' @param maxMissWindow Max missing genotypes in a window (was \code{maxMissWindow}).
#' @param maxGap Maximum gap between SNPs in bp.
#' @param minLengthBps Minimum run length in bp.
#' @param minDensity Minimum SNP density in SNPs/kbp.
#' @param maxOppRun Max opposite genotypes in the full run.
#' @param maxMissRun Max missing genotypes in the full run.
#' @return A \code{\link{RUNS}} object (was a plain data.frame).
#' @export
#' @keywords internal
slidingRUNS.run <- function(genotypeFile, mapFile,
                            windowSize   = 15,
                            threshold    = 0.05,
                            minSNP       = 3,
                            ROHet        = FALSE,
                            maxOppWindow = 1,
                            maxMissWindow = 1,
                            maxGap       = 10^6,
                            minLengthBps = 1000,
                            minDensity   = 1/1000,
                            maxOppRun    = NULL,
                            maxMissRun   = NULL) {
  .Deprecated("scanRUNS")
  scanRUNS(
    genoFile     = genotypeFile,
    mapFile      = mapFile,
    method       = "sliding",
    windowSize   = windowSize,
    threshold    = threshold,
    minSNP       = minSNP,
    ROHet        = ROHet,
    maxOpp       = maxOppWindow,
    maxMiss      = maxMissWindow,
    maxGap       = maxGap,
    minLengthBps = minLengthBps,
    minDensity   = minDensity,
    maxOppRun    = maxOppRun,
    maxMissRun   = maxMissRun
  )
}

#' Deprecated: use \code{\link{scanRUNS}} instead
#'
#' @param genotypeFile Path to the PED genotype file (was first argument).
#' @param mapFile Path to the MAP file (was second argument).
#' @param ROHet If \code{TRUE}, detect heterozygosity runs.
#' @param maxOppRun Max opposite genotypes in the full run.
#' @param maxMissRun Max missing genotypes in the full run.
#' @param minSNP Minimum SNPs in a run.
#' @param minLengthBps Minimum run length in bp.
#' @param maxGap Maximum gap between SNPs in bp.
#' @return A \code{\link{RUNS}} object (was a plain data.frame).
#' @export
#' @keywords internal
consecutiveRUNS.run <- function(genotypeFile, mapFile,
                                ROHet        = FALSE,
                                maxOppRun    = 0,
                                maxMissRun   = 0,
                                minSNP       = 15,
                                minLengthBps = 1000,
                                maxGap       = 10^6) {
  .Deprecated("scanRUNS")
  scanRUNS(
    genoFile     = genotypeFile,
    mapFile      = mapFile,
    method       = "consecutive",
    ROHet        = ROHet,
    maxOppRun    = maxOppRun,
    maxMissRun   = maxMissRun,
    minSNP       = minSNP,
    minLengthBps = minLengthBps,
    maxGap       = maxGap
  )
}

# Suppress R CMD check NOTEs for column names used in data.table / ggplot2 NSE
utils::globalVariables(c(
  "CLASS", "group", "chrom", "freq", "id",
  "lengthBps", "nSNP", "MB",
  "CHROMOSOME", "CHR_LENGTH",
  "POPULATION", "IND", "COUNT", "START", "END",
  "SNP_NAME", "POSITION", "BREED", "PERCENTAGE",
  "value", "variable", "Froh_genome",
  "n_ROH", "total_length_bp"
))

.onAttach <- function(libname, pkgname) {
  version = packageVersion("detectRUNS")
  packageStartupMessage(paste("Using detectRUNS", version))
}

.onUnload <- function (libpath) {
  library.dynam.unload("detectRUNS", libpath)
}
