###########################################################
### Binary save / load for RUNS scan results
###########################################################


#' Save RUNS scan results to a compact binary file
#'
#' Serialises the \code{$runs} table from a \code{scanRUNS()} result to a
#' self-contained binary file (format ROHB v1).  The file is smaller and
#' faster to read than CSV and can be loaded back with \code{loadRUNS()}.
#'
#' @param result Named list returned by \code{scanRUNS()}.  Must contain
#'   at least \code{$runs} and \code{$chrom_map}.
#' @param path   Output file path (e.g. \code{"results.roh"}).
#'
#' @return Invisibly returns \code{path}.
#'
#' @seealso \code{\link{loadRUNS}}, \code{\link{scanRUNS}}
#'
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                 minLengthBps = 100000)
#' saveRUNS(roh, tempfile(fileext = ".roh"))
#' }
saveRUNS <- function(result, path) {
    if (!inherits(result, "RUNS"))
        stop("'result' must be a RUNS object from scanRUNS() or as_RUNS()")
    if (!is.character(path) || length(path) != 1L)
        stop("'path' must be a single file path string")
    saveRDS(result, path, compress = TRUE)
    invisible(path)
}


#' Load RUNS scan results from a binary file
#'
#' Reads a binary file written by \code{saveRUNS()} and returns a named list
#' whose \code{$runs} element is a \code{data.table} in the same format as
#' \code{scanRUNS()$runs} (columns: group, id, chrom, nSNP, from, to,
#' lengthBps).
#'
#' @param path Path to a \code{.roh} binary file.
#'
#' @return A named list with element \code{$runs} (a \code{data.table}).
#'
#' @seealso \code{\link{saveRUNS}}, \code{\link{scanRUNS}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata", "Kijas2016_Sheep_subset.bed",
#'                         package = "detectRUNS")
#' roh     <- scanRUNS(bedFile, method = "sliding", minSNP = 15,
#'                     minLengthBps = 100000)
#' rohFile <- tempfile(fileext = ".roh")
#' saveRUNS(roh, rohFile)
#' roh2 <- loadRUNS(rohFile)
#' nrow(roh2$runs)
#' }
loadRUNS <- function(path) {
    if (!file.exists(path))
        stop(paste("File not found:", path))
    # Try new RDS format (full RUNS object); fall back to legacy C binary format
    tryCatch(
        readRDS(path),
        error = function(e) {
            result      <- C_load_roh(path)
            result$runs <- data.table::as.data.table(result$runs)
            result
        }
    )
}
